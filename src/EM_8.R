#EM_8.R

library(data.table)
library(optimx)

# 0. CONFIGURAZIONE DEI MECCANISMI

# Ogni meccanismo (es. ID = inverted deletion, DD = direct
# deletion) è descritto da 3 colonne del dataset:
#
# - mechanism_col : colonna logica (TRUE/FALSE) che dice se
#                   quella riga usa questo meccanismo
# - motif_pos_col : posizioni dei motivi (microomologie)
# - motif_len_col : lunghezza del motivo/pattern ripetuto
#
# Per aggiungere un nuovo meccanismo basta aggiungere un
# nuovo elemento a questa lista. NON bisogna modificare il
# resto dell'EM.

MECHANISMS <- list(

  ID = list(
    mechanism_col = "SD_inverted_deletion",
    motif_pos_col = "SD_ID_motif_pos",
    motif_len_col = "SD_ID_repeat_pat_len",
    len_scale = 0.5
  ),

  DD = list(
    mechanism_col = "SD_direct_deletion",
    motif_pos_col = "SD_DD_motif_pos",
    motif_len_col = "SD_DD_repeat_pat_len",
    len_scale = 0.5
  ),

  MMEJ = list(
    mechanism_col = "MMEJ",
    motif_pos_col = NULL,
    motif_len_col = "MMEJ_motif_len",
    len_scale = 1
  ),

  DI = list(
    mechanism_col = "SD_direct_insertion",
    motif_pos_col = "SD_DI_motif_pos",
    motif_len_col = "SD_DI_repeat_pat_len",
    len_scale = 0.5
  ),

  II = list(
    mechanism_col = "SD_inverted_insertion",
    motif_pos_col = "SD_II_motif_pos",
    motif_len_col = "SD_II_repeat_pat_len",
    len_scale = 0.5
  )

)

# 0b. CONFIGURAZIONE DEI BACKGROUND
#
# Componenti "di base" (nessun motivo, posizione uniforme).
# type = tipo di variante che la componente può spiegare:
#   "DEL" = delezione, "INS" = inserzione, "SUB" = sostituzione
#   (anche delins di lunghezza diversa);
#   NULL = qualsiasi tipo (comportamento "STANDARD": un solo NHEJ
#   generico che spiega tutto).
# Una componente non può spiegare varianti di tipo diverso dal
# suo: la sua likelihood è 0 e quindi pi stimato = 0 se nel
# file non ci sono varianti compatibili.

BACKGROUNDS <- list(
  NHEJ = list(type = "DEL"),
  RI   = list(type = "INS"),
  RS   = list(type = "SUB")
)

BACKGROUNDS_STANDARD <- list(NHEJ = list(type = NULL))


# 1. FUNZIONI BASE

# motif_strength()
# g = sum_i exp( beta*(motif_l_i - L0) - |distance_i|/lambda )
# (versione di riferimento; il codice usa calculate_g_all_fast)

motif_strength <- function(distance, motif_l, lambda, beta, L0 = 2) {
  stopifnot(length(distance) == length(motif_l), lambda > 0)
  sum(exp(beta * (motif_l - L0) - abs(distance) / lambda))
}

# position_probability()
# p = g_obs / (g_obs + sum(g_cf)); 0 se il totale è 0 o NA

position_probability <- function(g_obs, g_cf) {
  total <- g_obs + sum(g_cf)
  if (is.na(total) || total == 0) return(0)
  g_obs / total
}

# indel_length_prob(): lunghezza geometrica di parametro lambda_I
indel_length_prob <- function(L, lambda_I) dgeom(L - 1, prob = lambda_I)

# trim_variant()
# Toglie prefisso e suffisso comuni tra ANC e DER e restituisce
# le lunghezze di ciò che resta (a_len per ANC, d_len per DER).
# Funziona anche con basi di ancoraggio in stile VCF.

trim_variant <- function(anc, der) {
  a <- strsplit(toupper(anc), "")[[1]]
  d <- strsplit(toupper(der), "")[[1]]
  n <- min(length(a), length(d))

  p <- 0                                   # prefisso comune
  while (p < n && a[p + 1] == d[p + 1]) p <- p + 1

  s <- 0                                   # suffisso comune (senza sovrapporsi al prefisso)
  while (s < n - p && a[length(a) - s] == d[length(d) - s]) s <- s + 1

  list(a_len = length(a) - p - s, d_len = length(d) - p - s)
}

# get_variant_type()
# Dopo il trimming:
# - ANC vuoto, DER non vuoto   -> INS
# - DER vuoto, ANC non vuoto   -> DEL
# - entrambi non vuoti         -> SUB (anche di lunghezza diversa, es. 2 nt -> 3 nt)

get_variant_type <- function(anc, der) {
  t <- trim_variant(anc, der)
  if (t$a_len == 0 && t$d_len > 0) return("INS")
  if (t$d_len == 0 && t$a_len > 0) return("DEL")
  "SUB"
}

# compute_L()
# - stessa lunghezza        -> distanza di Hamming (>= 1)
# - delins di lunghezza diversa -> max tra i due resti dopo il trimming
# - indel puro              -> differenza di lunghezza

compute_L <- function(anc, der) {
  la <- nchar(anc); ld <- nchar(der)

  if (la == ld) {
    a_chars <- strsplit(anc, "")[[1]]
    d_chars <- strsplit(der, "")[[1]]
    L <- sum(a_chars != d_chars)
    return(if (L == 0) 1 else L)  # sicurezza: ANC identico a DER
  }

  t <- trim_variant(anc, der)
  if (t$a_len > 0 && t$d_len > 0)
    return(max(t$a_len, t$d_len))

  abs(la - ld)
}

max_L <- function(anc, der) {
  LA <- nchar(anc)
  LD <- nchar(der)
  if (LA >= LD) return(LA)
  else return(LD)
}

# extract_first_from_brackets
# "[True, False]" -> "True", "[3, 8]" -> "3". Se non ci sono
# quadre restituisce il valore invariato. NA, "" e "-" restano tali.

extract_first_from_brackets <- function(x) {
  if (is.na(x)) return(NA_character_)

  xc <- trimws(as.character(x))
  if (xc == "" || xc == "-") return(xc)

  if (startsWith(xc, "[") && endsWith(xc, "]")) {
    inner <- substr(xc, 2, nchar(xc) - 1)
    first <- strsplit(inner, ",")[[1]][1]
    first <- trimws(first)
    first <- gsub("^['\"]|['\"]$", "", first)
    return(first)
  }

  xc
}

# 2. PRE-PARSING DEI MOTIVI

parse_motif_positions <- function(x) {
  if (is.na(x)) return(numeric(0))

  val <- extract_first_from_brackets(x)
  if (is.na(val) || val == "" || val == "-") return(numeric(0))

  p <- suppressWarnings(as.numeric(val))
  if (is.na(p)) return(numeric(0))
  p
}


prepare_motif_data <- function(data, mechanism) {

  mechanism_col <- mechanism$mechanism_col
  motif_pos_col <- mechanism$motif_pos_col
  motif_len_col <- mechanism$motif_len_col
  len_scale <- if (is.null(mechanism$len_scale)) 1 else mechanism$len_scale
  has_position <- !is.null(motif_pos_col) && nzchar(motif_pos_col)

  # toglie eventuali [ ] e spezza sulla virgola
  split_commas <- function(x) {
    if (is.na(x)) return(character(0))
    x <- gsub("^\\s*\\[|\\]\\s*$", "", as.character(x))
    trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  }

  mech_raw <- as.character(data[[mechanism_col]])
  is_M <- vapply(mech_raw, function(x)
    isTRUE(any(toupper(split_commas(x)) == "TRUE")), logical(1))

  idx <- integer(0); distance <- numeric(0); motif_l <- numeric(0)

  for (i in which(is_M)) {

    flags <- toupper(split_commas(mech_raw[i])) == "TRUE"
    lens  <- suppressWarnings(as.numeric(split_commas(data[[motif_len_col]][i]))) * len_scale
    K <- length(flags)

    if (has_position) {
      groups <- trimws(strsplit(as.character(data[[motif_pos_col]][i]), ";", fixed = TRUE)[[1]])
      length(groups) <- K                      # completa con NA se mancano gruppi
    }

    for (k in seq_len(K)) {
      if (!flags[k]) next
      ml <- if (k <= length(lens)) lens[k] else NA_real_
      if (is.na(ml)) next

      if (has_position) {
        if (is.na(groups[k])) next
        d <- suppressWarnings(as.numeric(split_commas(groups[k])))
        d <- d[!is.na(d)]
      } else {
        d <- 0
      }
      if (length(d) == 0) next

      idx      <- c(idx, rep(i, length(d)))
      distance <- c(distance, d)
      motif_l  <- c(motif_l, rep(ml, length(d)))
    }
  }

  list(is_M = is_M, motif_row_idx = idx, motif_distance = distance, motif_length = motif_l)
}

# 3. CALCOLO VELOCE DI g
# w = exp(beta*(motif_length-L0) - |motif_distance|/lambda), sommato per riga.

calculate_g_all_fast <- function(motif_info, n_data, lambda, beta, L0 = 2) {

  g <- numeric(n_data)
  if (length(motif_info$motif_row_idx) == 0) return(g)

  w <- exp(beta * (motif_info$motif_length - L0) - abs(motif_info$motif_distance) / lambda)

  s <- rowsum(w, motif_info$motif_row_idx, reorder = FALSE)
  idx <- as.integer(rownames(s))
  g[idx] <- s[, 1]
  g
}

# 4. LIKELIHOOD

likelihood_mechanism <- function(L, lambda_I, g_obs, g_cf)
  indel_length_prob(L, lambda_I) * position_probability(g_obs, g_cf)

likelihood_NHEJ <- function(L, lambda_I_NHEJ, n_positions)
  indel_length_prob(L, lambda_I_NHEJ) * (1 / n_positions)

# 5. PRE-CALCOLO
#
# Prepara, prima dell'EM, la tabella "stats" (una riga per ogni
# allineamento candidato), motif_info e cf_index_list.
# Le likelihood di background NON sono calcolate qui: le ricalcola
# l'EM a ogni iterazione. Viene però salvato var_type (DEL/INS/SUB)
# che serve a decidere quali background possono spiegare la variante.

precompute_variant_stats_multi <- function(data, mechanisms, lambda_init, beta_init,
                                           lambda_I_init, L0 = 2) {

  required_cols <- c("variant_id", "alignment_ID", "ANC", "DER")
  for (mech_name in names(mechanisms)) {
    mech <- mechanisms[[mech_name]]
    required_cols <- c(required_cols, mech$mechanism_col, mech$motif_pos_col, mech$motif_len_col)
  }
  required_cols <- unique(required_cols)

  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0)
    stop(paste("Mancano le colonne:", paste(missing_cols, collapse = ", ")))

  data$.data_idx <- seq_len(nrow(data))

  motif_info <- list()
  for (mech_name in names(mechanisms))
    motif_info[[mech_name]] <- prepare_motif_data(data, mechanisms[[mech_name]])

  for (mech_name in names(mechanisms))
    data[[paste0("g_", mech_name)]] <- calculate_g_all_fast(
      motif_info[[mech_name]], nrow(data), lambda_init[[mech_name]], beta_init[[mech_name]], L0
    )

  variant_ids <- unique(data$variant_id)
  out_list <- vector("list", length(variant_ids))

  for (i in seq_along(variant_ids)) {
    id <- variant_ids[i]
    group <- data[variant_id == id]

    cf <- group[alignment_ID == 0]
    align_rows <- group[alignment_ID != 0]

    if (nrow(align_rows) == 0) {
      warning(paste("variant_id", id, "non ha alignment candidati."))
      next
    }

    first_mech <- names(mechanisms)[1]
    g_first <- group[[paste0("g_", first_mech)]]
    n_positions <- sum(!is.na(g_first))

    align_list <- vector("list", nrow(align_rows))

    for (j in seq_len(nrow(align_rows))) {
      row <- align_rows[j]
      L <- compute_L(row$ANC, row$DER)

      result_row <- data.frame(
        variant_id = id,
        alignment_ID = row$alignment_ID,
        indel_length = L,
        var_type = get_variant_type(row$ANC, row$DER),   # DEL / INS / SUB
        n_positions = n_positions,
        data_idx = row$.data_idx
      )

      for (mech_name in names(mechanisms)) {
        g_col <- paste0("g_", mech_name)
        g_obs <- row[[g_col]]
        g_cf <- cf[[g_col]]

        p_pos <- position_probability(g_obs, g_cf)
        LM <- likelihood_mechanism(L, lambda_I_init[[mech_name]], g_obs, g_cf)

        result_row[[paste0("g_", mech_name)]] <- g_obs
        result_row[[paste0("n_counterfactuals_", mech_name)]] <- length(g_cf)
        result_row[[paste0("p_position_", mech_name)]] <- p_pos
        result_row[[paste0("L_", mech_name)]] <- LM
      }

      align_list[[j]] <- result_row
    }

    out_list[[i]] <- do.call(rbind, align_list)
  }

  stats <- do.call(rbind, out_list)

  variant_levels <- unique(stats$variant_id)
  stats$variant_group <- match(stats$variant_id, variant_levels)

  cf_index_list <- lapply(variant_levels, function(id)
    data[variant_id == id & alignment_ID == 0, .data_idx])

  list(stats = stats, motif_info = motif_info, cf_index_list = cf_index_list, variant_levels = variant_levels)
}

# 6. OTTIMIZZAZIONE lambda E beta (M-STEP, un meccanismo)
#
# Q(lambda,beta) = sum_i P_old_i * log( p_position_i(lambda,beta) )
# Ottimizzazione vincolata (lambda>0, beta>=0) con L-BFGS-B.

optimize_lambda_beta_mechanism <- function(data, stats, P_old, motif_info, cf_index_list,
                                           lambda_init, beta_init, L0 = 2) {

  n_data <- nrow(data)
  obs_data_idx <- stats$data_idx
  variant_group <- stats$variant_group

  Q_function <- function(par) {
    lambda <- par[1]; beta <- par[2]

    if (!is.finite(lambda) || !is.finite(beta) || lambda <= 0 || beta < 0) return(1e100)

    g_all <- calculate_g_all_fast(motif_info, n_data, lambda, beta, L0)
    g_obs <- g_all[obs_data_idx]

    cf_sum_by_variant <- vapply(cf_index_list, function(idx)
      if (length(idx) == 0) 0 else sum(g_all[idx]), numeric(1))
    g_cf_sum <- cf_sum_by_variant[variant_group]

    total <- g_obs + g_cf_sum
    p_position <- numeric(nrow(stats))
    valid_total <- is.finite(g_obs) & is.finite(total) & total > 0
    p_position[valid_total] <- g_obs[valid_total] / total[valid_total]

    valid <- P_old > 0 & is.finite(P_old) & p_position > 0 & is.finite(p_position)
    if (!any(valid)) return(1e100)

    Q <- sum(P_old[valid] * log(p_position[valid]))
    -Q
  }

  opt <- optimx(
    par = c(lambda = lambda_init, beta = beta_init),
    fn = Q_function,
    method = "L-BFGS-B",
    lower = c(lambda = 1e-6, beta = 0),
    upper = c(lambda = 1e4, beta = 10),
    control = list(trace = 0)
  )

  best <- opt[1, ]

  list(lambda = best$lambda, beta = best$beta, Q = -best$value,
       convergence = best$convcode, optimx_result = opt)
}

# 7. EM MULTI-MECCANISMO
#
# Ogni allineamento osservato può essere generato da uno dei
# meccanismi con motivo (pi_<mech>) oppure da una delle componenti
# di background (pi_<bg>, ognuna con il suo lambda_I).
# Vincolo: sum(pi) + sum(pi_bg) = 1.
# Un background spiega solo varianti del proprio tipo (maschera
# `compat`): per le altre la likelihood è 0.

em_multi_mechanism <- function(data, stats, motif_info, cf_index_list, mechanisms, backgrounds,
                               pi_init, pi_bg_init, lambda_init, beta_init, lambda_I_init,
                               lambda_I_bg_init, max_iter = 100, tol = 1e-5,
                               L0 = 2, verbose = TRUE) {

  mechanism_names <- names(mechanisms)
  bg_names <- names(backgrounds)
  n_mechanisms <- length(mechanism_names)

  if (length(pi_init) != n_mechanisms) stop("Numero di pi_init non compatibile con i meccanismi.")
  if (length(lambda_init) != n_mechanisms) stop("Numero di lambda_init non compatibile.")
  if (length(beta_init) != n_mechanisms) stop("Numero di beta_init non compatibile.")
  if (length(lambda_I_init) != n_mechanisms) stop("Numero di lambda_I_init non compatibile.")
  if (length(pi_bg_init) != length(bg_names)) stop("Numero di pi_bg_init non compatibile.")
  if (abs(sum(pi_init) + sum(pi_bg_init) - 1) > 1e-8) stop("pi meccanismi + pi background deve fare 1.")

  # ---- Parametri correnti ----
  pi <- pi_init;           names(pi) <- mechanism_names
  pi_bg <- pi_bg_init;     names(pi_bg) <- bg_names
  lambda <- lambda_init
  beta <- beta_init
  lambda_I <- lambda_I_init
  lambda_I_bg <- lambda_I_bg_init; names(lambda_I_bg) <- bg_names

  # ---- Quantità comuni ----
  p_position_bg <- 1 / stats$n_positions   # posizione uniforme per tutti i background
  L_vec <- stats$indel_length
  v_id <- stats$variant_id
  first_of_variant <- !duplicated(v_id)

  # Maschera di compatibilità tipo variante <-> componente background
  compat <- list()
  for (b in bg_names) {
    tp <- backgrounds[[b]]$type
    compat[[b]] <- if (is.null(tp)) rep(TRUE, nrow(stats)) else stats$var_type == tp
  }

  history <- data.frame()
  eps <- 1e-6

  # Stima geometrica pesata, con guardia contro 0/0
  # (componente vuota: lambda_I resta al valore precedente)
  update_lambda_I <- function(P, old) {
    den <- sum(P * L_vec, na.rm = TRUE)
    if (!is.finite(den) || den <= 0) return(old)
    min(max(sum(P, na.rm = TRUE) / den, eps), 1 - eps)
  }

  for (iter in seq_len(max_iter)) {

    # ---------------- E-STEP ----------------
    g_all <- list()
    for (m in mechanism_names)
      g_all[[m]] <- calculate_g_all_fast(motif_info[[m]], nrow(data), lambda[[m]], beta[[m]], L0)

    L_mechanisms <- list()
    for (m in mechanism_names) {
      g_obs <- g_all[[m]][stats$data_idx]

      cf_sum_by_variant <- vapply(cf_index_list, function(idx)
        if (length(idx) == 0) 0 else sum(g_all[[m]][idx]), numeric(1))
      g_cf_sum <- cf_sum_by_variant[stats$variant_group]

      total_position <- g_obs + g_cf_sum
      p_position <- numeric(nrow(stats))
      ok <- is.finite(g_obs) & is.finite(total_position) & total_position > 0
      p_position[ok] <- g_obs[ok] / total_position[ok]

      stats[[paste0("p_position_", m)]] <- p_position
      stats[[paste0("g_", m)]] <- g_obs

      L_mechanisms[[m]] <- dgeom(L_vec - 1, prob = lambda_I[[m]]) * p_position
      stats[[paste0("L_", m)]] <- L_mechanisms[[m]]
    }

    # Likelihood background: 0 se il tipo di variante non è compatibile
    L_bg <- list()
    for (b in bg_names) {
      L_bg[[b]] <- as.numeric(compat[[b]]) *
        dgeom(L_vec - 1, prob = lambda_I_bg[[b]]) * p_position_bg
      stats[[paste0("L_", b)]] <- L_bg[[b]]
    }

    # ---------- RESPONSIBILITIES ----------
    num_mech <- list(); num_bg <- list()
    for (m in mechanism_names) num_mech[[m]] <- pi[[m]] * L_mechanisms[[m]]
    for (b in bg_names)        num_bg[[b]]   <- pi_bg[[b]] * L_bg[[b]]

    variant_contribution <- Reduce(`+`, c(num_mech, num_bg))
    variant_tot <- ave(variant_contribution, v_id, FUN = sum)
    bad <- is.na(variant_tot) | variant_tot == 0

    P <- list(); P_bg <- list()
    for (m in mechanism_names) {
      P[[m]] <- ifelse(bad, NA_real_, num_mech[[m]] / variant_tot)
      stats[[paste0("P_", m)]] <- P[[m]]
    }
    for (b in bg_names) {
      P_bg[[b]] <- ifelse(bad, NA_real_, num_bg[[b]] / variant_tot)
      stats[[paste0("P_", b)]] <- P_bg[[b]]
    }

    Loglik <- sum(log(variant_tot[first_of_variant]), na.rm = TRUE)

    # ---------------- M-STEP ----------------
    # 1. lambda e beta di ogni meccanismo
    lambda_new <- lambda; beta_new <- beta
    for (m in mechanism_names) {
      opt <- optimize_lambda_beta_mechanism(
        data = data, stats = stats, P_old = P[[m]], motif_info = motif_info[[m]],
        cf_index_list = cf_index_list, lambda_init = lambda[[m]], beta_init = beta[[m]], L0 = L0)
      lambda_new[[m]] <- opt$lambda
      beta_new[[m]] <- opt$beta
    }

    # 2. prior pi (meccanismi e background)
    pi_new <- pi
    for (m in mechanism_names) {
      Pv <- ave(P[[m]], v_id, FUN = sum)
      pi_new[[m]] <- mean(Pv[first_of_variant], na.rm = TRUE)
    }
    pi_bg_new <- pi_bg
    for (b in bg_names) {
      Pv <- ave(P_bg[[b]], v_id, FUN = sum)
      pi_bg_new[[b]] <- mean(Pv[first_of_variant], na.rm = TRUE)
    }

    # 3. lambda_I (geometrica pesata, forma chiusa)
    lambda_I_new <- lambda_I
    for (m in mechanism_names) lambda_I_new[[m]] <- update_lambda_I(P[[m]], lambda_I[[m]])
    lambda_I_bg_new <- lambda_I_bg
    for (b in bg_names) lambda_I_bg_new[[b]] <- update_lambda_I(P_bg[[b]], lambda_I_bg[[b]])

    # ---------------- CONVERGENZA ----------------
    delta_values <- c(abs(pi_new - pi), abs(pi_bg_new - pi_bg),
                      abs(lambda_I_bg_new - lambda_I_bg))
    for (m in mechanism_names)
      delta_values <- c(delta_values,
                        abs(lambda_new[[m]] - lambda[[m]]),
                        abs(beta_new[[m]] - beta[[m]]),
                        abs(lambda_I_new[[m]] - lambda_I[[m]]))
    delta <- max(delta_values, na.rm = TRUE)

    # ---------------- HISTORY ----------------
    history_row <- data.frame(iter = iter, Loglik = Loglik)
    for (m in mechanism_names) history_row[[paste0("pi_", m)]] <- pi_new[[m]]
    for (b in bg_names)        history_row[[paste0("pi_", b)]] <- pi_bg_new[[b]]
    for (m in mechanism_names) {
      history_row[[paste0("lambda_", m)]]   <- lambda_new[[m]]
      history_row[[paste0("beta_", m)]]     <- beta_new[[m]]
      history_row[[paste0("lambda_I_", m)]] <- lambda_I_new[[m]]
    }
    for (b in bg_names) history_row[[paste0("lambda_I_", b)]] <- lambda_I_bg_new[[b]]
    history_row$delta <- delta
    history <- rbind(history, history_row)

    # ---------------- STAMPA ----------------
    if (verbose) {
      cat(sprintf("Iter %3d:", iter))
      for (m in mechanism_names)
        cat(sprintf(" pi_%s=%.6f lambda_%s=%.6f beta_%s=%.6f lambda_I_%s=%.6f",
                    m, pi_new[[m]], m, lambda_new[[m]], m, beta_new[[m]], m, lambda_I_new[[m]]))
      for (b in bg_names)
        cat(sprintf(" pi_%s=%.6f lambda_I_%s=%.6f", b, pi_bg_new[[b]], b, lambda_I_bg_new[[b]]))
      cat(sprintf(" delta=%.2e Loglik=%.6f\n", delta, Loglik))
    }

    # ---------------- AGGIORNAMENTO PARAMETRI ----------------
    pi <- pi_new; pi_bg <- pi_bg_new
    lambda <- lambda_new; beta <- beta_new
    lambda_I <- lambda_I_new; lambda_I_bg <- lambda_I_bg_new

    if (delta < tol) {
      if (verbose) cat("tolleranza raggiunta in iter =", iter, "\n")
      break
    }
  }

  list(results = stats, pi = pi, pi_bg = pi_bg, lambda = lambda, beta = beta,
       lambda_I = lambda_I, lambda_I_bg = lambda_I_bg, history = history,
       n_iter = iter, converged = delta < tol)
}

# 8. FUNZIONE PRINCIPALE
#
# `mechanisms`: sottoinsieme di MECHANISMS
# `backgrounds`: sottoinsieme di BACKGROUNDS (default: NHEJ generico "legacy")

run_em_multi <- function(data, mechanisms, backgrounds = BACKGROUNDS_STANDARD,
                         lambda_init = NULL, beta_init = NULL, pi_init = NULL, pi_bg_init = NULL,
                         lambda_I_init = NULL, lambda_I_bg_init = NULL,
                         max_iter = 100, tol = 1e-8, L0 = 2, verbose = TRUE) {

  mechanism_names <- names(mechanisms)
  bg_names <- names(backgrounds)
  n_mech <- length(mechanism_names); n_bg <- length(bg_names)
  if (n_mech == 0) stop("Nessun meccanismo selezionato.")
  if (n_bg == 0) stop("Nessuna componente di background selezionata.")

  if (is.null(lambda_init))   { lambda_init   <- rep(50,  n_mech); names(lambda_init)   <- mechanism_names }
  if (is.null(beta_init))     { beta_init     <- rep(0.7, n_mech); names(beta_init)     <- mechanism_names }
  if (is.null(lambda_I_init)) { lambda_I_init <- rep(0.3, n_mech); names(lambda_I_init) <- mechanism_names }
  if (is.null(lambda_I_bg_init)) { lambda_I_bg_init <- rep(0.3, n_bg); names(lambda_I_bg_init) <- bg_names }

  # Prior iniziali: 0.5 ai meccanismi, il resto ai background
  if (is.null(pi_init)) { pi_init <- rep(0.5 / n_mech, n_mech); names(pi_init) <- mechanism_names }
  if (sum(pi_init) >= 1) stop("La somma dei prior dei meccanismi deve essere < 1.")

  if (is.null(pi_bg_init)) pi_bg_init <- rep((1 - sum(pi_init)) / n_bg, n_bg)
  names(pi_bg_init) <- bg_names

  cat("Pre-calcolo delle quantità fisse...\n")
  pre <- precompute_variant_stats_multi(
    data = data, mechanisms = mechanisms, lambda_init = lambda_init,
    beta_init = beta_init, lambda_I_init = lambda_I_init, L0 = L0)
  stats <- pre$stats
  cat("Pre-calcolo completato:", nrow(stats), "righe di alignment.\n\n")

  cat("Esecuzione EM su:", mechanism_names, "+", bg_names, "\n\n")

  em_multi_mechanism(
    data = data, stats = stats, motif_info = pre$motif_info, cf_index_list = pre$cf_index_list,
    mechanisms = mechanisms, backgrounds = backgrounds,
    pi_init = pi_init, pi_bg_init = pi_bg_init, lambda_init = lambda_init, beta_init = beta_init,
    lambda_I_init = lambda_I_init, lambda_I_bg_init = lambda_I_bg_init,
    max_iter = max_iter, tol = tol, L0 = L0, verbose = verbose)
}

# 9. CLI
#
# Esempi:
#   Solo NHEJ generico (come prima):
#     Rscript EM_7.R -i in.tsv -o out.tsv --mechanisms DI
#   Background separati per tipo:
#     Rscript EM_7.R -i in.tsv -o out.tsv --mechanisms DI --background NHEJ,RI,RS
#   Tutto insieme sui dati reali:
#     Rscript EM_7.R -i in.tsv -o out.tsv --mechanisms ID,DD,MMEJ,DI,II --background NHEJ,RI,RS

if (sys.nframe() == 0) {

  args <- commandArgs(trailingOnly = TRUE)

  get_arg <- function(flag, default = NULL, numeric = FALSE) {
    pos <- which(args == flag)
    if (length(pos) == 0) return(default)
    if (pos[1] >= length(args)) return(default)
    val <- args[pos[1] + 1]
    if (numeric) val <- as.numeric(val)
    val
  }

  input_file <- get_arg("-i")
  output_file <- get_arg("-o")
  mechanisms_arg <- get_arg("--mechanisms")
  background_arg <- get_arg("--background")   # es. NHEJ,RI,RS
  pi_arg <- get_arg("-p")
  tol <- get_arg("-t", 1e-8, TRUE)
  max_iter <- get_arg("-m", 100, TRUE)

  if (is.null(input_file) || is.null(output_file))
    stop("Uso:\nRscript EM_7.R -i input.tsv -o output.tsv --mechanisms ID,DD [--background NHEJ,RI,RS] [-p pi1,pi2] [-t tol] [-m max_iter]")

  if (is.null(mechanisms_arg)) stop("Devi specificare --mechanisms.")

  # ---- Meccanismi ----
  selected_names <- trimws(unlist(strsplit(mechanisms_arg, ",")))

  unknown <- setdiff(selected_names, names(MECHANISMS))
  if (length(unknown) > 0)
    stop(paste("Meccanismi non riconosciuti:", paste(unknown, collapse = ", "),
               "\nMeccanismi disponibili:", paste(names(MECHANISMS), collapse = ", ")))

  mechanisms <- MECHANISMS[selected_names]

  # ---- Background ----
  if (is.null(background_arg)) {
    backgrounds <- BACKGROUNDS_STANDARD           # un solo NHEJ generico, come prima
  } else {
    bg_names_sel <- trimws(unlist(strsplit(background_arg, ",")))
    unknown_bg <- setdiff(bg_names_sel, names(BACKGROUNDS))
    if (length(unknown_bg) > 0)
      stop(paste("Background non riconosciuti:", paste(unknown_bg, collapse = ", "),
                 "\nDisponibili:", paste(names(BACKGROUNDS), collapse = ", ")))
    backgrounds <- BACKGROUNDS[bg_names_sel]
  }

  cat("Input :", input_file, "\n")
  cat("Output:", output_file, "\n")
  cat("Meccanismi:", paste(selected_names, collapse = ", "),
      "+", paste(names(backgrounds), collapse = ", "), "\n\n")

  # ---- Importazione dati ----
  cat("Importazione del TSV...\n")
  data <- fread(input_file, sep = "\t", header = TRUE)
  cat("Importazione completata:", nrow(data), "righe x", ncol(data), "colonne\n\n")

  # ---- Parametri iniziali ----
  n_mechanisms <- length(selected_names)

  lambda_init <- rep(50, n_mechanisms); names(lambda_init) <- selected_names
  beta_init <- rep(0.7, n_mechanisms); names(beta_init) <- selected_names
  lambda_I_init <- rep(0.3, n_mechanisms); names(lambda_I_init) <- selected_names

  # Prior iniziali dei meccanismi (-p 0.25,0.25); il resto va ai background
  if (!is.null(pi_arg)) {
    pi_init <- as.numeric(unlist(strsplit(pi_arg, ",")))
    if (length(pi_init) != n_mechanisms)
      stop(paste("Con -p devi specificare", n_mechanisms, "valori."))
    names(pi_init) <- selected_names
  } else {
    pi_init <- rep(0.5 / n_mechanisms, n_mechanisms)
    names(pi_init) <- selected_names
  }

  # ---- Esecuzione EM ----
  em_out <- run_em_multi(
    data = data, mechanisms = mechanisms, backgrounds = backgrounds,
    lambda_init = lambda_init, beta_init = beta_init, pi_init = pi_init,
    lambda_I_init = lambda_I_init, max_iter = max_iter, tol = tol, L0 = 2, verbose = TRUE
  )

  # ---- Stampa risultati finali ----
  cat("\n==============================\nPRIOR FINALI\n==============================\n")
  for (m in selected_names)      cat(m, "=", em_out$pi[[m]], "\n")
  for (b in names(backgrounds))  cat(b, "=", em_out$pi_bg[[b]], "\n")

  cat("\n==============================\nPARAMETRI FINALI\n==============================\n")
  for (m in selected_names) {
    cat("\n", m, "\n")
    cat("lambda =", em_out$lambda[[m]], "\n")
    cat("beta =", em_out$beta[[m]], "\n")
    cat("lambda_I =", em_out$lambda_I[[m]], "\n")
  }
  for (b in names(backgrounds)) {
    cat("\n", b, "\n")
    cat("lambda_I =", em_out$lambda_I_bg[[b]], "\n")
  }

  cat("\nIterazioni:", em_out$n_iter, "\n")
  cat("Convergenza:", em_out$converged, "\n\n")

  # ---- Scrittura risultati ----
  cat("Scrittura dei risultati...\n")
  fwrite(em_out$results, file = output_file, sep = "\t", quote = FALSE, na = "NA")
  cat("Output salvato in:", output_file, "\n")

  # ---- Riepilogo ultima iterazione ----
  last_iter_file <- paste0(output_file, "_paramaters_estimated")
  last_iter <- tail(em_out$history, 1)

  cat("Scrittura della history...\n")
  fwrite(last_iter, file = last_iter_file, sep = "\t", quote = FALSE, na = "NA")
  cat("last iter salvata in:", last_iter_file, "\n")
}

