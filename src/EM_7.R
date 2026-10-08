#EM_7.R

library(data.table)
library(optimx)

# 0. CONFIGURAZIONE DEI MECCANISMI

# Ogni meccanismo (es. ID = inverted deletion, DD = direct
# deletion) è descritto da 3 colonne del dataset:
#
# - mechanism_col : colonna logica (TRUE/FALSE) che dice se
#                   quella riga usa questo meccanismo
# - motif_pos_col : posizioni dei motivi (microomologie)
#                   trovati, come stringa "1,4,7"
# - motif_len_col : lunghezza del motivo/pattern ripetuto
#
# Per aggiungere un nuovo meccanismo basta aggiungere un
# nuovo elemento a questa lista con questi 3 campi.
# NON bisogna modificare il resto dell'EM: tutte le funzioni
# sotto sono generiche e leggono i nomi delle colonne da qui.


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

  DI= list(
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
  ),

  DS = list(
    mechanism_col = "SD_direct_substitution",
    motif_pos_col = "SD_DS_motif_pos",
    motif_len_col = "SD_DS_repeat_pat_len",
    len_scale = 0.5
  ),
  IS = list(
    mechanism_col = "SD_inverted_substitution",
    motif_pos_col = "SD_IS_motif_pos",
    motif_len_col = "SD_IS_repeat_pat_len",
    len_scale = 0.5
  )

)


# 1. FUNZIONI BASE

# motif_strength()
#
# Calcola la "forza" complessiva di un insieme di motivi
# associati a una posizione/riga, sommando un contributo
# esponenziale per ciascun motivo:
#
#   g = sum_i exp( beta*(motif_l_i - L0) - |distance_i|/lambda )
#
# - motivi più lunghi (motif_l alto)  -> contributo maggiore
# - motivi più lontani (distance alta) -> contributo minore
# - lambda controlla quanto velocemente il contributo
#   decade con la distanza
# - beta controlla quanto pesa la lunghezza del motivo
#
# Non è chiamata direttamente nel resto del codice (che usa
# la versione vettorizzata calculate_g_all_fast), ma resta
# utile come riferimento/per calcoli singoli.
# ---------------------------------------------------------

motif_strength <- function(distance, motif_l, lambda, beta, L0 = 2) {
  stopifnot(length(distance) == length(motif_l), lambda > 0)
  sum(exp(beta * (motif_l - L0) - abs(distance) / lambda))
}

# position_probability()
#
# Dato il "punteggio" g della posizione osservata (g_obs) e
# i punteggi g di tutte le posizioni counterfactual (g_cf,
# cioè le altre posizioni possibili ma non osservate),
# restituisce la probabilità relativa che il meccanismo
# abbia agito proprio nella posizione osservata:
#
#   p = g_obs / (g_obs + sum(g_cf))
#
# Se il totale è 0 o NA (nessun segnale di motivo da
# nessuna parte) restituisce 0.

position_probability <- function(g_obs, g_cf) {
  total <- g_obs + sum(g_cf)
  if (is.na(total) || total == 0) return(0)
  g_obs / total
}

# indel_length_prob()
#
# Probabilità della lunghezza L dell'indel secondo una
# distribuzione geometrica di parametro lambda_I (stessa
# forma sia per i meccanismi con motivo sia per NHEJ, solo
# con lambda_I diversi).

indel_length_prob <- function(L, lambda_I) dgeom(L - 1, prob = lambda_I)

# compute_L()

# se ANC e DER hanno lunghezza diversa -> vero indel,L = differenza di lunghezza
# - se ANC e DER hanno la stessa lunghezza -> sostituzione/MNV, L = numero di basi che differiscono
#cioè la distanza di Hamming,sempre >= 1 perché ANC != DER per definizione di variante
compute_L <- function(anc, der) {
  la <- nchar(anc); ld <- nchar(der)
  if (la != ld) return(abs(la - ld))

  a_chars <- strsplit(anc, "")[[1]]
  d_chars <- strsplit(der, "")[[1]]
  L <- sum(a_chars != d_chars)
  if (L == 0) L <- 1  # sicurezza: ANC identico a DER (errore nei dati)
  L
}

#se ANC e DER hanno lunghezza diversa -> vero indel,L = lunghezza del maggiore fra ANC e DER
#- se ANC e DER hanno la stessa lunghezza -> sostituzione/MNV, L = lunghezza di ANC che è uguale a lunghezza di DER
max_L <- function(anc,der) {
  LA <- nchar(anc)
  LD <- nchar(der)
  if (LA >= LD) return (LA)
  else return (LD)
}


# extract_first_from_brackets
# Questa funzione riconosce la presenza delle parentesi
# quadre e restituisce SEMPRE E SOLO il primo elemento della
# lista (es. "[True, False]" -> "True", "[3, 8]" -> "3"),
# ripulito da spazi ed eventuali apici/virgolette.
#
# Se il valore non è racchiuso tra quadre, viene restituito invariato,
# Gestisce NA, stringa vuota e "-" restituendoli così come
# sono (verranno poi gestiti dalle funzioni chiamanti).
extract_first_from_brackets <- function(x) {
  if (is.na(x)) return(NA_character_)
 
  xc <- trimws(as.character(x))
  if (xc == "" || xc == "-") return(xc)
 
  if (startsWith(xc, "[") && endsWith(xc, "]")) {
    inner <- substr(xc, 2, nchar(xc) - 1)
    first <- strsplit(inner, ",")[[1]][1]
    first <- trimws(first)
    first <- gsub("^['\"]|['\"]$", "", first)  # rimuove eventuali apici/virgolette
    return(first)
  }
 
  # Nessuna quadra: valore già "semplice", lo lascio invariato
  xc
}

# 2. PRE-PARSING DEI MOTIVI
## Converte una cella (es. "[3, 8]" nel nuovo formato, oppure
# "1,4,7" nel vecchio formato) nella posizione del motivo da
# usare per quella riga.
#
# Con il nuovo formato a lista, prende SEMPRE E SOLO il primo
# valore della lista (tramite extract_first_from_brackets),
# restituendolo come vettore numerico di lunghezza 1.
#
# Gestisce i casi vuoti/mancanti ("", NA, "-") restituendo
# un vettore vuoto.

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

  # per MMEJ non ci interessa la distanza dle motivo, quindi 
  # in quel caso motif_pos_col è NULL, e non si legge nessuna posizione.
  has_position <- !is.null(motif_pos_col) && nzchar(motif_pos_col)

  # mechanism_col può arrivare come "[True, False]": estraggo il primo valore
  # della lista e confronto (case-insensitive) con "TRUE".
  mech_raw <- as.character(data[[mechanism_col]])
  mech_first <- vapply(mech_raw, extract_first_from_brackets, character(1))
  is_M <- toupper(trimws(mech_first)) == "TRUE"
  is_M[is.na(is_M)] <- FALSE

  idx <- integer(0); distance <- numeric(0); motif_l <- numeric(0)

  for (i in which(is_M)) {

    if (has_position) {
      d <- parse_motif_positions(data[[motif_pos_col]][i])
    } else {
      # il motivo "conta" indipendentemente dalla posizione,
      # quindi uso distance = 0. Nella formula di
      # motif_strength il termine -|distance|/lambda si
      # annulla sempre, e g dipende solo dalla lunghezza del
      # motivo (beta*(motif_len - L0)).
      d <- 0
    }

    # motif_len_col può arrivare come "[3, 8]": prendo sempre
    # il primo valore della lista.
    ml_raw <- extract_first_from_brackets(data[[motif_len_col]][i])
    ml <- suppressWarnings(as.numeric(ml_raw)) * len_scale

    if (length(d) == 0 || is.na(ml)) next  # riga senza motivi utilizzabili

    # "Espande" la riga i: un elemento per ogni posizione di motivo trovata
    idx <- c(idx, rep(i, length(d)))
    distance <- c(distance, d)
    motif_l <- c(motif_l, rep(ml, length(d)))
  }

  list(is_M = is_M, motif_row_idx = idx, motif_distance = distance, motif_length = motif_l)
}

# 3. CALCOLO VELOCE DI g
# Versione vettorizzata di motif_strength(), applicata a
# TUTTE le righe del dataset in un colpo solo, per una data
# coppia di parametri (lambda, beta) di un meccanismo.
#
# Per ogni motivo calcola il peso
#   w = exp(beta*(motif_length-L0) - |motif_distance|/lambda)
# poi somma i pesi per riga con rowsum() e assegna il risultato al vettore
# g (lungo quanto il dataset). Le righe senza motivi restano a 0.
# Questa funzione viene richiamata a ogni iterazione dell'EM e a ogni valutazione della funzione da ottimizzare
# (Q_function).

calculate_g_all_fast <- function(motif_info, n_data, lambda, beta, L0 = 2) {

  g <- numeric(n_data)  # default 0 per tutte le righe
  if (length(motif_info$motif_row_idx) == 0) return(g)

  w <- exp(beta * (motif_info$motif_length - L0) - abs(motif_info$motif_distance) / lambda)

  # Somma i pesi w per riga originale (motif_row_idx può ripetersi
  # più volte se una riga ha più motivi)
  s <- rowsum(w, motif_info$motif_row_idx, reorder = FALSE)
  idx <- as.integer(rownames(s))
  g[idx] <- s[, 1]
  g
}

# 4. LIKELIHOOD

# Likelihood di un meccanismo con motivo: prob. della lunghezza
# dell'indel (geometrica) * prob. della posizione osservata
# (rispetto alle counterfactual)
likelihood_mechanism <- function(L, lambda_I, g_obs, g_cf)
  indel_length_prob(L, lambda_I) * position_probability(g_obs, g_cf)

# Likelihood NHEJ (nessun motivo/microomologia): stessa forma
# per la lunghezza, ma posizione uniforme su tutte le
# n_positions disponibili (nessuna preferenza di posizione)
likelihood_NHEJ <- function(L, lambda_I_NHEJ, n_positions)
  indel_length_prob(L, lambda_I_NHEJ) * (1 / n_positions)

# precompute_variant_stats_multi()
#
# Prepara, PRIMA di iniziare l'EM, tutte le quantità che
# dipendono dai parametri iniziali:
# per ogni "variante" (variant_id) raggruppa gli allineamenti
# candidati (alignment_ID != 0) e i counterfactual
# (alignment_ID == 0, cioè le posizioni alternative non
# osservate), e per ognuno calcola:
# - la lunghezza dell'indel (L)
# - g osservato e g dei counterfactual, per ogni meccanismo
# - p_position (prob. di posizione) e L_<mech> (likelihood)
#   iniziali per ogni meccanismo
# - L_NHEJ iniziale.
# Il risultato (stats) è una tabella "long" con una riga per
# ogni alignment candidato, che verrà poi aggiornata a ogni
# iterazione dell'EM (in em_multi_mechanism) con i nuovi
# parametri.
#
# Restituisce anche motif_info (riusato per ricalcolare g in
# ogni iterazione) e cf_index_list (indici dei counterfactual
# per ogni variante, per sommare velocemente i loro g).

precompute_variant_stats_multi <- function(data, mechanisms, lambda_init, beta_init,
                                            lambda_I_init, lambda_I_NHEJ_init, L0 = 2) {

  # ---- Controllo che tutte le colonne necessarie esistano ----

  required_cols <- c("variant_id", "alignment_ID", "ANC", "DER")
  for (mech_name in names(mechanisms)) {
    mech <- mechanisms[[mech_name]]
    required_cols <- c(required_cols, mech$mechanism_col, mech$motif_pos_col, mech$motif_len_col)
  }
  required_cols <- unique(required_cols)

  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0)
    stop(paste("Mancano le colonne:", paste(missing_cols, collapse = ", ")))

  # Indice di riga originale (serve per rimappare g_all -> righe)
  data$.data_idx <- seq_len(nrow(data))

  # ---- Parsing motivi + calcolo g iniziale, per ogni meccanismo ----

  motif_info <- list()
  for (mech_name in names(mechanisms))
    motif_info[[mech_name]] <- prepare_motif_data(data, mechanisms[[mech_name]])

  for (mech_name in names(mechanisms))
    data[[paste0("g_", mech_name)]] <- calculate_g_all_fast(
      motif_info[[mech_name]], nrow(data), lambda_init[[mech_name]], beta_init[[mech_name]], L0
    )

  # ---- Costruzione della tabella "stats", variante per variante ----

  variant_ids <- unique(data$variant_id)
  out_list <- vector("list", length(variant_ids))

  for (i in seq_along(variant_ids)) {
    id <- variant_ids[i]
    group <- data[variant_id == id]

    cf <- group[alignment_ID == 0]          # posizioni counterfactual (non osservate)
    align_rows <- group[alignment_ID != 0]  # allineamenti candidati osservati

    if (nrow(align_rows) == 0) {
      warning(paste("variant_id", id, "non ha alignment candidati."))
      next
    }

    # Numero di posizioni valide per questa variante: preso dal
    # primo meccanismo, perché tutti i meccanismi condividono la
    # stessa struttura di counterfactual (stesse posizioni possibili)
    first_mech <- names(mechanisms)[1]
    g_first <- group[[paste0("g_", first_mech)]]
    n_positions <- sum(!is.na(g_first))

    align_list <- vector("list", nrow(align_rows))

    # Un blocco per ogni allineamento candidato di questa variante
    for (j in seq_len(nrow(align_rows))) {
      row <- align_rows[j]
      L <- max_L(row$ANC, row$DER)  # lunghezza dell'indel, se ANC = DER allora distanza di hamming 

      # Colonne comuni a tutti i meccanismi
      result_row <- data.frame(
        variant_id = id,
        alignment_ID = row$alignment_ID,
        indel_length = L,
        n_positions = n_positions,
        data_idx = row$.data_idx
      )

      # Colonne specifiche per ogni meccanismo (g, p_position, L)
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

      # Colonna NHEJ (nessun motivo coinvolto)
      result_row$L_NHEJ <- likelihood_NHEJ(L, lambda_I_NHEJ_init, n_positions)

      align_list[[j]] <- result_row
    }

    out_list[[i]] <- do.call(rbind, align_list)
  }

  # Unisce tutte le varianti in un'unica tabella "long"
  stats <- do.call(rbind, out_list)

  # Assegna a ogni variant_id un indice numerico progressivo
  # (usato più avanti per operazioni vettorizzate per gruppo)
  variant_levels <- unique(stats$variant_id)
  stats$variant_group <- match(stats$variant_id, variant_levels)

  # Per ogni variante, gli indici (nel dataset originale) delle
  # righe counterfactual: servono per sommare velocemente i loro
  # g ad ogni iterazione, senza rifare il subset ogni volta
  cf_index_list <- lapply(variant_levels, function(id)
    data[variant_id == id & alignment_ID == 0, .data_idx])

  list(stats = stats, motif_info = motif_info, cf_index_list = cf_index_list, variant_levels = variant_levels)
}

# 6. OTTIMIZZAZIONE lambda E beta (M-STEP, un meccanismo)

# optimize_lambda_beta_mechanism()
#
# Dato un meccanismo e le responsibility correnti (P_old,
# cioè quanto l'EM "crede" che ogni allineamento sia dovuto
# a questo meccanismo, calcolate nell'E-step), trova i valori
# di lambda e beta che massimizzano la Q-function pesata:
#
#   Q(lambda,beta) = sum_i P_old_i * log( p_position_i(lambda,beta) )
#
# cioè la parte della log-verosimiglianza attesa (Q dell'EM)
# che dipende da lambda e beta, pesata dalle responsibility.
# L'ottimizzazione è vincolata (lambda>0, beta>=0) e usa
# L-BFGS-B tramite optimx.


optimize_lambda_beta_mechanism <- function(data, stats, P_old, motif_info, cf_index_list,
                                            lambda_init, beta_init, L0 = 2) {

  n_data <- nrow(data)
  obs_data_idx <- stats$data_idx
  variant_group <- stats$variant_group

  # Funzione obiettivo da MINIMIZZARE (quindi -Q): valutata più
  # volte da optimx durante la ricerca dell'ottimo
  Q_function <- function(par) {
    lambda <- par[1]; beta <- par[2]

    # Parametri fuori dai vincoli -> penalità enorme
    if (!is.finite(lambda) || !is.finite(beta) || lambda <= 0 || beta < 0) return(1e100)

    # Ricalcola g per TUTTE le righe con i parametri correnti (par)
    g_all <- calculate_g_all_fast(motif_info, n_data, lambda, beta, L0)
    g_obs <- g_all[obs_data_idx]

    # Somma dei g dei counterfactual, per ogni variante, poi
    # riportata a livello di ogni singolo allineamento (stats)
    cf_sum_by_variant <- vapply(cf_index_list, function(idx)
      if (length(idx) == 0) 0 else sum(g_all[idx]), numeric(1))
    g_cf_sum <- cf_sum_by_variant[variant_group]

    # Probabilità di posizione per ogni allineamento
    total <- g_obs + g_cf_sum
    p_position <- numeric(nrow(stats))
    valid_total <- is.finite(g_obs) & is.finite(total) & total > 0
    p_position[valid_total] <- g_obs[valid_total] / total[valid_total]

    # Scarta le righe non utilizzabili nel logaritmo (p=0, NA, ecc.)
    valid <- P_old > 0 & is.finite(P_old) & p_position > 0 & is.finite(p_position)
    if (!any(valid)) return(1e100)

    # Q pesata dalle responsibility; si minimizza -Q
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

# Cuore dell'algoritmo: alterna E-step e M-step finché i
# parametri non convergono (o si raggiunge max_iter).
#
# Può lavorare con un numero qualsiasi di meccanismi definiti
# in `mechanisms` (es. solo ID, solo DD, oppure ID+DD), più
# la componente NHEJ che è sempre presente come alternativa
# "di base" (nessun motivo).

# Modello concettuale: ogni allineamento osservato può essere
# generato da uno dei meccanismi (con probabilità pi_<mech>)
# oppure da NHEJ (probabilità pi_NHEJ = 1 - sum(pi)).

em_multi_mechanism <- function(data, stats, motif_info, cf_index_list, mechanisms,
                                pi_init, lambda_init, beta_init, lambda_I_init,
                                lambda_I_NHEJ_init, max_iter = 100, tol = 1e-5,
                                L0 = 2, verbose = TRUE) {

  mechanism_names <- names(mechanisms)
  n_mechanisms <- length(mechanism_names)

  # ---- Controlli di coerenza sugli input ----

  if (length(pi_init) != n_mechanisms) stop("Numero di pi_init non compatibile con i meccanismi.")
  if (length(lambda_init) != n_mechanisms) stop("Numero di lambda_init non compatibile.")
  if (length(beta_init) != n_mechanisms) stop("Numero di beta_init non compatibile.")
  if (length(lambda_I_init) != n_mechanisms) stop("Numero di lambda_I_init non compatibile.")

  # ---- Parametri correnti (verranno aggiornati a ogni iterazione) ----

  pi <- pi_init; names(pi) <- mechanism_names
  pi_NHEJ <- 1 - sum(pi)
  if (pi_NHEJ <= 0) stop("La somma dei pi iniziali deve essere < 1.")

  lambda <- lambda_init   # decadimento con la distanza, per meccanismo
  beta <- beta_init       # peso della lunghezza del motivo, per meccanismo
  lambda_I <- lambda_I_init           # parametro geometrico lunghezza indel, per meccanismo
  lambda_I_NHEJ <- lambda_I_NHEJ_init # idem per NHEJ

  # ---- Quantità comuni, calcolate una sola volta ----

  p_position_NHEJ <- 1 / stats$n_positions  # NHEJ: nessuna preferenza di posizione
  L_vec <- stats$indel_length
  v_id <- stats$variant_id
  first_of_variant <- !duplicated(v_id)  # una riga rappresentativa per variante (per somme/medie per variante)

  history <- data.frame(iter = integer(0), Loglik = numeric(0))
  eps <- 1e-6  # margine per tenere lambda_I dentro (0,1)

  # ===================== ITERAZIONI EM =====================

  for (iter in seq_len(max_iter)) {

    # ---------------- E-STEP ----------------
    # Calcola, con i parametri CORRENTI, la likelihood di ogni
    # meccanismo e di NHEJ per ogni allineamento, poi ne deriva
    # le responsibility (probabilità a posteriori di ciascun
    # meccanismo per ogni allineamento).

    # g aggiornato per ogni meccanismo, con lambda/beta correnti
    g_all <- list()
    for (mech_name in mechanism_names)
      g_all[[mech_name]] <- calculate_g_all_fast(
        motif_info[[mech_name]], nrow(data), lambda[[mech_name]], beta[[mech_name]], L0
      )

    L_mechanisms <- list()

    for (mech_name in mechanism_names) {
      g_obs <- g_all[[mech_name]][stats$data_idx]

      # Somma dei g dei counterfactual per variante, riportata a
      # livello di ogni allineamento
      cf_sum_by_variant <- vapply(cf_index_list, function(idx)
        if (length(idx) == 0) 0 else sum(g_all[[mech_name]][idx]), numeric(1))
      g_cf_sum <- cf_sum_by_variant[stats$variant_group]

      total_position <- g_obs + g_cf_sum
      p_position <- numeric(nrow(stats))
      valid_total <- is.finite(g_obs) & is.finite(total_position) & total_position > 0
      p_position[valid_total] <- g_obs[valid_total] / total_position[valid_total]

      stats[[paste0("p_position_", mech_name)]] <- p_position
      stats[[paste0("g_", mech_name)]] <- g_obs

      # Likelihood del meccanismo = P(lunghezza) * P(posizione)
      p_length <- dgeom(L_vec - 1, prob = lambda_I[[mech_name]])
      L_mech <- p_length * p_position

      L_mechanisms[[mech_name]] <- L_mech
      stats[[paste0("L_", mech_name)]] <- L_mech
    }

    # Likelihood NHEJ: P(lunghezza) * P(posizione uniforme)
    p_length_NHEJ <- dgeom(L_vec - 1, prob = lambda_I_NHEJ)
    L_NHEJ <- p_length_NHEJ * p_position_NHEJ
    stats$L_NHEJ <- L_NHEJ

    # ---------- POSTERIOR / RESPONSIBILITIES ----------
    # Per ogni allineamento: numeratore = pi_k * L_k (per ogni
    # meccanismo k e per NHEJ). Le responsibility si ottengono
    # normalizzando sulla somma di TUTTI i meccanismi e su TUTTI
    # gli allineamenti della stessa variante (perché la variante,
    # nel suo complesso, deve avere probabilità totale 1).

    num_mechanisms <- list()
    for (mech_name in mechanism_names)
      num_mechanisms[[mech_name]] <- pi[[mech_name]] * L_mechanisms[[mech_name]]

    num_NHEJ <- pi_NHEJ * L_NHEJ

    # Contributo totale (non normalizzato) di ogni singolo allineamento
    variant_contribution <- num_NHEJ
    for (mech_name in mechanism_names)
      variant_contribution <- variant_contribution + num_mechanisms[[mech_name]]

    # Somma dei contributi di tutti gli allineamenti della stessa
    # variante = denominatore di normalizzazione (per variante)
    variant_tot <- ave(variant_contribution, v_id, FUN = sum)

    # Responsibility di ogni meccanismo, per ogni allineamento
    P <- list()
    for (mech_name in mechanism_names) {
      P_mech <- ifelse(is.na(variant_tot) | variant_tot == 0, NA_real_,
                        num_mechanisms[[mech_name]] / variant_tot)
      P[[mech_name]] <- P_mech
      stats[[paste0("P_", mech_name)]] <- P_mech
    }

    P_NHEJ <- ifelse(is.na(variant_tot) | variant_tot == 0, NA_real_, num_NHEJ / variant_tot)
    stats$P_NHEJ <- P_NHEJ

    # Log-verosimiglianza osservata: una sola volta per variante
    # (variant_tot è lo stesso per tutte le righe della variante,
    # quindi si prende solo first_of_variant)
    Loglik <- sum(log(variant_tot[first_of_variant]), na.rm = TRUE)

    # ---------------- M-STEP ----------------
    # Aggiorna tutti i parametri usando le responsibility appena
    # calcolate nell'E-step.

    # 1. lambda e beta di ogni meccanismo (ottimizzazione numerica)
    opt_results <- list()
    for (mech_name in mechanism_names)
      opt_results[[mech_name]] <- optimize_lambda_beta_mechanism(
        data = data, stats = stats, P_old = P[[mech_name]],
        motif_info = motif_info[[mech_name]], cf_index_list = cf_index_list,
        lambda_init = lambda[[mech_name]], beta_init = beta[[mech_name]], L0 = L0
      )

    lambda_new <- lambda; beta_new <- beta
    for (mech_name in mechanism_names) {
      lambda_new[[mech_name]] <- opt_results[[mech_name]]$lambda
      beta_new[[mech_name]] <- opt_results[[mech_name]]$beta
    }

    # 2. Aggiornamento dei prior pi (media delle responsibility
    #    per variante, sommate sugli allineamenti della variante)
    pi_new <- numeric(n_mechanisms); names(pi_new) <- mechanism_names
    for (mech_name in mechanism_names) {
      P_per_variant <- ave(P[[mech_name]], v_id, FUN = sum)
      pi_new[[mech_name]] <- mean(P_per_variant[first_of_variant], na.rm = TRUE)
    }

    P_NHEJ_per_variant <- ave(P_NHEJ, v_id, FUN = sum)
    pi_NHEJ_new <- mean(P_NHEJ_per_variant[first_of_variant], na.rm = TRUE)

    # 3. Aggiornamento lambda_I di ogni meccanismo: stima di
    #    massima verosimiglianza pesata (M-step chiuso, in forma
    #    analitica, per la parte geometrica della lunghezza)
    lambda_I_new <- lambda_I
    for (mech_name in mechanism_names) {
      P_mech <- P[[mech_name]]
      lambda_I_new[[mech_name]] <- sum(P_mech, na.rm = TRUE) / sum(P_mech * L_vec, na.rm = TRUE)
      lambda_I_new[[mech_name]] <- min(max(lambda_I_new[[mech_name]], eps), 1 - eps)  # tieni in (0,1)
    }

    # 4. Aggiornamento lambda_I di NHEJ, stessa logica
    lambda_I_NHEJ_new <- sum(P_NHEJ, na.rm = TRUE) / sum(P_NHEJ * L_vec, na.rm = TRUE)
    lambda_I_NHEJ_new <- min(max(lambda_I_NHEJ_new, eps), 1 - eps)

    # ---------------- CONVERGENZA ----------------
    # delta = variazione massima assoluta tra vecchi e nuovi
    # parametri (pi, lambda, beta, lambda_I, per ogni meccanismo,
    # più pi_NHEJ e lambda_I_NHEJ). Se scende sotto tol, l'EM
    # si ferma.

    delta_values <- abs(pi_new - pi)
    for (mech_name in mechanism_names)
      delta_values <- c(delta_values,
        abs(lambda_new[[mech_name]] - lambda[[mech_name]]),
        abs(beta_new[[mech_name]] - beta[[mech_name]]),
        abs(lambda_I_new[[mech_name]] - lambda_I[[mech_name]])
      )
    delta_values <- c(delta_values, abs(pi_NHEJ_new - pi_NHEJ), abs(lambda_I_NHEJ_new - lambda_I_NHEJ))

    delta <- max(delta_values, na.rm = TRUE)

    # ---------------- HISTORY ----------------
    # Salva una riga per iterazione con tutti i parametri e la
    # log-verosimiglianza, utile per controllare la convergenza
    # a posteriori (es. plottando Loglik o delta per iterazione).

    history_row <- data.frame(iter = iter, Loglik = Loglik)
    for (mech_name in mechanism_names)
      history_row[[paste0("pi_", mech_name)]] <- pi_new[[mech_name]]
    history_row$pi_NHEJ <- pi_NHEJ_new

    for (mech_name in mechanism_names) {
      history_row[[paste0("lambda_", mech_name)]] <- lambda_new[[mech_name]]
      history_row[[paste0("beta_", mech_name)]] <- beta_new[[mech_name]]
      history_row[[paste0("lambda_I_", mech_name)]] <- lambda_I_new[[mech_name]]
    }
    history_row$lambda_I_NHEJ <- lambda_I_NHEJ_new
    history_row$delta <- delta

    history <- rbind(history, history_row)

    # ---------------- STAMPA (log a schermo) ----------------

    if (verbose) {
      cat(sprintf("Iter %3d:", iter))
      for (mech_name in mechanism_names)
        cat(sprintf(" pi_%s=%.6f lambda_%s=%.6f beta_%s=%.6f lambda_I_%s=%.6f",
                    mech_name, pi_new[[mech_name]], mech_name, lambda_new[[mech_name]],
                    mech_name, beta_new[[mech_name]], mech_name, lambda_I_new[[mech_name]]))
      cat(sprintf(" pi_NHEJ=%.6f lambda_I_NHEJ=%.6f delta=%.2e Loglik=%.6f\n",
                  pi_NHEJ_new, lambda_I_NHEJ_new, delta, Loglik))
    }

    # ---------------- AGGIORNAMENTO PARAMETRI ----------------
    # I parametri "new" diventano quelli correnti per la
    # prossima iterazione (o per l'output finale, se si esce
    # dal ciclo).

    pi <- pi_new; pi_NHEJ <- pi_NHEJ_new
    lambda <- lambda_new; beta <- beta_new
    lambda_I <- lambda_I_new; lambda_I_NHEJ <- lambda_I_NHEJ_new

    # ---------------- TEST DI CONVERGENZA ----------------

    if (delta < tol) {
      if (verbose) cat("tolleranza raggiunta in iter =", iter, "\n")
      break
    }
  }

  # ===================== RISULTATI FINALI =====================
  # stats contiene, per ogni allineamento, le colonne g_/L_/P_
  # finali per ogni meccanismo: è la tabella che viene poi
  # scritta su disco.

  list(results = stats, pi = pi, pi_NHEJ = pi_NHEJ, lambda = lambda, beta = beta,
       lambda_I = lambda_I, lambda_I_NHEJ = lambda_I_NHEJ, history = history,
       n_iter = iter, converged = delta < tol)
}

# 8. FUNZIONE PRINCIPALE

# run_em_multi()
#
# Punto di ingresso "ad alto livello": riempie i parametri
# iniziali mancanti con dei default ragionevoli, esegue il
# pre-calcolo (sezione 5) e poi lancia l'EM (sezione 7).
#
# `mechanisms` può contenere: ID, oppure DD, oppure ID + DD
# (o qualunque sottoinsieme di MECHANISMS).

run_em_multi <- function(data, mechanisms, lambda_init = NULL, beta_init = NULL, pi_init = NULL,
                          lambda_I_init = NULL, lambda_I_NHEJ_init = 0.3, max_iter = 100,
                          tol = 1e-8, L0 = 2, verbose = TRUE) {

  mechanism_names <- names(mechanisms)
  n_mechanisms <- length(mechanism_names)

  if (n_mechanisms == 0) stop("Nessun meccanismo selezionato.")

  # ---- Valori iniziali di default, se non forniti ----

  if (is.null(lambda_init)) { lambda_init <- rep(50, n_mechanisms); names(lambda_init) <- mechanism_names }
  if (is.null(beta_init)) { beta_init <- rep(0.7, n_mechanisms); names(beta_init) <- mechanism_names }
  if (is.null(lambda_I_init)) { lambda_I_init <- rep(0.3, n_mechanisms); names(lambda_I_init) <- mechanism_names }

  # Prior iniziali: distribuiamo 0.5 tra i meccanismi selezionati
  # e lasciamo 0.5 a NHEJ (es. ID+DD -> pi_ID=0.25, pi_DD=0.25, pi_NHEJ=0.50)
  if (is.null(pi_init)) { pi_init <- rep(0.5 / n_mechanisms, n_mechanisms); names(pi_init) <- mechanism_names }

  if (sum(pi_init) >= 1) stop("La somma dei prior dei meccanismi deve essere < 1.")

  # ---- Pre-calcolo delle quantità fisse (sezione 5) ----

  cat("Pre-calcolo delle quantità fisse...\n")

  pre <- precompute_variant_stats_multi(
    data = data, mechanisms = mechanisms, lambda_init = lambda_init,
    beta_init = beta_init, lambda_I_init = lambda_I_init,
    lambda_I_NHEJ_init = lambda_I_NHEJ_init, L0 = L0
  )

  stats <- pre$stats; motif_info <- pre$motif_info; cf_index_list <- pre$cf_index_list

  cat("Pre-calcolo completato:", nrow(stats), "righe di alignment.\n\n")

  cat("Esecuzione EM su:")
  for (mech_name in mechanism_names) cat(mech_name, " ")
  cat("+ NHEJ\n\n")

  # ---- Esecuzione dell'EM vero e proprio (sezione 7) ----

  em_multi_mechanism(
    data = data, stats = stats, motif_info = motif_info, cf_index_list = cf_index_list,
    mechanisms = mechanisms, pi_init = pi_init, lambda_init = lambda_init, beta_init = beta_init,
    lambda_I_init = lambda_I_init, lambda_I_NHEJ_init = lambda_I_NHEJ_init,
    max_iter = max_iter, tol = tol, L0 = L0, verbose = verbose
  )
}

# 9. CLI (interfaccia da riga di comando)

# Questo blocco viene eseguito SOLO quando lo script è
# lanciato direttamente con Rscript (sys.nframe() == 0), non
# quando le funzioni sopra vengono sourceate/importate.
#
# Esempi:
#   ID vs NHEJ:       Rscript EM_multi.R -i input.tsv -o output.tsv --mechanisms ID
#   DD vs NHEJ:       Rscript EM_multi.R -i input.tsv -o output.tsv --mechanisms DD
#   ID vs DD vs NHEJ: Rscript EM_multi.R -i input.tsv -o output.tsv --mechanisms ID,DD

if (sys.nframe() == 0) {

  args <- commandArgs(trailingOnly = TRUE)

  # Legge il valore associato a un flag da riga di comando
  # (es. get_arg("-i") legge il valore dopo "-i")
  get_arg <- function(flag, default = NULL, numeric = FALSE) {
    pos <- which(args == flag)
    if (length(pos) == 0) return(default)
    if (pos[1] >= length(args)) return(default)
    val <- args[pos[1] + 1]
    if (numeric) val <- as.numeric(val)
    val
  }

  # ---- Lettura argomenti da riga di comando ----

  input_file <- get_arg("-i")
  output_file <- get_arg("-o")
  mechanisms_arg <- get_arg("--mechanisms")
  pi_arg <- get_arg("-p")
  tol <- get_arg("-t", 1e-8, TRUE)
  max_iter <- get_arg("-m", 100, TRUE)

  if (is.null(input_file) || is.null(output_file))
    stop("Uso:\nRscript EM_multi.R -i input.tsv -o output.tsv --mechanisms ID,DD [-p pi1,pi2] [-t tol] [-m max_iter]")

  if (is.null(mechanisms_arg)) stop("Devi specificare --mechanisms.")

  # ---- Parsing e validazione dei meccanismi richiesti ----

  selected_names <- trimws(unlist(strsplit(mechanisms_arg, ",")))

  unknown <- setdiff(selected_names, names(MECHANISMS))
  if (length(unknown) > 0)
    stop(paste("Meccanismi non riconosciuti:", paste(unknown, collapse = ", "),
               "\nMeccanismi disponibili:", paste(names(MECHANISMS), collapse = ", ")))

  mechanisms <- MECHANISMS[selected_names]

  cat("Input :", input_file, "\n")
  cat("Output:", output_file, "\n")
  cat("Meccanismi:")
  cat(paste(selected_names, collapse = ", "), "+ NHEJ\n\n")

  # ---- Importazione dati ----

  cat("Importazione del TSV...\n")
  data <- fread(input_file, sep = "\t", header = TRUE)
  cat("Importazione completata:", nrow(data), "righe x", ncol(data), "colonne\n\n")

  # ---- Parametri iniziali (di default, salvo -p per i prior) ----

  n_mechanisms <- length(selected_names)

  lambda_init <- rep(50, n_mechanisms); names(lambda_init) <- selected_names
  beta_init <- rep(0.7, n_mechanisms); names(beta_init) <- selected_names
  lambda_I_init <- rep(0.3, n_mechanisms); names(lambda_I_init) <- selected_names

  # Prior iniziali: -p 0.25,0.25 -> pi_ID=0.25, pi_DD=0.25, pi_NHEJ=0.50
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
    data = data, mechanisms = mechanisms, lambda_init = lambda_init, beta_init = beta_init,
    pi_init = pi_init, lambda_I_init = lambda_I_init, lambda_I_NHEJ_init = 0.3,
    max_iter = max_iter, tol = tol, L0 = 2, verbose = TRUE
  )

  # ---- Stampa risultati finali a schermo ----

  cat("\n==============================\nPRIOR FINALI\n==============================\n")
  for (mech_name in selected_names) cat(mech_name, "=", em_out$pi[[mech_name]], "\n")
  cat("NHEJ =", em_out$pi_NHEJ, "\n")

  cat("\n==============================\nPARAMETRI FINALI\n==============================\n")
  for (mech_name in selected_names) {
    cat("\n", mech_name, "\n")
    cat("lambda =", em_out$lambda[[mech_name]], "\n")
    cat("beta =", em_out$beta[[mech_name]], "\n")
    cat("lambda_I =", em_out$lambda_I[[mech_name]], "\n")
  }

  cat("\nNHEJ\n")
  cat("lambda_I_NHEJ =", em_out$lambda_I_NHEJ, "\n")

  cat("\nIterazioni:", em_out$n_iter, "\n")
  cat("Convergenza:", em_out$converged, "\n\n")

  # ---- Scrittura risultati su file ----

  cat("Scrittura dei risultati...\n")
  fwrite(em_out$results, file = output_file, sep = "\t", quote = FALSE, na = "NA")
  cat("Output salvato in:", output_file, "\n")

    # ---- Secondo output: riepilogo dell'ultima iterazione ----
  # em_out$history ha una riga per ogni iterazione dell'EM; qui ne
  # teniamo solo l'ultima (tail(..., 1)), cioè Loglik, delta e tutti
  # i parametri (pi, lambda, beta, lambda_I) così come aggiornati
  # all'ultima iterazione. Il nome del file è derivato da output_file
  # aggiungendo un suffisso, così non serve un nuovo argomento CLI.

  last_iter_file <- paste0(output_file, "_last_iter")
  last_iter <- tail(em_out$history, 1)

  cat("Scrittura della history...\n")
  fwrite(last_iter, file = last_iter_file, sep = "\t", quote = FALSE, na = "NA")
  cat("last iter salvata in:", last_iter_file, "\n")
}



