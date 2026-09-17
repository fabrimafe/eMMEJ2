#terzo_EM
library(data.table)

# 1. MOTIF STRENGTH
motif_strength <- function(distance, motif_l, lambda, beta, L0 = 2) {

    stopifnot(length(distance) == length(motif_l))
    stopifnot(lambda > 0)

    w <- exp(
        beta * (motif_l - L0) -
        abs(distance) / lambda
    )
    g <- sum(w)

    return(g)
}

# 2. POSITION PROBABILITY
position_probability <- function(g_obs, g_cf) {

    total <- g_obs + sum(g_cf)

    if (is.na(total) || total == 0) {
        return(0)
    }
    g_obs / total
}

# 3. PROBABILITÀ LUNGHEZZA INDEL
indel_length_prob <- function(L, lambda_I) {
    dgeom(
        L - 1,
        prob = lambda_I
    )
}

# n. posizioni valide per variant_id
calculate_n_positions <- function(group) {
    if (!"g_M" %in% colnames(group)) {
        stop("La colonna non è presente nel gruppo.")
    }
    n_valid <- sum(!is.na(group$g_M))
    return(n_valid)
}

# 4. PARSING POSIZIONI MOTIVI
parse_motif_positions <- function(x) {

    if (is.na(x) || x == "" || x == "-") {
        return(numeric(0))
    }

    positions <- as.numeric(unlist(strsplit(as.character(x), ",")))
    positions <- positions[!is.na(positions)]

    return(positions)
}

# 5. MOTIF STRENGTH DI UNA SINGOLA RIGA
calculate_M_strength <- function(row, lambda, beta, L0 = 2) {

    is_M <- trimws(row[["SD_inverted_deletion"]])

    if (is.na(is_M) || is_M != "TRUE") {
        return(0)
    }

    distance <- parse_motif_positions(row[["SD_ID_motif_pos"]])

    if (length(distance) == 0) {
        return(0)
    }

    motif_lengths <- as.numeric(row[["SD_ID_repeat_pat_len"]])

    if (is.na(motif_lengths)) {
        return(0)
    }

    motif_l <- rep(motif_lengths, length(distance))

    g <- motif_strength(distance = distance, motif_l = motif_l, lambda = lambda, beta = beta, L0 = L0)

    return(g)
}
# RAPPORTO OBSERVED / COUNTERFACTUAL
obs_cf_ratio <- function(g_obs, g_cf) {
    cf_mean <- mean(g_cf, na.rm = TRUE)
    if (is.na(g_obs) || is.na(cf_mean) || cf_mean == 0) {
        return(NA_real_)
    }
    ratio <- g_obs / cf_mean
    return(ratio)
}

# Z-SCORE
motif_z_score <- function(g_obs, g_cf) {
    cf_mean <- mean(g_cf, na.rm = TRUE)
    cf_sd <- sd(g_cf, na.rm = TRUE)
    if (is.na(g_obs) || is.na(cf_mean) || is.na(cf_sd) || cf_sd == 0) {
        return(NA_real_)
    }
    z <- (g_obs - cf_mean) / cf_sd
    return(z)
}

# LIKELIHOOD DI M
likelihood_M <- function(L, lambda_I_M, g_obs, g_cf) {

    p_length_M <- indel_length_prob(L = L, lambda_I = lambda_I_M)
    p_position_M <- position_probability(g_obs = g_obs, g_cf = g_cf)

    L_M <- p_length_M * p_position_M

    return(L_M)
}

# LIKELIHOOD DI NHEJ
likelihood_NHEJ <- function(L, lambda_I_NHEJ, n_positions) {
    p_length_NHEJ <- indel_length_prob(L = L, lambda_I = lambda_I_NHEJ)
    p_position_NHEJ <- 1 / n_positions
    L_NHEJ <- p_length_NHEJ * p_position_NHEJ
    return(L_NHEJ)
}

# POSTERIOR CON PRIOR ESPLICITO (generalizzazione, prior non
# più fisso a 0.5/0.5 ma parametro pi_M / pi_NHEJ)
posterior_M_NHEJ_prior <- function(L_M, L_NHEJ, pi_M, pi_NHEJ) {

    num_M <- pi_M * L_M
    num_NHEJ <- pi_NHEJ * L_NHEJ
    total <- num_M + num_NHEJ

    if (is.na(total) || total == 0) {
        return(list(P_M = NA_real_, P_NHEJ = NA_real_))
    }

    list(P_M = num_M / total, P_NHEJ = num_NHEJ / total)
}
# =========================================================
# PRE-CALCOLO MULTI-ALIGNMENT (v2)
# Discriminante = alignment_ID, non observed:
#  - alignment_ID == 0  -> pool controfattuale
#  - alignment_ID != 0  -> alignment candidati (uno per ciascun ID > 0)
# n_positions = totale righe del gruppo (cf + alignment)
# =========================================================
precompute_variant_stats_multi <- function(data, lambda_M, beta_M,
                                            lambda_I_M, lambda_I_NHEJ, L0 = 2) {

    required_cols <- c("variant_id", "alignment_ID", "ANC", "DER",
                        "SD_inverted_deletion", "SD_ID_motif_pos", "SD_ID_repeat_pat_len")
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
        stop(paste("Mancano le colonne:", paste(missing_cols, collapse = ", ")))
    }

    data$g_M <- apply(data, 1, calculate_M_strength,
                       lambda = lambda_M, beta = beta_M, L0 = L0)

    variant_ids <- unique(data$variant_id)
    out_list <- vector("list", length(variant_ids))

    for (i in seq_along(variant_ids)) {

        current_id <- variant_ids[i]
        group <- data[data$variant_id == current_id, ]

        # pool controfattuale: alignment_ID == 0
        cf <- group[group$alignment_ID == 0, ]
        g_cf <- cf$g_M

        # n_positions = TUTTE le righe del gruppo (cf + alignment candidati)
        n_positions <- sum(!is.na(group$g_M))

        # alignment candidati: alignment_ID != 0 (indipendentemente da observed)
        align_rows <- group[group$alignment_ID != 0, ]

        if (nrow(align_rows) == 0) {
            warning(paste("variant_id", current_id, "non ha alignment candidati (alignment_ID != 0)."))
            next
        }
        align_list <- vector("list", nrow(align_rows))

        for (j in seq_len(nrow(align_rows))) {
            row <- align_rows[j, ]
            g_obs <- row$g_M
            indel_length <- nchar(row$ANC) - nchar(row$DER)
            p_position_M <- position_probability(g_obs = g_obs, g_cf = g_cf)

            L_M <- likelihood_M(L = indel_length, lambda_I_M = lambda_I_M,
                                 g_obs = g_obs, g_cf = g_cf)
            L_NHEJ <- likelihood_NHEJ(L = indel_length, lambda_I_NHEJ = lambda_I_NHEJ,
                                      n_positions = n_positions)

            align_list[[j]] <- data.frame(
                variant_id = current_id,
                alignment_ID = row$alignment_ID,
                indel_length = indel_length,
                g_obs = g_obs,
                n_counterfactuals = length(g_cf),
                n_positions = n_positions,
                p_position_M = p_position_M,
                L_M = L_M,
                L_NHEJ = L_NHEJ
            )
        }

        out_list[[i]] <- do.call(rbind, align_list)
    }

    results <- do.call(rbind, out_list)
    return(results)
}

# EM SUI PRIOR (pi_M, pi_NHEJ) + lambda_I_M / lambda_I_NHEJ
# generalizzato per K alignment per variante:
# la posterior normalizza su TUTTE le ipotesi (mecc. x alignment)
# della STESSA variante, non riga per riga.

em_prior_lambda_M_NHEJ_multi <- function(stats,
                                          pi_M_init = 0.5,
                                          lambda_I_M_init = 0.3,
                                          lambda_I_NHEJ_init = 0.3,
                                          max_iter = 100, tol = 1e-5, verbose = TRUE) {

    stopifnot(all(c("variant_id", "p_position_M", "n_positions", "indel_length") %in% colnames(stats)))
    stopifnot(pi_M_init > 0 && pi_M_init < 1)
    stopifnot(lambda_I_M_init > 0 && lambda_I_M_init < 1)
    stopifnot(lambda_I_NHEJ_init > 0 && lambda_I_NHEJ_init < 1)

    pi_M <- pi_M_init
    pi_NHEJ <- 1 - pi_M_init
    lambda_I_M <- lambda_I_M_init
    lambda_I_NHEJ <- lambda_I_NHEJ_init

    p_position_NHEJ <- 1 / stats$n_positions
    L_vec <- stats$indel_length
    v_id <- stats$variant_id
    first_of_variant <- !duplicated(v_id)   # per contare ogni variante una sola volta

    history <- data.frame(iter = integer(0), pi_M = numeric(0), pi_NHEJ = numeric(0),
                           lambda_I_M = numeric(0), lambda_I_NHEJ = numeric(0),
                           delta = numeric(0), Loglik = numeric(0))

    eps <- 1e-6

    for (iter in seq_len(max_iter)) {

        # ---- ricalcolo L_M, L_NHEJ con i lambda correnti ----
        p_length_M    <- dgeom(L_vec - 1, prob = lambda_I_M)
        p_length_NHEJ <- dgeom(L_vec - 1, prob = lambda_I_NHEJ)
        L_M    <- p_length_M    * stats$p_position_M
        L_NHEJ <- p_length_NHEJ * p_position_NHEJ

        # ---- E-step: posterior normalizzata su TUTTI gli alignment della variante ----
        joint_M    <- pi_M    * L_M
        joint_NHEJ <- pi_NHEJ * L_NHEJ

        # somma per variante (broadcast: ogni riga della stessa variante ottiene lo stesso totale)
        variant_evidence <- ave(joint_M + joint_NHEJ, v_id, FUN = sum)

        P_M    <- ifelse(is.na(variant_evidence) | variant_evidence == 0, NA_real_, joint_M    / variant_evidence)
        P_NHEJ <- ifelse(is.na(variant_evidence) | variant_evidence == 0, NA_real_, joint_NHEJ / variant_evidence)

        # Loglik: una sola volta per variante (i totali sono ripetuti sulle righe della stessa variante)
        Loglik <- sum(log(variant_evidence[first_of_variant]), na.rm = TRUE)

        # ---- M-step ----
        # pi_M: media, PER VARIANTE, della massa posterior totale assegnata a M
        #       (somma di P_M sugli alignment di quella variante), non media riga per riga
        pM_per_variant    <- ave(P_M,    v_id, FUN = sum)
        pNHEJ_per_variant <- ave(P_NHEJ, v_id, FUN = sum)

        pi_M_new    <- mean(pM_per_variant[first_of_variant],    na.rm = TRUE)
        pi_NHEJ_new <- mean(pNHEJ_per_variant[first_of_variant], na.rm = TRUE)

        # lambda: media pesata su TUTTE le coppie (variante, alignment)
        lambda_I_M_new    <- sum(P_M,    na.rm = TRUE) / sum(P_M    * L_vec, na.rm = TRUE)
        lambda_I_NHEJ_new <- sum(P_NHEJ, na.rm = TRUE) / sum(P_NHEJ * L_vec, na.rm = TRUE)

        lambda_I_M_new    <- min(max(lambda_I_M_new, eps), 1 - eps)
        lambda_I_NHEJ_new <- min(max(lambda_I_NHEJ_new, eps), 1 - eps)

        delta <- max(
            abs(pi_M_new - pi_M),
            abs(lambda_I_M_new - lambda_I_M),
            abs(lambda_I_NHEJ_new - lambda_I_NHEJ)
        )

        history <- rbind(
            history,
            data.frame(iter = iter, pi_M = pi_M_new, pi_NHEJ = pi_NHEJ_new,
                       lambda_I_M = lambda_I_M_new, lambda_I_NHEJ = lambda_I_NHEJ_new,
                       delta = delta, Loglik = Loglik)
        )

        if (verbose) {
            cat(sprintf("Iter %3d: pi_M = %.6f  lambda_I_M = %.6f  lambda_I_NHEJ = %.6f  delta = %.2e  Loglik = %.6f\n",
                        iter, pi_M_new, lambda_I_M_new, lambda_I_NHEJ_new, delta, Loglik))
        }

        pi_M <- pi_M_new
        pi_NHEJ <- pi_NHEJ_new
        lambda_I_M <- lambda_I_M_new
        lambda_I_NHEJ <- lambda_I_NHEJ_new

        if (delta < tol) {
            if (verbose) cat("tolleranza raggiunta in iter = ", iter, "\n")
            break
        }
    }

    stats$P_M <- P_M
    stats$P_NHEJ <- P_NHEJ
    stats$L_M <- L_M
    stats$L_NHEJ <- L_NHEJ

    return(list(
        results = stats,
        pi_M = pi_M,
        pi_NHEJ = pi_NHEJ,
        lambda_I_M = lambda_I_M,
        lambda_I_NHEJ = lambda_I_NHEJ,
        history = history,
        n_iter = iter,
        converged = delta < tol
    ))
}

# FUNZIONE PRINCIPALE (precalcolo + EM)

run_em_M_vs_NHEJ_multi <- function(data, lambda_M, beta_M,
                              pi_M_init = 0.5,
                              lambda_I_M_init = 0.3, lambda_I_NHEJ_init = 0.3,
                              max_iter = 100, tol = 1e-8,
                              L0 = 2, verbose = TRUE) {

    cat("Pre-calcolo delle quantità fisse per variante (posizione, lunghezza)...\n")
    stats <- precompute_variant_stats_multi(
        data = data,
        lambda_M = lambda_M,
        beta_M = beta_M,
        lambda_I_M = lambda_I_M_init,      # servono solo come placeholder per calcolare L_M/L_NHEJ iniziali diagnostici
        lambda_I_NHEJ = lambda_I_NHEJ_init,
        L0 = L0
    )
    cat("Pre-calcolo completato:", nrow(stats), "variant_id.\n\n")

    cat("Esecuzione EM su pi_M, pi_NHEJ, lambda_I_M, lambda_I_NHEJ...\n")
    em_out <- em_prior_lambda_M_NHEJ_multi(
        stats = stats,
        pi_M_init = pi_M_init,
        lambda_I_M_init = lambda_I_M_init,
        lambda_I_NHEJ_init = lambda_I_NHEJ_init,
        max_iter = max_iter,
        tol = tol,
        verbose = verbose
    )

    return(em_out)
}

# CLI (opzionale, eseguito solo se lo script è chiamato da Rscript)
if (sys.nframe() == 0) {

    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) < 4) {
        stop("Uso: Rscript em_M_vs_NHEJ.R -i input.tsv -o output.tsv [-p pi_M_init] [-t tol] [-m max_iter]")
    }

    get_arg <- function(flag, default = NULL, numeric = FALSE) {
        pos <- which(args == flag)
        if (length(pos) == 0) return(default)
        val <- args[pos + 1]
        if (numeric) val <- as.numeric(val)
        val
    }

    input_file <- get_arg("-i")
    output_file <- get_arg("-o")
    pi_M_init <- get_arg("-p", default = 0.5, numeric = TRUE)
    tol <- get_arg("-t", default = 1e-8, numeric = TRUE)
    max_iter <- get_arg("-m", default = 100, numeric = TRUE)

    if (is.null(input_file) || is.null(output_file)) {
        stop("Devi specificare sia -i che -o.")
    }

    cat("Input :", input_file, "\n")
    cat("Output:", output_file, "\n\n")

    cat("Importazione del TSV...\n")
    data <- fread(input_file, sep = "\t", header = TRUE)
    cat("Importazione completata:", nrow(data), "righe x", ncol(data), "colonne\n\n")

    em_out <- run_em_M_vs_NHEJ_multi(
        data = data,
        lambda_M = 50,
        beta_M = 0.7,
        lambda_I_M_init = 0.3,
        lambda_I_NHEJ_init = 0.3,
        pi_M_init = pi_M_init,
        max_iter = max_iter,
        tol = tol
    )

    cat("\nPrior finali: pi_M =", em_out$pi_M, " pi_NHEJ =", em_out$pi_NHEJ, "\n")
    cat("Iterazioni:", em_out$n_iter, " Convergenza:", em_out$converged, "\n\n")

    cat("Scrittura dei risultati...\n")
    fwrite(em_out$results, file = output_file, sep = "\t", quote = FALSE, na = "NA")
    cat("Output salvato in:", output_file, "\n")
}

