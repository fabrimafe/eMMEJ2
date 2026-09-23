# =========================================================
# EM M vs NHEJ + optimx per lambda_M e beta_M
# =========================================================

library(data.table)
library(optimx)

# 1. MOTIF STRENGTH
motif_strength <- function(distance, motif_l, lambda, beta, L0 = 2) {
    stopifnot(length(distance) == length(motif_l))
    stopifnot(lambda > 0)
    w <- exp(beta * (motif_l - L0) - abs(distance) / lambda)
    return(sum(w))
}

# 2. POSITION PROBABILITY
position_probability <- function(g_obs, g_cf) {
    total <- g_obs + sum(g_cf)
    if (is.na(total) || total == 0) return(0)
    return(g_obs / total)
}

# 3. PROBABILITÀ LUNGHEZZA INDEL
indel_length_prob <- function(L, lambda_I) {
    return(dgeom(L - 1, prob = lambda_I))
}

# 4. NUMERO POSIZIONI VALIDE
calculate_n_positions <- function(group) {
    if (!"g_M" %in% colnames(group))
        stop("La colonna non è presente nel gruppo.")
    return(sum(!is.na(group$g_M)))
}

# 5. PARSING POSIZIONI MOTIVI
parse_motif_positions <- function(x) {
    if (is.na(x) || x == "" || x == "-") return(numeric(0))
    positions <- as.numeric(unlist(strsplit(as.character(x), ",")))
    return(positions[!is.na(positions)])
}

# 6. MOTIF STRENGTH DI UNA SINGOLA RIGA
calculate_M_strength <- function(row, lambda, beta, L0 = 2) {
    is_M <- trimws(row[["SD_inverted_deletion"]])
    if (is.na(is_M) || is_M != "TRUE") return(0)

    distance <- parse_motif_positions(row[["SD_ID_motif_pos"]])
    if (length(distance) == 0) return(0)

    motif_lengths <- as.numeric(row[["SD_ID_repeat_pat_len"]])
    if (is.na(motif_lengths)) return(0)

    motif_l <- rep(motif_lengths, length(distance))
    return(motif_strength(distance, motif_l, lambda, beta, L0))
}

# 7. RAPPORTO OBSERVED / COUNTERFACTUAL
obs_cf_ratio <- function(g_obs, g_cf) {
    cf_mean <- mean(g_cf, na.rm = TRUE)
    if (is.na(g_obs) || is.na(cf_mean) || cf_mean == 0)
        return(NA_real_)
    return(g_obs / cf_mean)
}

# 8. Z-SCORE
motif_z_score <- function(g_obs, g_cf) {
    cf_mean <- mean(g_cf, na.rm = TRUE)
    cf_sd <- sd(g_cf, na.rm = TRUE)
    if (is.na(g_obs) || is.na(cf_mean) || is.na(cf_sd) || cf_sd == 0)
        return(NA_real_)
    return((g_obs - cf_mean) / cf_sd)
}

# 9. LIKELIHOOD DI M
likelihood_M <- function(L, lambda_I_M, g_obs, g_cf) {
    p_length_M <- indel_length_prob(L, lambda_I_M)
    p_position_M <- position_probability(g_obs, g_cf)
    return(p_length_M * p_position_M)
}

# 10. LIKELIHOOD DI NHEJ
likelihood_NHEJ <- function(L, lambda_I_NHEJ, n_positions) {
    p_length_NHEJ <- indel_length_prob(L, lambda_I_NHEJ)
    p_position_NHEJ <- 1 / n_positions
    return(p_length_NHEJ * p_position_NHEJ)
}

# 11. POSTERIOR CON PRIOR ESPLICITO
posterior_M_NHEJ_prior <- function(L_M, L_NHEJ, pi_M, pi_NHEJ) {
    num_M <- pi_M * L_M
    num_NHEJ <- pi_NHEJ * L_NHEJ
    total <- num_M + num_NHEJ

    if (is.na(total) || total == 0)
        return(list(P_M = NA_real_, P_NHEJ = NA_real_))

    return(list(P_M = num_M / total, P_NHEJ = num_NHEJ / total))
}

# 12. PRE-CALCOLO MULTI-ALIGNMENT
precompute_variant_stats_multi <- function(data, lambda_M, beta_M,
                                           lambda_I_M, lambda_I_NHEJ, L0 = 2) {
    required_cols <- c("variant_id", "alignment_ID", "ANC", "DER",
                       "SD_inverted_deletion", "SD_ID_motif_pos",
                       "SD_ID_repeat_pat_len")
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0)
        stop(paste("Mancano le colonne:", paste(missing_cols, collapse = ", ")))

    data$g_M <- apply(data, 1, calculate_M_strength,
                      lambda = lambda_M, beta = beta_M, L0 = L0)

    variant_ids <- unique(data$variant_id)
    out_list <- vector("list", length(variant_ids))

    for (i in seq_along(variant_ids)) {
        current_id <- variant_ids[i]
        group <- data[data$variant_id == current_id, ]

        # Counterfactual pool
        cf <- group[group$alignment_ID == 0, ]
        g_cf <- cf$g_M

        # Tutte le righe del gruppo
        n_positions <- sum(!is.na(group$g_M))

        # Alignment candidati
        align_rows <- group[group$alignment_ID != 0, ]

        if (nrow(align_rows) == 0) {
            warning(paste("variant_id", current_id,
                          "non ha alignment candidati (alignment_ID != 0)."))
            next
        }

        align_list <- vector("list", nrow(align_rows))

        for (j in seq_len(nrow(align_rows))) {
            row <- align_rows[j, ]
            g_obs <- row$g_M
            indel_length <- nchar(row$ANC) - nchar(row$DER)

            p_position_M <- position_probability(g_obs, g_cf)
            L_M <- likelihood_M(indel_length, lambda_I_M, g_obs, g_cf)
            L_NHEJ <- likelihood_NHEJ(indel_length, lambda_I_NHEJ, n_positions)

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

    return(do.call(rbind, out_list))
}

# 13. OTTIMIZZAZIONE lambda_M E beta_M CON optimx

optimize_lambda_beta_M <- function(data, stats, P_M_old,
                                   lambda_M_init, beta_M_init, L0 = 2) {

    Q_function <- function(par) {

        lambda_M <- par[1]
        beta_M <- par[2]

        if (!is.finite(lambda_M) || !is.finite(beta_M) ||
            lambda_M <= 0 || beta_M < 0)
            return(1e100)

        # Ricalcolo g_M sulle righe ORIGINALI

        g_all <- vapply(seq_len(nrow(data)),function(i) {
                calculate_M_strength(data[i, ],lambda = lambda_M,beta = beta_M,L0 = L0
                )
            },
            numeric(1)
        )
        # Calcolo p_position_M per le righe di stats
        
        p_position_M <- numeric(nrow(stats))

        for (i in seq_len(nrow(stats))) {

            id <- stats$variant_id[i]
            alignment <- stats$alignment_ID[i]

            # tutte le righe del variant nel data originale
            idx_data <- which(data$variant_id == id)

            # counterfactual: alignment_ID == 0
            cf_idx <- idx_data[data$alignment_ID[idx_data] == 0]

            # alignment candidato corrispondente
            obs_idx <- idx_data[
                data$alignment_ID[idx_data] == alignment
            ]

            if (length(obs_idx) == 0 || length(cf_idx) == 0) {
                p_position_M[i] <- 0
                next
            }
            g_obs <- g_all[obs_idx[1]]
            g_cf <- g_all[cf_idx]
            total <- g_obs + sum(g_cf)
            if (is.na(total) || total <= 0 || is.na(g_obs)) {
                p_position_M[i] <- 0
            } else {
                p_position_M[i] <- g_obs / total
            }
        }

        # -------------------------------------------------
        # Q(lambda_M, beta_M)
        # -------------------------------------------------
        valid <- P_M_old > 0 & is.finite(P_M_old) &
                 p_position_M > 0 &
                 is.finite(p_position_M)

        if (!any(valid))
            return(1e100)

        Q <- sum(
            P_M_old[valid] *
            log(p_position_M[valid])
        )

        # optimx minimizza
        return(-Q)
    }

    opt <- optimx(
        par = c(lambda_M = lambda_M_init,beta_M = beta_M_init),
        fn = Q_function,method = "L-BFGS-B",
        lower = c(lambda_M = 1e-6,beta_M = 0
        ),
        upper = c(lambda_M = 1e4,beta_M = 10),
        control = list(trace = 0)
    )

    best <- opt[1, ]

    return(list(
        lambda_M = best$lambda_M,
        beta_M = best$beta_M,
        Q = -best$value,
        convergence = best$convcode,
        optimx_result = opt
    ))
}
# 14. EM
em_prior_lambda_M_NHEJ_multi <- function(data,stats,
                                          pi_M_init = 0.5,
                                          lambda_M_init = 50,
                                          beta_M_init = 0.7,
                                          lambda_I_M_init = 0.3,
                                          lambda_I_NHEJ_init = 0.3,
                                          max_iter = 100, tol = 1e-5,
                                          verbose = TRUE) {

    stopifnot(all(c("variant_id", "p_position_M",
                    "n_positions", "indel_length") %in% colnames(stats)))
    stopifnot(pi_M_init > 0 && pi_M_init < 1)
    stopifnot(lambda_I_M_init > 0 && lambda_I_M_init < 1)
    stopifnot(lambda_I_NHEJ_init > 0 && lambda_I_NHEJ_init < 1)
    stopifnot(lambda_M_init > 0)

    pi_M <- pi_M_init
    pi_NHEJ <- 1 - pi_M_init
    lambda_M <- lambda_M_init
    beta_M <- beta_M_init
    lambda_I_M <- lambda_I_M_init
    lambda_I_NHEJ <- lambda_I_NHEJ_init

    p_position_NHEJ <- 1 / stats$n_positions
    L_vec <- stats$indel_length
    v_id <- stats$variant_id
    first_of_variant <- !duplicated(v_id)

    history <- data.frame(
        iter = integer(0),
        pi_M = numeric(0),
        pi_NHEJ = numeric(0),
        lambda_M = numeric(0),
        beta_M = numeric(0),
        lambda_I_M = numeric(0),
        lambda_I_NHEJ = numeric(0),
        delta = numeric(0),
        Loglik = numeric(0)
    )

    eps <- 1e-6

    for (iter in seq_len(max_iter)) {

        # ---- E-step: ricalcolo g_M ----
        g_all <- vapply(seq_len(nrow(data)),
    function(i) {calculate_M_strength(data[i, ],lambda = lambda_M,beta = beta_M,L0 = 2
        )
    },
    numeric(1)
)

        # ---- E-step: p_position_M ----
        p_position_M <- numeric(nrow(stats))
        variant_ids <- unique(stats$variant_id)

        for (id in variant_ids) {
        idx_data  <- which(data$variant_id == id)
        cf_idx    <- idx_data[data$alignment_ID[idx_data] == 0]
        g_cf      <- g_all[cf_idx]

        idx_stats <- which(stats$variant_id == id)
        for (j in idx_stats) {
            align_id <- stats$alignment_ID[j]
            obs_idx  <- idx_data[data$alignment_ID[idx_data] == align_id]
            g_obs    <- g_all[obs_idx[1]]
            total    <- g_obs + sum(g_cf)

            if (is.na(total) || total <= 0 || is.na(g_obs))
                p_position_M[j] <- 0
            else
                p_position_M[j] <- g_obs / total
    }
}

        # ---- Likelihood ----
        p_length_M <- dgeom(L_vec - 1, prob = lambda_I_M)
        p_length_NHEJ <- dgeom(L_vec - 1, prob = lambda_I_NHEJ)

        L_M <- p_length_M * p_position_M
        L_NHEJ <- p_length_NHEJ * p_position_NHEJ

        # ---- E-step: posterior ----
        num_M <- pi_M * L_M
        num_NHEJ <- pi_NHEJ * L_NHEJ

        variant_tot <- ave(num_M + num_NHEJ, v_id, FUN = sum)

        P_M <- ifelse(
            is.na(variant_tot) | variant_tot == 0,
            NA_real_,
            num_M / variant_tot
        )

        P_NHEJ <- ifelse(
            is.na(variant_tot) | variant_tot == 0,
            NA_real_,
            num_NHEJ / variant_tot
        )

        Loglik <- sum(
            log(variant_tot[first_of_variant]),
            na.rm = TRUE
        )

        # =================================================
        # M-STEP
        # =================================================

        # ---- Ottimizzazione lambda_M e beta_M ----
        opt_M <- optimize_lambda_beta_M(
            data = data,
            stats = stats,
            P_M_old = P_M,
            lambda_M_init = lambda_M,
            beta_M_init = beta_M,
            L0 = 2
        )

        lambda_M_new <- opt_M$lambda_M
        beta_M_new <- opt_M$beta_M

        # ---- Aggiornamento pi ----
        pM_per_variant <- ave(P_M, v_id, FUN = sum)
        pNHEJ_per_variant <- ave(P_NHEJ, v_id, FUN = sum)

        pi_M_new <- mean(
            pM_per_variant[first_of_variant],
            na.rm = TRUE
        )

        pi_NHEJ_new <- mean(
            pNHEJ_per_variant[first_of_variant],
            na.rm = TRUE
        )

        # ---- Aggiornamento lambda_I ----
        lambda_I_M_new <- sum(P_M, na.rm = TRUE) /
                          sum(P_M * L_vec, na.rm = TRUE)

        lambda_I_NHEJ_new <- sum(P_NHEJ, na.rm = TRUE) /
                             sum(P_NHEJ * L_vec, na.rm = TRUE)

        lambda_I_M_new <- min(max(lambda_I_M_new, eps), 1 - eps)
        lambda_I_NHEJ_new <- min(max(lambda_I_NHEJ_new, eps), 1 - eps)

        # ---- Convergenza ----
        delta <- max(
            abs(pi_M_new - pi_M),
            abs(lambda_M_new - lambda_M),
            abs(beta_M_new - beta_M),
            abs(lambda_I_M_new - lambda_I_M),
            abs(lambda_I_NHEJ_new - lambda_I_NHEJ)
        )

        history <- rbind(
            history,
            data.frame(
                iter = iter,
                pi_M = pi_M_new,
                pi_NHEJ = pi_NHEJ_new,
                lambda_M = lambda_M_new,
                beta_M = beta_M_new,
                lambda_I_M = lambda_I_M_new,
                lambda_I_NHEJ = lambda_I_NHEJ_new,
                delta = delta,
                Loglik = Loglik
            )
        )

        if (verbose) {
            cat(sprintf(
                "Iter %3d: pi_M = %.6f  lambda_M = %.6f  beta_M = %.6f  lambda_I_M = %.6f  lambda_I_NHEJ = %.6f  delta = %.2e  Loglik = %.6f\n",
                iter, pi_M_new, lambda_M_new, beta_M_new,
                lambda_I_M_new, lambda_I_NHEJ_new, delta, Loglik
            ))
        }

        # ---- Aggiornamento parametri ----
        pi_M <- pi_M_new
        pi_NHEJ <- pi_NHEJ_new
        lambda_M <- lambda_M_new
        beta_M <- beta_M_new
        lambda_I_M <- lambda_I_M_new
        lambda_I_NHEJ <- lambda_I_NHEJ_new

        if (delta < tol) {
            if (verbose)
                cat("tolleranza raggiunta in iter =", iter, "\n")
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
        lambda_M = lambda_M,
        beta_M = beta_M,
        lambda_I_M = lambda_I_M,
        lambda_I_NHEJ = lambda_I_NHEJ,
        history = history,
        n_iter = iter,
        converged = delta < tol
    ))
}

# =========================================================
# 15. FUNZIONE PRINCIPALE
# =========================================================

run_em_M_vs_NHEJ_multi <- function(data, lambda_M, beta_M,
                                    pi_M_init = 0.5,
                                    lambda_I_M_init = 0.3,
                                    lambda_I_NHEJ_init = 0.3,
                                    max_iter = 100, tol = 1e-8,
                                    L0 = 2, verbose = TRUE) {

    cat("Pre-calcolo delle quantità fisse per variante (posizione, lunghezza)...\n")

    stats <- precompute_variant_stats_multi(
        data = data,
        lambda_M = lambda_M,
        beta_M = beta_M,
        lambda_I_M = lambda_I_M_init,
        lambda_I_NHEJ = lambda_I_NHEJ_init,
        L0 = L0
    )

    cat("Pre-calcolo completato:", nrow(stats), "righe di alignment.\n\n")

    cat("Esecuzione EM su pi_M, pi_NHEJ, lambda_M, beta_M, lambda_I_M, lambda_I_NHEJ...\n")

    em_out <- em_prior_lambda_M_NHEJ_multi(
        data = data,
        stats = stats,
        pi_M_init = pi_M_init,
        lambda_M_init = lambda_M,
        beta_M_init = beta_M,
        lambda_I_M_init = lambda_I_M_init,
        lambda_I_NHEJ_init = lambda_I_NHEJ_init,
        max_iter = max_iter,
        tol = tol,
        verbose = verbose
    )

    return(em_out)
}

# =========================================================
# 16. CLI
# =========================================================

if (sys.nframe() == 0) {

    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) < 4)
        stop("Uso: Rscript em_M_vs_NHEJ.R -i input.tsv -o output.tsv [-p pi_M_init] [-t tol] [-m max_iter]")

    get_arg <- function(flag, default = NULL, numeric = FALSE) {
        pos <- which(args == flag)
        if (length(pos) == 0) return(default)
        val <- args[pos + 1]
        if (numeric) val <- as.numeric(val)
        return(val)
    }

    input_file <- get_arg("-i")
    output_file <- get_arg("-o")
    pi_M_init <- get_arg("-p", default = 0.5, numeric = TRUE)
    tol <- get_arg("-t", default = 1e-8, numeric = TRUE)
    max_iter <- get_arg("-m", default = 100, numeric = TRUE)

    if (is.null(input_file) || is.null(output_file))
        stop("Devi specificare sia -i che -o.")

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

    cat("\nPrior finali: pi_M =", em_out$pi_M,
        " pi_NHEJ =", em_out$pi_NHEJ, "\n")
    cat("lambda_M finale =", em_out$lambda_M, "\n")
    cat("beta_M finale =", em_out$beta_M, "\n")
    cat("lambda_I_M finale =", em_out$lambda_I_M, "\n")
    cat("lambda_I_NHEJ finale =", em_out$lambda_I_NHEJ, "\n")
    cat("Iterazioni:", em_out$n_iter,
        " Convergenza:", em_out$converged, "\n\n")

    cat("Scrittura dei risultati...\n")

    fwrite(
        em_out$results,
        file = output_file,
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )

    cat("Output salvato in:", output_file, "\n")
}
