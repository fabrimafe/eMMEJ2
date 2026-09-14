library(data.table)


# 1. MOTIF STRENGTH

motif_strength <- function(distance,motif_l,lambda,beta,L0 = 2) {

    motif_l = motif_l/2
    stopifnot(length(distance) == length(motif_l))
    stopifnot(lambda > 0)

    # Peso di ogni singolo motivo
    w <- exp(
        beta * (motif_l - L0) -
        abs(distance) / lambda
    )

    # Somma dei pesi dei motivi
    g <- sum(w)

    return(g)
}

# 2. POSITION PROBABILITY
# g_obs = motif strength nella posizione osservata
# g_cf  = motif strength nelle 100 counterfactual

position_probability <- function(g_obs, g_cf) {

    total <- g_obs + sum(g_cf)

    # se il motivo non viene trovato allora trattiamola come 0,
    # non come NaN.
    if (is.na(total) || total == 0) {
        return(0)
    }
    g_obs / total
}

# 3. PROBABILITÀ basandosi sulla LUNGHEZZA DELL'INDEL
# dgeom() in R parte da 0, mentre la lunghezza
# dell'indel parte da 1.

indel_length_prob <- function(L, lambda_I) {

    dgeom(
        L - 1,
        prob = lambda_I
    )
}

calculate_n_positions <- function(group) {
	# group = sottoinsieme di data relativo a un singolo variant_id
	if (!"g_M" %in% colnames(group)) {
		stop("La colonna non è presente nel gruppo.")
	} 
	n_valid <- sum(!is.na(group$g_M))
	
	return(n_valid)
} 

# 4. PARSING DELLE POSIZIONI DEI MOTIVI
# "40,50,70,100"
# diventa:
# c(40, 50, 70, 100)

parse_motif_positions <- function(x) {

    if (is.na(x) || x == "" || x == "-") {

        return(numeric(0))
    }

    positions <- as.numeric(unlist(strsplit(as.character(x),",")))

    positions <- positions[
        !is.na(positions)
    ]

    return(positions)
}

# 5. CALCOLO DELLA MOTIF STRENGTH DI UNA SINGOLA RIGA
# M = SD_inverted_deletion
# SD_ID_motif_pos
#       -> distance
#
# SD_ID_repeat_pat_len
#       -> length del motivo

calculate_M_strength <- function(row,lambda,beta,L0 = 2) {

    # SD_inverted_deletion is TRUE?
    # trimws() rimuove eventuali spazi introdotti da apply() quando converte l'intero data.frame in matrice di caratteri

    is_M <- trimws(row[["SD_inverted_deletion"]])
    #is_M <- row[["SD_inverted_deletion"]]

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

    # Se ci sono 3 motivi e la lunghezza è 4:
    # distance = c(40,70,100)
    # length   = c(4,4,4)

    motif_l <- rep(motif_lengths,length(distance))
   
    g <- motif_strength(distance = distance,motif_l = motif_l,lambda = lambda,beta = beta,L0 = L0)

    return(g)
}

# CALCOLO DEL RAPPORTO OBSERVED / COUNTERFACTUAL, g_obs / mean(g_cf)
obs_cf_ratio <- function(g_obs,g_cf) {
	cf_mean <- mean(g_cf,na.rm = TRUE)
	if (is.na(g_obs) || is.na(cf_mean) || cf_mean == 0) {
		return(NA_real_)
	}
	ratio <- g_obs / cf_mean
	return(ratio)
}

#Z-SCORE
motif_z_score <- function(g_obs, g_cf) {
	cf_mean <- mean(
		g_cf,
		na.rm = TRUE
	) 
	cf_sd <- sd(g_cf,na.rm = TRUE)
	if (is.na(g_obs) || is.na(cf_mean) || is.na(cf_sd) || cf_sd == 0) {
		return(NA_real_)
	}
	z <- (g_obs - cf_mean) / cf_sd
	return(z)
}



# LIKELIHOOD DI M
# L_M = P(indel length | M) * P(position | M)

likelihood_M <- function(L,lambda_I_M,g_obs,g_cf) {

    # Probabilità della lunghezza

    p_length_M <- indel_length_prob(L = L,lambda_I = lambda_I_M)

    p_position_M <- position_probability(g_obs = g_obs,g_cf = g_cf)


    # Likelihood complessiva di M

    L_M <- p_length_M * p_position_M


    return(L_M)
}


# 7. LIKELIHOOD DI NHEJ

# Per NHEJ, nel modello attuale, consideriamo soltanto
# la distribuzione della lunghezza dell'indel.

likelihood_NHEJ <- function(L,lambda_I_NHEJ, n_positions) {
	p_length_NHEJ <- indel_length_prob(L = L, lambda_I = lambda_I_NHEJ)
	p_position_NHEJ <- 1 / n_positions
	L_NHEJ <- p_length_NHEJ * p_position_NHEJ
	return(L_NHEJ)

}

# 8. CONVERSIONE DELLE LIKELIHOOD IN PROBABILITÀ
posterior_M_NHEJ <- function(L_M,L_NHEJ) {

    total <- L_M + L_NHEJ

    # Copre sia denominator == 0 sia NA/NaN (es. propagati da
    # indel_length mancante o da altri calcoli falliti)
    if (is.na(total) || total == 0) {

        return(list(P_M = NA,P_NHEJ = NA))
    }

    P_M <- L_M / total
    P_NHEJ <- L_NHEJ / total

    return(list(P_M = P_M,P_NHEJ = P_NHEJ))
}
# 9. FUNZIONE PRINCIPALE
calculate_M_vs_NHEJ <- function(data,lambda_M,beta_M,lambda_I_M,lambda_I_NHEJ,L0 = 2) {
   
    required_cols <- c("variant_id","observed","ANC","DER","SD_inverted_deletion","SD_ID_motif_pos","SD_ID_repeat_pat_len")

    missing_cols <- setdiff(required_cols,colnames(data))


    if (length(missing_cols) > 0) {

        stop(paste("Mancano le colonne:",paste(missing_cols,collapse = ", ")))
    }
     
    data$g_M <- apply(data,1,calculate_M_strength,lambda = lambda_M,beta = beta_M,L0 = L0)
       # Lista dei variant_id

    variant_ids <- unique(
        data$variant_id
    )

    results <- vector(
        "list",
        length(variant_ids)
    )

    # CICLO SUI VARIANT_ID

    for (i in seq_along(variant_ids)) {

        current_id <- variant_ids[i]

        # Selezioniamo le righe della variante
        group <- data[
            data$variant_id == current_id,
        ]
        # Riga osservata
        n_positions <- calculate_n_positions(group)

	if (n_positions != 101) {
		warning(
			paste(
			      "variant_id",current_id,"ha",n_positions,"posizioni invece di 101."
			)
		)
	} 


	obs <- group[
            group$observed == 1,
        ]


        if (nrow(obs) != 1) {

            warning(
                paste(
                    "variant_id",
                    current_id,
                    "non ha esattamente una riga observed = 1"
                )
            )

            next
        }

        # Counterfactual
       
        cf <- group[group$observed == 0,]

        # g_obs
       
        g_obs <- obs$g_M

        # g_cf

        g_cf <- cf$g_M

        # Lunghezza dell'indel
      
        indel_length <- (nchar(obs$ANC) - nchar(obs$DER))
       
        # Media delle counterfactual
	cf_mean <- mean(
		g_cf,
		na.rm = TRUE
	) 

	# Somma delle counterfactual
	cf_sum <- sum(
		g_cf,
		na.rm = TRUE
	)

	# SD delle counterfactual

	cf_sd <- sd(
		g_cf,
		na.rm = TRUE
	)

	# Rapporto observed / counterfactual mean
	ratio <- obs_cf_ratio(
		g_obs = g_obs,
		g_cf = g_cf
	)

	# Z-score
	z <- motif_z_score(
		g_obs = g_obs,
		g_cf = g_cf
	)

	p_position_M <- position_probability(g_obs = g_obs,g_cf = g_cf)

        # P(length | M)

        p_length_M <- indel_length_prob(L = indel_length,lambda_I = lambda_I_M)

        # Likelihood M

        L_M <- likelihood_M(L = indel_length,lambda_I_M = lambda_I_M,g_obs = g_obs,g_cf = g_cf)

        # Likelihood NHEJ

        L_NHEJ <- likelihood_NHEJ(L = indel_length,lambda_I_NHEJ = lambda_I_NHEJ,n_positions = n_positions)

        # Posterior M vs NHEJ Assumiamo prior uguali:
        # P(M) = P(NHEJ)
        
        posterior <- posterior_M_NHEJ(L_M = L_M,L_NHEJ = L_NHEJ)

        #risultato

        results[[i]] <- data.frame(variant_id = current_id,
				   indel_length = indel_length,
				   g_obs = g_obs,
				   g_cf_sum = sum(g_cf),
				   n_counterfactuals = length(g_cf),
				   cf_mean = cf_mean,
				   obs_cf_ratio = ratio,
				   cf_sd = cf_sd,
				   z_score = z,
				   likelihood_M = L_M,
				   likelihood_NHEJ = L_NHEJ,
				   P_M = posterior$P_M,
				   P_NHEJ = posterior$P_NHEJ
        )
    }

    # Uniamo i risultati
    results <- do.call(
        rbind,
        results
    )

    return(results)
}

# ARGOMENTI

args <- commandArgs(
    trailingOnly = TRUE
)

if (length(args) != 4) {

    stop(
        "Uso: Rscript script.R -i input.tsv -o output.tsv"
    )
}

# Posizione di -i

i_pos <- which(
    args == "-i"
)

# Posizione di -o

o_pos <- which(
    args == "-o"
)


# Controlliamo che esistano entrambi

if (length(i_pos) != 1 ||
    length(o_pos) != 1) {

    stop(
        "Devi specificare sia -i che -o."
    )
}


# Controlliamo che dopo -i e -o ci sia il percorso

if (i_pos == length(args) ||
    o_pos == length(args)) {

    stop(
        "-i e -o devono essere seguiti dal percorso del file."
    )
}


# Recuperiamo i percorsi

input_file <- args[
    i_pos + 1
]

output_file <- args[
    o_pos + 1
]

# 11. IMPORTAZIONE TSV

cat("Input :",input_file,"\n")

cat( "Output:",output_file,"\n\n")

cat("Importazione del TSV...\n")

data <- fread(
    input_file,
    sep = "\t",
    header = TRUE
)

cat(
    "Importazione completata:",
    nrow(data),
    "righe x",
    ncol(data),
    "colonne\n\n"
)

# 12. CALCOLO

cat("Calcolo M vs NHEJ...\n")

results <- calculate_M_vs_NHEJ(data = data,lambda_M = 50,beta_M = 0.7,lambda_I_M = 0.3,lambda_I_NHEJ = 0.3)

cat("Calcolo completato:",nrow(results),"variant_id analizzati.\n\n")
# 13. OUTPUT

cat("Scrittura dei risultati...\n")


fwrite(results,file = output_file,sep = "\t",quote = FALSE,na = "NA")

cat("Output salvato in:",output_file,"\n")

