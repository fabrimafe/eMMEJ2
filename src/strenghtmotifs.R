#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)
})

# 
strengthmotif_1 <- function(dt, colonna, group_col = "variant_id") {
  setDT(dt)
  
  if (is.character(dt[[colonna]])) {
    dt[, (colonna) := as.logical(get(colonna))]
  }
  
  dt[, true_count := as.numeric(get(colonna))]
  
  risultato <- dt[, .(
    N_variants = .N,
    N_true = sum(true_count, na.rm = TRUE),
    N_true_variant_id = sum(true_count, na.rm = TRUE) / .N
  ), by = group_col]
  
  return(risultato)
}


strengthmotif_2 <- function(dt, colonna, group_col = "variant_id") {
  setDT(dt)
  
  
  dt[, n_motivi := sapply(get(colonna), function(x) {
    if (is.na(x) || trimws(x) == "") return(0)
    length(strsplit(x, ",")[[1]])
  })]
  
  # Raggruppa per group_col
  risultato <- dt[, .(
    N_variants = .N,
    N_motifs = sum(n_motivi, na.rm = TRUE),
    N_motifsxvariants = sum(n_motivi, na.rm = TRUE) / .N
  ), by = group_col]
  
  return(risultato)
}

#lapply prende la colonna che specifico io nel comando e applica una funzione per ogni valore x che trova,
#se non esiste allora assegna un valore vuoto "", se esiste fra str.split e divide per tutti i valori che 
#sono separati da una virgola, unlist alla fine, tutti i vettori (uno per riga) vengono uniti in un unico
#vettore lungo, con tutte le distanze di tutte le righe messe insieme (appiattite).

strengthmotif_3 <- function(dt, colonna, group_col = "variant_id") {
  setDT(dt)
  
  risultato <- dt[, {
    distances <- unlist(lapply(get(colonna), function(x) {
      if (is.na(x) || trimws(x) == "") return(numeric(0))
      as.numeric(strsplit(x, ",")[[1]])
    }))
    .(N_totale = .N,
      N_distances = length(distances),
      mean_distances = mean(distances, na.rm = TRUE))
  }, by = group_col]
  
  return(risultato)
}


#ARGOMENTI ---
#commandArgs() recupera tutti gli argomenti passati al comando quando lanci Rscript
#args è ora un vettore che tutte le flag e osa c'è dentro le flag
args <- commandArgs(trailingOnly = TRUE)
#qui costuisco opt che è la lista di opzioni, in cui se è NULL non vengono specificati e devo esserlo, 
#infatti dopo si controlla, group lo specifico io a priori, e se è FALSE allora è una flag booleana
parse_args_manual <- function(args) {
  opt <- list(input = NULL, output = NULL, colonna = NULL,
              group = "variant_id", motif = FALSE, numberM = FALSE, distance = FALSE)

#questo è un ciclo che va avanti di varie posizioni in base a che valore c'è, va da 1 alla lunghezza di args,
#se a (valore dentro args per iterazione i), contiene -i o --input allora prende il valore successivo che è quello che metto io nel comando,
#e poi va avanti di due passando a quello dopo, questo succede per i NULL
#per le flag booleano non abbiamo un valore da specificare dopo, quindi se c'è allora bene, se non c'è va avanti di uno e non di due
  
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-i", "--input")) {
      opt$input <- args[i + 1]; i <- i + 2
    } else if (a %in% c("-o", "--output")) {
      opt$output <- args[i + 1]; i <- i + 2
    } else if (a %in% c("-c", "--colonna")) {
      opt$colonna <- args[i + 1]; i <- i + 2
    } else if (a %in% c("-g", "--group")) {
      opt$group <- args[i + 1]; i <- i + 2
    } else if (a %in% c("-m", "--motif")) {
      opt$motif <- TRUE; i <- i + 1
    } else if (a %in% c("-n", "--numberM")) {
	opt$numberM <- TRUE; i <- i + 1
    } else if (a %in% c("-d", "--distance")) {
	opt$distance <- TRUE; i <- i + 1
    } else {
      stop(paste("Argomento non riconosciuto:", a), call. = FALSE)
    }
  }
  return(opt)
}

opt <- parse_args_manual(args)

#controlla se le colonne esistono e sono state inserite correttamente
#if non esiste opt$colonna x allora stoppa e printa la frase
if (is.null(opt$input)) {
  stop("input is missing  in -i / --input", call. = FALSE)
}
if (is.null(opt$output)) {
  stop("output is missing  in -o / --output", call. = FALSE)
}
if (is.null(opt$colonna)) {
  stop("column is missing in -c / --colonna", call. = FALSE)
}

# file input 
dt <- fread(opt$input, sep = "\t")

# ---default column ---, se el colonne non esistono allora stoppa
if (!(opt$colonna %in% names(dt))) {
  stop(paste0("coulmn '", opt$colonna, "' doesn't exist in input.\n",
              "columns: ", paste(names(dt), collapse = ", ")), call. = FALSE)
}
if (!(opt$group %in% names(dt))) {
  stop(paste0("group_by '", opt$group, "'doesn't exist in input.\n",
              "columns: ", paste(names(dt), collapse = ", ")), call. = FALSE)
}

# --- doing function ---
if (opt$motif) {
  cat("Eseguo strengthmotif_1 su colonna:", opt$colonna, "raggruppato per:", opt$group, "\n")
  risultato <- strengthmotif_1(dt, colonna = opt$colonna, group_col = opt$group)
  
  fwrite(risultato, opt$output, sep = "\t")
  cat("Saved in:", opt$output, "\n")
  
} else if (opt$numberM) {
  cat("count of motif in column:", opt$colonna, "raggruppato per:", opt$group, "\n")
  risultato <- strengthmotif_2(dt, colonna = opt$colonna, group_col = opt$group)
  
  fwrite(risultato, opt$output, sep = "\t")
  cat("Saved in:", opt$output, "\n")

} else if (opt$distance) {
  cat("mean distance in column: ", opt$colonna, "raggruppato per:", opt$group, "\n")
   risultato <- strengthmotif_3(dt, colonna = opt$colonna, group_col = opt$group)  
  
  fwrite(risultato, opt$output, sep = "\t")
  cat("Saved in:", opt$output, "\n")
  
} else {
  cat("No function activated. Use -m, -n ora - d .\n")
}


