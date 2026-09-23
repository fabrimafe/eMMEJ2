
library(data.table)
library(optimx)

# =========================================================
# 1. FUNZIONI BASE

# calcola la forza dei motivi presenti i nuna determinata posizione e calcola la loro forza complessiva,
#questa foza dipende dalla lunghezza del motivo, distanza dall'indel,
#lambda che controlla quanto diminuiscel'effetto sulla distanza, e beta che controlla l'effetto della lunghezza del motivo 
#controlla anche che ogni motivo abbia una lunghezza associata

motif_strength <- function(distance,motif_l,lambda,beta,L0=2){
    stopifnot(length(distance)==length(motif_l),lambda>0)
    sum(exp(beta*(motif_l-L0)-abs(distance)/lambda))
}
#trasforma la forza dei motivi in una probabilità relativa della posizione osservata, rispetto alla forza nei counterfactuals
position_probability <- function(g_obs,g_cf){
    total<-g_obs+sum(g_cf)
    if(is.na(total)||total==0) return(0)
    g_obs/total
}
#la probabilità di un'indel legata solamente alla lunghezza dell'indel
indel_length_prob <- function(L,lambda_I) dgeom(L-1,prob=lambda_I)
#conta le posizioni valide nel gruppo
calculate_n_positions <- function(group){
    if(!"g_M"%in%colnames(group)) stop("La colonna non è presente nel gruppo.")
    sum(!is.na(group$g_M))
}
#fa il parsing delle distanca, insomma pulisce le distanze
parse_motif_positions <- function(x){
    if(is.na(x)||x==""||x=="-") return(numeric(0))
    p<-as.numeric(unlist(strsplit(as.character(x),",")))
    p[!is.na(p)]
}
#controlla se M è presente, e se è presente calcola la forza del motivo, in base alla distanza del motivo e alla lunghezza
calculate_M_strength <- function(row,lambda,beta,L0=2){
    is_M<-trimws(row[["SD_inverted_deletion"]])
    if(is.na(is_M)||is_M!="TRUE") return(0)
    distance<-parse_motif_positions(row[["SD_ID_motif_pos"]])
    if(length(distance)==0) return(0)
    motif_l<-as.numeric(row[["SD_ID_repeat_pat_len"]])
    if(is.na(motif_l)) return(0)
    motif_strength(distance,rep(motif_l,length(distance)),lambda,beta,L0)
}
#non la usiamo più
obs_cf_ratio <- function(g_obs,g_cf){
    cf_mean<-mean(g_cf,na.rm=TRUE)
    if(is.na(g_obs)||is.na(cf_mean)||cf_mean==0) return(NA_real_)
    g_obs/cf_mean
}
#valuta quanto è distante il motivo osservato dalla distribuzione dei cf
motif_z_score <- function(g_obs,g_cf){
    cf_mean<-mean(g_cf,na.rm=TRUE)
    cf_sd<-sd(g_cf,na.rm=TRUE)
    if(is.na(g_obs)||is.na(cf_mean)||is.na(cf_sd)||cf_sd==0) return(NA_real_)
    (g_obs-cf_mean)/cf_sd
}
# calcola la L di M, in base alla lunghezza dell'indel e a quanto è forte il motivo rispetto alla finestra attorno 
likelihood_M <- function(L,lambda_I_M,g_obs,g_cf){
    indel_length_prob(L,lambda_I_M)*position_probability(g_obs,g_cf)
}
#L_NHEJ quanto è probabile che NHEJ avvenga lì solo in base alla lunghezza indel
likelihood_NHEJ <- function(L,lambda_I_NHEJ,n_positions){
    indel_length_prob(L,lambda_I_NHEJ)*(1/n_positions)
}
#calcola posterior e prior ma sta operazione viene fatta dentro em_prior_lambda_M_NHEJ_multi()
posterior_M_NHEJ_prior <- function(L_M,L_NHEJ,pi_M,pi_NHEJ){
    num_M<-pi_M*L_M
    num_NHEJ<-pi_NHEJ*L_NHEJ
    total<-num_M+num_NHEJ
    if(is.na(total)||total==0)
        return(list(P_M=NA_real_,P_NHEJ=NA_real_))
    list(P_M=num_M/total,P_NHEJ=num_NHEJ/total)
}

# =========================================================
# 2. PRE-PARSING MOTIVI
#is_M vettore logico che contiene true or false in base se trova il meccanismo 
#crea dei valori ovvero l'indice, il numero della variante
#la distanza
#la lunghezza del motivo
#per ogni riga M estrae lunghezza motivo e estrae la distanza
prepare_motif_data <- function(data){
    is_M<-trimws(as.character(data$SD_inverted_deletion))=="TRUE"
    idx<-integer(0); distance<-numeric(0); motif_l<-numeric(0)
    for(i in which(is_M)){
        d<-parse_motif_positions(data$SD_ID_motif_pos[i])
        ml<-as.numeric(data$SD_ID_repeat_pat_len[i])
        if(length(d)==0||is.na(ml)) next
        idx<-c(idx,rep(i,length(d)))
        distance<-c(distance,d)
        motif_l<-c(motif_l,rep(ml,length(d)))
    }
    list(is_M=is_M,motif_row_idx=idx,motif_distance=distance,motif_length=motif_l)
}
#estrae una lista in cui abbiamo indice, distanza e la lunghezza del motif
# =========================================================
# 3. CALCOLO VELOCE g_M
#calcolare g_M per tutte le righe del dataset senza analizzare nuovamente le stringhe.
#prende tutti i mootivi già preparati e calcola i pesi in una volta e somma per riga
#inizialmente crea g vuoto poi calcola contemporaneamente tutti i pesi senza fare un ciclo R esplicito

calculate_g_all_fast <- function(motif_info,n_data,lambda,beta,L0=2){
    g<-numeric(n_data)
    if(length(motif_info$motif_row_idx)==0) return(g)
    w<-exp(beta*(motif_info$motif_length-L0)-abs(motif_info$motif_distance)/lambda)
    s<-rowsum(w,motif_info$motif_row_idx,reorder=FALSE)
    idx<-as.integer(rownames(s))
    g[idx]<-s[,1]
    g
}
# 4. PRE-CALCOLO
#prende le colonne, vede se tutte le colonne ci sono e prepara i motivi
# calculate_g_all_fast calcola g_M e aggiunge il g_M dataset per ogni riga 
#ad ogni riga assegna il suo indice originale
#crea la lista della variante e poi processa una variante alla volta
#separa alignments = 0 e quando != 0
#g_cf contiene le forze dei motivi nelle posizioni counterfactual.contiene le forze dei motivi nelle posizioni counterfactual
# n_position conta le posizioni valide
precompute_variant_stats_multi <- function(data,lambda_M,beta_M,lambda_I_M,lambda_I_NHEJ,L0=2){
    required_cols<-c("variant_id","alignment_ID","ANC","DER",
                     "SD_inverted_deletion","SD_ID_motif_pos","SD_ID_repeat_pat_len")
    missing_cols<-setdiff(required_cols,colnames(data))
    if(length(missing_cols)>0)
        stop(paste("Mancano le colonne:",paste(missing_cols,collapse=", ")))

    motif_info<-prepare_motif_data(data)
    data$g_M<-calculate_g_all_fast(motif_info,nrow(data),lambda_M,beta_M,L0)
    data$.data_idx<-seq_len(nrow(data))

    variant_ids<-unique(data$variant_id)
    out_list<-vector("list",length(variant_ids))

    for(i in seq_along(variant_ids)){
        id<-variant_ids[i]
        group<-data[variant_id==id]
        cf<-group[alignment_ID==0]
        g_cf<-cf$g_M
        n_positions<-sum(!is.na(group$g_M))
        align_rows<-group[alignment_ID!=0]

        if(nrow(align_rows)==0){
            warning(paste("variant_id",id,"non ha alignment candidati."))
            next
        }

        align_list<-vector("list",nrow(align_rows))

        for(j in seq_len(nrow(align_rows))){
            row<-align_rows[j]
            g_obs<-row$g_M
            L<-nchar(row$ANC)-nchar(row$DER)
            p_pos<-position_probability(g_obs,g_cf)
            LM<-likelihood_M(L,lambda_I_M,g_obs,g_cf)
            LN<-likelihood_NHEJ(L,lambda_I_NHEJ,n_positions)

            align_list[[j]]<-data.frame(
                variant_id=id,
                alignment_ID=row$alignment_ID,
                indel_length=L,
                g_obs=g_obs,
                n_counterfactuals=length(g_cf),
                n_positions=n_positions,
                p_position_M=p_pos,
                L_M=LM,
                L_NHEJ=LN,
                data_idx=row$.data_idx
            )
        }
        out_list[[i]]<-do.call(rbind,align_list)
    }
    #per ongi allinemanto candidato, prende la forza del motivo nella posizione osservata e  poi calcola la lunghezza dell'indel
    #la p_pos che è la position probability
    #LM e LN
    #per ogni alignment viene salvata una riga con varie informazioni
    stats<-do.call(rbind,out_list)
    variant_levels<-unique(stats$variant_id)
    stats$variant_group<-match(stats$variant_id,variant_levels)
    #ogni ID delle varianti diventa un numero 
    cf_index_list<-lapply(
        variant_levels,
        function(id)data[variant_id==id&alignment_ID==0,.data_idx]
    )

    list(
        stats=stats,
        motif_info=motif_info,
        cf_index_list=cf_index_list,
        variant_levels=variant_levels
    )
}
# =========================================================
# 5. OTTIMIZZAZIONE lambda_M E beta_M
#a qusto punto ottimizziamo lambda e beta in modo che lmbda e beta massimizza Q ma non aggiorna le responsibilities

optimize_lambda_beta_M <- function(data,stats,P_M_old,motif_info,
                                   cf_index_list,lambda_M_init,beta_M_init,L0=2){

    n_data<-nrow(data)
    obs_data_idx<-stats$data_idx
    variant_group<-stats$variant_group
    #massimizzare Q andando a minimalizzare -Q, prova diversi paramente di lambda e beta
    Q_function<-function(par){
        lambda_M<-par[1]
        beta_M<-par[2]

        if(!is.finite(lambda_M)||!is.finite(beta_M)||
           lambda_M<=0||beta_M<0) return(1e100)
	#per ogni parametro calcola tutti i g_M e viene rpetuto da optimx
        g_all<-calculate_g_all_fast(
            motif_info,n_data,lambda_M,beta_M,L0
        )
	#g_obs prende direttamente il g della posizione osservata di ogni alignment.
        g_obs<-g_all[obs_data_idx]
	#aggiunge ad ogni riga di stats la somma Cf della sua variant
	#cf_sum_by_variant calocla la forza dei cf per ongi variant e assegna
        cf_sum_by_variant<-vapply(
            cf_index_list,
            function(idx){
                if(length(idx)==0) 0 else sum(g_all[idx])
            },
            numeric(1)
        )

        g_cf_sum<-cf_sum_by_variant[variant_group]
        total<-g_obs+g_cf_sum
        p_position_M<-numeric(nrow(stats))

        valid_total<-is.finite(g_obs)&is.finite(total)&total>0
        p_position_M[valid_total]<-g_obs[valid_total]/total[valid_total]
	#elimina i vlaori in cui abbiamo zero perchè dopo dobbiamo fare il log
        valid<-P_M_old>0&
               is.finite(P_M_old)&
               p_position_M>0&
               is.finite(p_position_M)

        if(!any(valid)) return(1e100)
	#funzione Q per lambda e beta e ogni osservazione viene pesata con le vecchie responsibilities

        Q<-sum(P_M_old[valid]*log(p_position_M[valid]))
        -Q
    }
    #cerca il massimo di Q minimizzando  -Q
    opt<-optimx(
        par=c(lambda_M=lambda_M_init,beta_M=beta_M_init),
        fn=Q_function,
        method="L-BFGS-B",
        lower=c(lambda_M=1e-6,beta_M=0),
        upper=c(lambda_M=1e4,beta_M=10),
        control=list(trace=0)
    )
    #quando prende il risultato migliore allora restituisce lambda_M, beta_m, Q, convergence;optimx results
    best<-opt[1,]

    list(
        lambda_M=best$lambda_M,
        beta_M=best$beta_M,
        Q=-best$value,
        convergence=best$convcode,
        optimx_result=opt
    )
}
# =========================================================
# 6. EM
#parti da valori iniziali e poi li assegna come nuovi
EM_4V2.R
em_prior_lambda_M_NHEJ_multi <- function(
    data,stats,motif_info,cf_index_list,
    pi_M_init=0.5,
    lambda_M_init=50,
    beta_M_init=0.7,
    lambda_I_M_init=0.3,
    lambda_I_NHEJ_init=0.3,
    max_iter=100,
    tol=1e-5,
    verbose=TRUE
){

    stopifnot(all(c(
        "variant_id","p_position_M",
        "n_positions","indel_length"
    )%in%colnames(stats)))

    stopifnot(pi_M_init>0&&pi_M_init<1)
    stopifnot(lambda_I_M_init>0&&lambda_I_M_init<1)
    stopifnot(lambda_I_NHEJ_init>0&&lambda_I_NHEJ_init<1)
    stopifnot(lambda_M_init>0)

    pi_M<-pi_M_init
    pi_NHEJ<-1-pi_M_init
    lambda_M<-lambda_M_init
    beta_M<-beta_M_init
    lambda_I_M<-lambda_I_M_init
    lambda_I_NHEJ<-lambda_I_NHEJ_init

    p_position_NHEJ<-1/stats$n_positions
    L_vec<-stats$indel_length
    v_id<-stats$variant_id
    first_of_variant<-!duplicated(v_id)
    #salva le lunghezze e gli ID
    history<-data.frame(
        iter=integer(0),
        pi_M=numeric(0),
        pi_NHEJ=numeric(0),
        lambda_M=numeric(0),
        beta_M=numeric(0),
        lambda_I_M=numeric(0),
        lambda_I_NHEJ=numeric(0),
        delta=numeric(0),
        Loglik=numeric(0)
    )

    eps<-1e-6

    for(iter in seq_len(max_iter)){

        # -------------------------------------------------
        # E-STEP: g_M
        # usa gli attuali lambda e beta per calcolare la forza dei motivi 

        g_all<-calculate_g_all_fast(
            motif_info,
            nrow(data),
            lambda_M,
            beta_M,
            2
        )

        # -------------------------------------------------
        # E-STEP: p_position_M
        # prende la posizione osservata e recupera i valori di g per ogni observed

        g_obs<-g_all[stats$data_idx]
	#calcola la forza dei motivi delle cf
        cf_sum_by_variant<-vapply(
            cf_index_list,
            function(idx){
                if(length(idx)==0) 0 else sum(g_all[idx])
            },
            numeric(1)
        )
	#e assegna la somma alla riga corrispondente
        g_cf_sum<-cf_sum_by_variant[
            stats$variant_group
        ]
	#p_position_M e questa è la componente motivo/posizione
        total<-g_obs+g_cf_sum
        p_position_M<-numeric(nrow(stats))

        valid_total<-is.finite(g_obs)&
                     is.finite(total)&
                     total>0

        p_position_M[valid_total]<-
            g_obs[valid_total]/total[valid_total]

        # -------------------------------------------------
        # LIKELIHOOD
        # -------------------------------------------------

        p_length_M<-dgeom(
            L_vec-1,
            prob=lambda_I_M
        )

        p_length_NHEJ<-dgeom(
            L_vec-1,
            prob=lambda_I_NHEJ
        )

        L_M<-p_length_M*p_position_M
        L_NHEJ<-p_length_NHEJ*p_position_NHEJ
        
        # POSTERIOR
        # applico i prior e moltiplico per le likelihood

        num_M<-pi_M*L_M
        num_NHEJ<-pi_NHEJ*L_NHEJ
	#per una stessa variant_id allora somma i ocntributi delle diverse posizione candidate
	#perchè come abbiamo detto, in caso di diversi allineamenti dobbiamo scoprire quale
	#allineamento ha una probabilità maggiore che avvenga e queto dipende dalla forza del motivo
        variant_tot<-ave(
            num_M+num_NHEJ,
            v_id,
            FUN=sum
        )
	#e quindi per ogni allinemaneto calcola la responsibility di ogni allignment divisa facendo fratto la somma dei contributi
        P_M<-ifelse(
            is.na(variant_tot)|variant_tot==0,
            NA_real_,
            num_M/variant_tot
        )

        P_NHEJ<-ifelse(
            is.na(variant_tot)|variant_tot==0,
            NA_real_,
            num_NHEJ/variant_tot
        )
	#calcola la log_likelihood complessiva del dataset
        Loglik<-sum(
            log(variant_tot[first_of_variant]),
            na.rm=TRUE
        )
        # M-STEP

        # Ottimizzazione lambda_M e beta_M
        # i pesi che ho ottenuto nell'E-step diventano i pesi con cui ottimizzi i parametri di posizione di M infatti dopo ottengo i nuovi valori di lmbda e beta
        opt_M<-optimize_lambda_beta_M(
            data=data,
            stats=stats,
            P_M_old=P_M,
            motif_info=motif_info,
            cf_index_list=cf_index_list,
            lambda_M_init=lambda_M,
            beta_M_init=beta_M,
            L0=2
        )

        lambda_M_new<-opt_M$lambda_M
        beta_M_new<-opt_M$beta_M
        # Aggiornamento pi
        # aggiorno pi_M e sommo le responsibilities di M all'interno di ogni variante 

        pM_per_variant<-ave(
            P_M,
            v_id,
            FUN=sum
        )
	# aggiorno pi_M e sommo le responsibilities di NHEJ all'interno di ogni variante
        pNHEJ_per_variant<-ave(
            P_NHEJ,
            v_id,
            FUN=sum
        )
	#pi_M_new è la media tra i valori di M in ogni variante e quindi il nuovo prior è la frazione media di responabilità attribuita a M
        pi_M_new<-mean(
            pM_per_variant[first_of_variant],
            na.rm=TRUE
        )

        pi_NHEJ_new<-mean(
            pNHEJ_per_variant[first_of_variant],
            na.rm=TRUE
        )
        # Aggiornamento lambda_I
        # le lunghezze vengono pesate in base alle responsibilities M
	# le osservazioni che sembrano maggiormente M allora contibuiscono maggiormente alla stima delle lunghezze M

        lambda_I_M_new<-
            sum(P_M,na.rm=TRUE)/
            sum(P_M*L_vec,na.rm=TRUE)
	#per NHEJ, la distribuzione NHEJ viene imparata tramite le responsibilities
        lambda_I_NHEJ_new<-
            sum(P_NHEJ,na.rm=TRUE)/
            sum(P_NHEJ*L_vec,na.rm=TRUE)
	#due limiti di sicurezza
        lambda_I_M_new<-
            min(max(lambda_I_M_new,eps),1-eps)

        lambda_I_NHEJ_new<-
            min(max(lambda_I_NHEJ_new,eps),1-eps)

        # Convergenza
        #umi calcola la  convergenza che è utile poe poi stoppare le iterazioni 
        delta<-max(
            abs(pi_M_new-pi_M),
            abs(lambda_M_new-lambda_M),
            abs(beta_M_new-beta_M),
            abs(lambda_I_M_new-lambda_I_M),
            abs(lambda_I_NHEJ_new-lambda_I_NHEJ)
        )

        history<-rbind(
            history,
            data.frame(
                iter=iter,
                pi_M=pi_M_new,
                pi_NHEJ=pi_NHEJ_new,
                lambda_M=lambda_M_new,
                beta_M=beta_M_new,
                lambda_I_M=lambda_I_M_new,
                lambda_I_NHEJ=lambda_I_NHEJ_new,pi_M_new
                delta=delta,
                Loglik=Loglik
            )
        )
	#stampa l'iterazione
        if(verbose){
            cat(sprintf(
                "Iter %3d: pi_M = %.6f  lambda_M = %.6f  beta_M = %.6f  lambda_I_M = %.6f  lambda_I_NHEJ = %.6f  delta = %.2e  Loglik = %.6f\n",
                iter,
                pi_M_new,
                lambda_M_new,
                beta_M_new,
                lambda_I_M_new,
                lambda_I_NHEJ_new,
                delta,
                Loglik
            ))
        }
	#e ad ogni iterazione i nuovi valri diventano quelli utilizzati nell'iterazione successiva
        pi_M<-pi_M_new
        pi_NHEJ<-pi_NHEJ_new
        lambda_M<-lambda_M_new
        beta_M<-beta_M_new
        lambda_I_M<-lambda_I_M_new
        lambda_I_NHEJ<-lambda_I_NHEJ_new
	#controllo convergenza
        if(delta<tol){
            if(verbose)
                cat("tolleranza raggiunta in iter =",iter,"\n")
            break
        }
    }

    stats$P_M<-P_M
    stats$P_NHEJ<-P_NHEJ
    stats$L_M<-L_M
    stats$L_NHEJ<-L_NHEJ

    list(
        results=stats,
        pi_M=pi_M,
        pi_NHEJ=pi_NHEJ,
        lambda_M=lambda_M,
        beta_M=beta_M,
        lambda_I_M=lambda_I_M,
        lambda_I_NHEJ=lambda_I_NHEJ,
        history=history,
        n_iter=iter,
        converged=delta<tol
    )
}

# 7. FUNZIONE PRINCIPALE
# =========================================================

run_em_M_vs_NHEJ_multi <- function(
    data,
    lambda_M,
    beta_M,
    pi_M_init=0.5,
    lambda_I_M_init=0.3,
    lambda_I_NHEJ_init=0.3,
    max_iter=100,
    tol=1e-8,
    L0=2,
    verbose=TRUE
){

    cat("Pre-calcolo delle quantità fisse...\n")

    pre<-precompute_variant_stats_multi(
        data=data,
        lambda_M=lambda_M,
        beta_M=beta_M,
        lambda_I_M=lambda_I_M_init,
        lambda_I_NHEJ=lambda_I_NHEJ_init,
        L0=L0
    )

    stats<-pre$stats
    motif_info<-pre$motif_info
    cf_index_list<-pre$cf_index_list

    cat(
        "Pre-calcolo completato:",
        nrow(stats),
        "righe di alignment.\n\n"
    )

    cat(
        "Esecuzione EM su pi_M, pi_NHEJ, lambda_M, beta_M, lambda_I_M, lambda_I_NHEJ...\n"
    )

    em_out<-em_prior_lambda_M_NHEJ_multi(
        data=data,
        stats=stats,
        motif_info=motif_info,
        cf_index_list=cf_index_list,
        pi_M_init=pi_M_init,
        lambda_M_init=lambda_M,
        beta_M_init=beta_M,
        lambda_I_M_init=lambda_I_M_init,
        lambda_I_NHEJ_init=lambda_I_NHEJ_init,
        max_iter=max_iter,
        tol=tol,
        verbose=verbose
    )

    em_out
}

# 8. CLI
# =========================================================

if(sys.nframe()==0){

    args<-commandArgs(trailingOnly=TRUE)

    if(length(args)<4)
        stop("Uso: Rscript em_M_vs_NHEJ.R -i input.tsv -o output.tsv [-p pi_M_init] [-t tol] [-m max_iter]")

    get_arg<-function(flag,default=NULL,numeric=FALSE){
        pos<-which(args==flag)
        if(length(pos)==0) return(default)
        val<-args[pos+1]
        if(numeric) val<-as.numeric(val)
        val
    }

    input_file<-get_arg("-i")
    output_file<-get_arg("-o")
    pi_M_init<-get_arg("-p",0.5,TRUE)
    tol<-get_arg("-t",1e-8,TRUE)
    max_iter<-get_arg("-m",100,TRUE)

    if(is.null(input_file)||is.null(output_file))
        stop("Devi specificare sia -i che -o.")

    cat("Input :",input_file,"\n")
    cat("Output:",output_file,"\n\n")

    cat("Importazione del TSV...\n")
    data<-fread(
        input_file,
        sep="\t",
        header=TRUE
    )

    cat(
        "Importazione completata:",
        nrow(data),
        "righe x",
        ncol(data),
        "colonne\n\n"
    )

    em_out<-run_em_M_vs_NHEJ_multi(
        data=data,
        lambda_M=50,
        beta_M=0.7,
        lambda_I_M_init=0.3,
        lambda_I_NHEJ_init=0.3,
        pi_M_init=pi_M_init,
        max_iter=max_iter,
        tol=tol
    )

    cat(
        "\nPrior finali: pi_M =",
        em_out$pi_M,
        " pi_NHEJ =",
        em_out$pi_NHEJ,
        "\n"
    )

    cat(
        "lambda_M finale =",
        em_out$lambda_M,
        "\n"
    )

    cat(
        "beta_M finale =",
        em_out$beta_M,
        "\n"
    )

    cat(
        "lambda_I_M finale =",
        em_out$lambda_I_M,
        "\n"
    )

    cat(
        "lambda_I_NHEJ finale =",
        em_out$lambda_I_NHEJ,
        "\n"
    )

    cat(
        "Iterazioni:",
        em_out$n_iter,
        " Convergenza:",
        em_out$converged,
        "\n\n"
    )

    cat("Scrittura dei risultati...\n")

    fwrite(
        em_out$results,
        file=output_file,
        sep="\t",
        quote=FALSE,
        na="NA"
    )

    cat(
        "Output salvato in:",
        output_file,
        "\n"
    )
}


