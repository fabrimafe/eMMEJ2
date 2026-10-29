"""
DSB simulator - versione "weighted motif choice".

Differenza rispetto alla versione precedente:
  prima  -> il meccanismo veniva simulato sul PRIMO motivo trovato nella finestra (str.find)
  adesso -> vengono enumerati TUTTI i motivi (o le configurazioni di motivi) validi nella
            finestra, a ciascuno viene assegnata una forza
                 exp( beta*(motif_l - L0) - |distance|/lambda )
            (stessa formula di motif_strength() usata nell'EM), le forze vengono normalizzate
            a 1 e il motivo usato per simulare viene ESTRATTO con quelle probabilita'.

Usare --selection first per tornare al comportamento "primo motivo trovato".
"""
import math
import random
import argparse
from pysam import FastaFile
import pandas as pd

def seq2complement(xseq):
    """Calculate complement of a sequence and NOT reverse complement"""
    x_complement = ""
    for i in xseq:
        if i == "T" or i == "t":
            x = "A"
        elif i == "A" or i == "a":
            x = "T"
        elif i == "G" or i == "g":
            x = "C"
        elif i == "C" or i == "c":
            x = "G"
        else:
            x = "N"
        x_complement = x_complement + x
    return(x_complement)

def seq2revcomplement(xseq):
    """Calculate reverse complement of a DNA sequence"""
    return(seq2complement(xseq)[::-1])

    #normalize vcf adding the nt anchor in substitution
def normalize_vcf(anc, der):
    
    while len(anc) > 1 and len(der) > 1 and anc[-1] == der[-1]:
        anc = anc[:-1]
        der = der[:-1]
    
    
    while len(anc) > 1 and len(der) > 1 and anc[0] == der[0]:
        anc = anc[1:]
        der = der[1:]
    
    return anc, der
    
    #removing nt when ANC and DER had the same nt at the beginning or at the end, and calculate the suffix removed,
    #the prefix removed and the total nt removed 
def sliding_vcf(anc, der):
    NT_removed_suffix = 0
    NT_removed_prefix = 0

    while len(anc) > 1 and len(der) > 1 and anc[-1] == der[-1]:
        anc = anc[:-1]
        der = der[:-1]
        NT_removed_suffix += 1

    while len(anc) > 1 and len(der) > 1 and anc[0] == der[0]:
        anc = anc[1:]
        der = der[1:]
        NT_removed_prefix += 1

    total_NT_removed = NT_removed_prefix + NT_removed_suffix 
    
    return anc, der, NT_removed_suffix, NT_removed_prefix, total_NT_removed
def pol_slip_del(refFA,chrom,pos,MH_length,max_length):
    down_MH = refFA.fetch(chrom,pos - MH_length,pos).upper()
    
    up_MH = refFA.fetch(chrom, pos, pos + MH_length).upper()
    
    if down_MH != up_MH:
        return('error')
    else:
        #print("uguali")
        ancestral_state_vcf = refFA.fetch(chrom,pos - 1,pos + len(up_MH)).upper()
        #print(ancestral_state_vcf)
        derived_state_vcf = refFA.fetch(chrom,pos - 1, pos).upper()
        #print(derived_state_vcf)
        
        variant_id=str(chrom)+"_"+str(pos)
        distance_motif = 0 
        indel_seq=refFA.fetch(chrom,pos,pos + len(up_MH)).upper()
        indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
        if (indelNs>-3):
            print("N found in indel")
            return('error')
        else:
            return(chrom,pos,variant_id,ancestral_state_vcf,derived_state_vcf,distance_motif)


def pol_slip_ins(refFA,chrom,pos,MH_length,max_length):
    down_MH = refFA.fetch(chrom,pos - MH_length,pos).upper()
    
    up_MH = refFA.fetch(chrom, pos, pos + MH_length).upper()
    window_fabio = refFA.fetch (chrom, pos - 10, pos + 10)
    print(window_fabio)
    if down_MH != up_MH:
        return('error')
    else:
        #print("uguali")
        derived_state_vcf = refFA.fetch(chrom,pos - 1,pos + len(up_MH)).upper()
        #print(ancestral_state_vcf)
        ancestral_state_vcf = refFA.fetch(chrom,pos - 1, pos).upper()
        #print(derived_state_vcf)
        
        variant_id=str(chrom)+"_"+str(pos)
        distance_motif = 0
        indel_seq=refFA.fetch(chrom,pos,pos + len(up_MH)).upper()

        indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
        if (indelNs>-3):
            print("N found in indel")
            return('error')
        else:
            return(chrom,pos,variant_id,ancestral_state_vcf,derived_state_vcf,distance_motif)


def fa2deletion(refFA,chrom,pos,min_length,max_length):
    """Function to create deletions on a fasta file giving outuput in vcf format. First and last 2kbp of chromosomes should not be used.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chrom: the chromosome to be mutated
    pos: the start (0-based) position of the mutation, i.e. if you want to delete position 121 and 122 you need to set 120
    min_length: minimum length of indel
    max_length: minimum length of indel
    Output has indel_length in addition to vcf-like annotation.
    """
    #indel_length=random.sample(range(min_length,max_length+1),1)[0]

    indel_length=max_length
    ancestral_state_vcf=refFA.fetch(chrom,pos-1,pos+indel_length).upper()
    derived_state_vcf=refFA.fetch(chrom,pos-1,pos).upper()
    
    indel_seq = refFA.fetch(chrom,pos,pos+indel_length).upper()
   
    variant_id=str(chrom)+"_"+str(pos)
    distance_motif = 0
    indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
    if (indelNs>-3):
        print("N found in indel")
        return('error')
    else:
        return(chrom,pos,variant_id,ancestral_state_vcf,derived_state_vcf,distance_motif)


def casual_sequence(lunghezza, seed=None):
    """Genera una sequenza genomica casuale composta da A, C, G, T."""
    if seed is not None:
        random.seed(seed)
    return "".join(random.choices("ACGT", k=lunghezza))
        
def random_substitution(refFA,chrom,pos,min_length,max_length):
    """Function to create substitutions on a fasta file giving outuput in vcf format. First and last 2kbp of chromosomes should not be used.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chrom: the chromosome to be mutated
    pos: the start (0-based) position of the mutation, i.e. if you want to replace position 121 and 122 you need to set 120
    min_length: minimum length of indel
    max_length: minimum length of indel
    Output has indel_length in addition to vcf-like annotation.
    """
    

    X_length=random.sample(range(2,max_length+1),1)[0]
    X_seq=indel_seq=refFA.fetch(chrom,pos,pos+X_length).upper()
   
    Z_length=random.sample(range(2,max_length+1),1)[0]
    #Y_length=random.sample(range(100,100+max_length+1),1)[0]
    #Z_seq=indel_seq=refFA.fetch(chrom,pos+X_length+Y_length,pos+X_length+Y_length+Z_length).upper()
    Z_seq=casual_sequence(Z_length)
    print("X_length: "+str(X_length))
    #print("Y_length: "+str(Y_length))
    print("Z_length: "+str(Z_length))

    seq_window=refFA.fetch(chrom,pos,pos+100)
    #print("seq_window: "+seq_window)
    
    distance_motif = 0
    #removing of the same nt at the beginning and at the end of ANC and DER, 
    #the number of the nt removed at the start of ANC and DER is used to slide the pos forward
    anc_slide, der_slide, NT_removed_suffix, NT_removed_prefix, total_NT_removed = sliding_vcf (X_seq,Z_seq)
    print("anc_slide: "+anc_slide)
    print("der_slide: "+der_slide)
    print("NT_removed_suffix: "+ str(NT_removed_suffix))
    print("NT_removed_prefix: " +str(NT_removed_prefix))
    pos_slide=pos+NT_removed_prefix
    variant_id=str(chrom)+"_"+str(pos_slide)
    print(variant_id)  
    if anc_slide == der_slide:
        return('error')
    else:
        indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
        if (indelNs>-3):
            print("N found in indel")
            return('error')
        else:
            return(chrom,pos_slide,variant_id,anc_slide, der_slide,distance_motif)
    

def fa2insertion(refFA,chrom,pos,min_length,max_length):
    """Function to create insertions on a fasta file giving outuput in vcf format. First and last 2kbp of chromosomes should not be used.
    Insertion sequences are generated by transposition, copying a region in the nearby 2000bp.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the mutation, i.e. if you want to delete position 121 and 122 you need to set 120
    min_length: minimum length of indel
    max_length: minimum length of indel
    """
    #indel_length=random.sample(range(min_length,max_length+1),1)[0]
    indel_length=max_length
    # --- taking the source_pos of the ref in a the pos surrounding window 
    window_start =  pos - 200      # non andare sotto 0 (inizio cromosoma)
    window_end = pos + 200

    # deleting indel_length from window_end range valido [window_start, window_end - indel_length]
    max_start = window_end - indel_length
    source_pos = random.randint(window_start, max_start)
    
    print("source_pos: " + str(source_pos))
    # get the random seq from the random source_pos 
    indel_seq = refFA.fetch(chrom, source_pos, source_pos + indel_length)
   
    ancestral_state_vcf=refFA.fetch(chrom,pos-1,pos).upper()
    derived_state_vcf=(refFA.fetch(chrom,pos-1,pos)+indel_seq).upper()
    variant_id=str(chrom)+"_"+str(pos)
    distance_motif=0
    indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
    if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
        print("N or NULL found in indel")
        return('error')
    else:
        return(chrom,pos,variant_id,ancestral_state_vcf,derived_state_vcf,distance_motif)




# =====================================================================================
#  FORZA DEI MOTIVI E SCELTA PESATA
# =====================================================================================

def log_motif_strength(distance, motif_l, lam, beta, L0=2):
    """log della forza di UN motivo (stessa formula della motif_strength() in R, senza la somma)."""
    return beta * (motif_l - L0) - abs(distance) / lam


def motif_strength(distance, motif_l, lam, beta, L0=2):
    """Forza di un motivo: exp(beta*(motif_l-L0) - |distance|/lambda)."""
    return math.exp(log_motif_strength(distance, motif_l, lam, beta, L0))


def choose_candidate(cands, W):
    """cands: lista di dict, ognuno con la chiave 'motifs' = lista di (distanza, lunghezza_motivo).
    La forza di una configurazione e' il PRODOTTO delle forze dei suoi motivi
    (calcolato come somma dei log, per evitare underflow).
    W: dict con lam, beta, L0, selection.
    Ritorna (candidato scelto, n candidati, probabilita' del candidato scelto)."""
    logw = [sum(log_motif_strength(d, l, W["lam"], W["beta"], W["L0"]) for d, l in c["motifs"])
            for c in cands]
    m = max(logw)
    w = [math.exp(x - m) for x in logw]       # softmax stabile
    tot = sum(w)
    probs = [x / tot for x in w]
    if W["selection"] == "first":
        idx = 0                               # comportamento vecchio: il primo trovato
    else:
        idx = random.choices(range(len(cands)), weights=probs, k=1)[0]
    return cands[idx], len(cands), probs[idx]


def F(refFA, chrom, start, end):
    """fetch in maiuscolo; stringa vuota per intervalli vuoti/invertiti (invece di un errore pysam)."""
    start = max(start, 0)
    if end <= start:
        return ""
    return refFA.fetch(chrom, start, end).upper()


def find_all(window, pattern):
    """Tutte le posizioni (anche sovrapposte) di pattern in window."""
    hits = []
    i = window.find(pattern)
    while i != -1:
        hits.append(i)
        i = window.find(pattern, i + 1)
    return hits


def has_N(seq):
    return any(c in seq for c in "Nn-")


def annot(refFA, chrom, bounds):
    """Stringa 'a|b|c|...' tagliata sui confini dati (solo per il log)."""
    return "|".join(F(refFA, chrom, bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1))


def report(c, n, p):
    print("candidates in window: " + str(n) + " | chosen: " +
          ", ".join(k + "=" + str(v) for k, v in c.items() if k != "motifs") +
          " | P(chosen)=" + format(p, ".4f"))


#  MECCANISMI CON MOTIVI (versione pesata)
#  Distanza usata per la forza = spaziatore (nt) tra il motivo e il suo "ancoraggio";
#  la colonna distance_motif in output resta quella di prima.


def fa2deletion_MMEJ(refFA, chrom, pos, MH_length, max_distance_MMEJ, W):
    """MMEJ: candidati = tutte le copie della MH entro la finestra. distanza = spaziatore taglio->MH2."""
    MH_seq = F(refFA, chrom, pos - MH_length, pos)
    seq_window = F(refFA, chrom, pos, pos + MH_length + max_distance_MMEJ)
    cands = [{"i_MH": i, "motifs": [(i, MH_length)]} for i in find_all(seq_window, MH_seq)]
    if not cands:
        print("no MH possible in range")
        return "error"
    c, n, p = choose_candidate(cands, W)
    report(c, n, p)
    i_MH = c["i_MH"]
    indel_length = i_MH + MH_length
    indel_seq = F(refFA, chrom, pos, pos + indel_length)
    print("MH motif: " + MH_seq)
    print("MH|DEL|DEL MH: ")
    print(annot(refFA, chrom, [pos - MH_length, pos, pos + i_MH, pos + indel_length]))
    ancestral_state_vcf = F(refFA, chrom, pos - 1, pos + indel_length)
    derived_state_vcf = F(refFA, chrom, pos - 1, pos)
    variant_id = str(chrom) + "_" + str(pos)
    if has_N(indel_seq):
        print("N found in indel")
        return "error"
    return (chrom, pos, variant_id, ancestral_state_vcf, derived_state_vcf, 0, n, p)


def fa2insertion_snapback(refFA, chrom, pos, MH_length, SD_motif_length, max_distance, W):
    """SD_inverted_insertion. Candidati = tutte le coppie (SDrevc, MHrevc) valide.
    Motivi pesati: SDrevc (distanza = Y, spaziatore del loop), MHrevc (distanza = Z, lunghezza inserzione)."""
    SD = SD_motif_length
    SD_motif = F(refFA, chrom, pos, pos + SD)
    MH_seq = F(refFA, chrom, pos - MH_length, pos)
    SD_revc = seq2revcomplement(SD_motif)
    MH_revc = seq2revcomplement(MH_seq)
    win1 = F(refFA, chrom, pos + SD + 1, pos + max_distance - 1)
    cands = []
    for k in find_all(win1, SD_revc):
        Y = k + 1                                   # 1 nt minimo per il loop
        win2 = F(refFA, chrom, pos + 2 * SD + Y, pos + max_distance)
        for Z in find_all(win2, MH_revc):
            if Z < 1:
                continue
            cands.append({"Y": Y, "Z": Z, "motifs": [(Y, SD), (Z, MH_length)]})
    if not cands:
        print("no SD/MH possible in range")
        return "error"
    c, n, p = choose_candidate(cands, W)
    report(c, n, p)
    Y, Z = c["Y"], c["Z"]
    print("MH|SD|snapback-loop|SDrevc|INSrevc|MHrevc: ")
    print(annot(refFA, chrom, [pos - MH_length, pos, pos + SD, pos + SD + Y, pos + 2 * SD + Y,
                               pos + 2 * SD + Y + Z, pos + 2 * SD + Y + Z + MH_length]))
    indel_seq = seq2revcomplement(F(refFA, chrom, pos + 2 * SD + Y, pos + 2 * SD + Y + Z))
    anchor = F(refFA, chrom, pos - 1, pos)
    ancestral_state_vcf = anchor
    derived_state_vcf = anchor + indel_seq
    variant_id = str(chrom) + "_" + str(pos)
    if has_N(indel_seq) or Z == 0 or anchor in ("N", "-"):
        print("N found in indel")
        return "error"
    distance_motif = Y + len(MH_seq)
    return (chrom, pos, variant_id, ancestral_state_vcf, derived_state_vcf, distance_motif, n, p)


def fa2deletion_snapback(refFA, chrom, pos, MH_length, SD_motif_length, max_distance, X_length, W):
    """SD_inverted_deletion (X fissa). Candidati = tutte le SDrevc con MHrevc subito dopo.
    Motivi pesati: SDrevc (distanza = Y), MHrevc (distanza = 0, adiacente)."""
    SD = SD_motif_length
    SD_motif = F(refFA, chrom, pos + X_length, pos + X_length + SD)
    MH_seq = F(refFA, chrom, pos - MH_length, pos)
    X_seq = F(refFA, chrom, pos, pos + X_length)
    SD_revc = seq2revcomplement(SD_motif)
    MH_revc = seq2revcomplement(MH_seq)
    win1 = F(refFA, chrom, pos + X_length + SD + 1, pos + max_distance - 1)
    cands = []
    for k in find_all(win1, SD_revc):
        Y = k + 1
        s = pos + X_length + 2 * SD + Y
        if F(refFA, chrom, s, s + MH_length) == MH_revc:
            cands.append({"Y": Y, "motifs": [(Y, SD), (0, MH_length)]})
    if not cands:
        print("no SD/MH possible in range")
        return "error"
    c, n, p = choose_candidate(cands, W)
    report(c, n, p)
    Y = c["Y"]
    print("MH|DEL X|SD|snapback-loop|SDrevc|MHrevc: ")
    print(annot(refFA, chrom, [pos - MH_length, pos, pos + X_length, pos + X_length + SD,
                               pos + X_length + SD + Y, pos + X_length + 2 * SD + Y,
                               pos + X_length + 2 * SD + Y + MH_length]))
    anchor = F(refFA, chrom, pos - 1, pos)
    derived_state_vcf = anchor
    ancestral_state_vcf = anchor + X_seq
    variant_id = str(chrom) + "_" + str(pos)
    if has_N(X_seq) or X_length == 0 or anchor in ("N", "-"):
        print("N found in indel")
        return "error"
    distance_motif = Y + len(MH_seq) + X_length
    return (chrom, pos, variant_id, ancestral_state_vcf, derived_state_vcf, distance_motif, n, p)


def fa2SD_direct_insertion(refFA, chrom, pos, MH_length, SD_motif_length, max_distance, W):
    """SD_direct_insertion. Candidati = tutte le coppie (MH2, SD2) valide.
    Motivi pesati: MH2 (distanza = Y, spaziatore dopo SD1), SD2 (distanza = Z, lunghezza inserzione)."""
    SD = SD_motif_length
    MH_seq = F(refFA, chrom, pos - MH_length, pos)
    SD_motif = F(refFA, chrom, pos, pos + SD)
    win1 = F(refFA, chrom, pos + SD, pos + max_distance - 1)
    cands = []
    for Y in find_all(win1, MH_seq):
        win2 = F(refFA, chrom, pos + SD + Y + MH_length, pos + max_distance)
        for Z in find_all(win2, SD_motif):
            if Z < 1:
                continue
            cands.append({"Y": Y, "Z": Z, "motifs": [(Y, MH_length), (Z, SD)]})
    if not cands:
        print("no MH/SD possible in range")
        return "error"
    c, n, p = choose_candidate(cands, W)
    report(c, n, p)
    Y, Z = c["Y"], c["Z"]
    a = pos + SD + Y + MH_length
    print("MH|SD|spacer|MH|INS|SD: ")
    print(annot(refFA, chrom, [pos - MH_length, pos, pos + SD, pos + SD + Y, a, a + Z, a + Z + SD]))
    indel_seq = F(refFA, chrom, a, a + Z)
    anchor = F(refFA, chrom, pos - 1, pos)
    ancestral_state_vcf = anchor
    derived_state_vcf = anchor + indel_seq
    variant_id = str(chrom) + "_" + str(pos)
    if has_N(indel_seq) or Z == 0 or anchor in ("N", "-"):
        print("N found in indel")
        return "error"
    distance_motif = Y + len(MH_seq)
    return (chrom, pos, variant_id, ancestral_state_vcf, derived_state_vcf, distance_motif, n, p)


def fa2SD_direct_deletion(refFA, chrom, pos, MH_length, SD_motif_length, max_distance, X_length, W):
    """SD_direct_deletion (X fissa). Candidati = tutte le MH2 con SD2 subito dopo.
    Motivi pesati: MH2 (distanza = spaziatore dopo SD1), SD2 (distanza = 0, adiacente)."""
    SD = SD_motif_length
    MH_seq = F(refFA, chrom, pos - MH_length, pos)
    SD_motif = F(refFA, chrom, pos + X_length, pos + X_length + SD)
    X_seq = F(refFA, chrom, pos, pos + X_length)
    win1 = F(refFA, chrom, pos + X_length + SD, pos + max_distance - 1)
    cands = []
    for i_MH in find_all(win1, MH_seq):
        s = pos + X_length + SD + i_MH + MH_length
        if F(refFA, chrom, s, s + SD) == SD_motif:
            cands.append({"i_MH": i_MH, "motifs": [(i_MH, MH_length), (0, SD)]})
    if not cands:
        print("no MH/SD possible in range")
        return "error"
    c, n, p = choose_candidate(cands, W)
    report(c, n, p)
    i_MH = c["i_MH"]
    print("MH|DEL X|SD|spacer|MH|SD: ")
    print(annot(refFA, chrom, [pos - MH_length, pos, pos + X_length, pos + X_length + SD,
                               pos + X_length + SD + i_MH, pos + X_length + SD + i_MH + MH_length,
                               pos + X_length + 2 * SD + i_MH + MH_length]))
    anchor = F(refFA, chrom, pos - 1, pos)
    derived_state_vcf = anchor
    ancestral_state_vcf = anchor + X_seq
    variant_id = str(chrom) + "_" + str(pos)
    if has_N(X_seq) or X_length == 0 or anchor in ("N", "-"):
        print("N found in indel")
        return "error"
    distance_motif = i_MH + len(X_seq) + len(MH_seq)
    return (chrom, pos, variant_id, ancestral_state_vcf, derived_state_vcf, distance_motif, n, p)


def fa2SD_direct_substitution(refFA, chrom, pos, MH_length, SD_motif_length, max_distance, X_length, W):
    """SD_direct_substitution. Candidati = tutte le coppie (MH2, SD2) valide; Z sostituisce X.
    Motivi pesati: MH2 (distanza = Y), SD2 (distanza = Z)."""
    SD = SD_motif_length
    MH_seq = F(refFA, chrom, pos - MH_length, pos)
    SD_motif = F(refFA, chrom, pos + X_length, pos + X_length + SD)
    win1 = F(refFA, chrom, pos + X_length + SD, pos + X_length + SD + max_distance - 1)
    cands = []
    for Y in find_all(win1, MH_seq):
        s = pos + X_length + SD + Y + MH_length
        win2 = F(refFA, chrom, s, s + max_distance)
        for Z in find_all(win2, SD_motif):
            if Z < 1:
                continue
            cands.append({"Y": Y, "Z": Z, "motifs": [(Y, MH_length), (Z, SD)]})
    if not cands:
        print("no MH/SD possible in range")
        return "error"
    c, n, p = choose_candidate(cands, W)
    report(c, n, p)
    Y, Z = c["Y"], c["Z"]
    s = pos + X_length + SD + Y + MH_length
    print("MH|X|SD|spacer|MH|Z|SD: ")
    print(annot(refFA, chrom, [pos - MH_length, pos, pos + X_length, pos + X_length + SD,
                               pos + X_length + SD + Y, s, s + Z, s + Z + SD]))
    indel_seq = F(refFA, chrom, s, s + Z)
    X_ANC_seq = F(refFA, chrom, pos, pos + X_length)
    anc_slide, der_slide, NT_removed_suffix, NT_removed_prefix, total_NT_removed = sliding_vcf(X_ANC_seq, indel_seq)
    if anc_slide == der_slide:
        return "error"
    pos_slided = pos + NT_removed_prefix
    variant_id = str(chrom) + "_" + str(pos_slided)
    if has_N(indel_seq) or Z == 0 or anc_slide in ("N", "-"):
        print("N found in indel")
        return "error"
    distance_motif = len(MH_seq) + Y + X_length
    return (chrom, pos_slided, variant_id, anc_slide, der_slide, distance_motif, n, p)


def fa2SD_snapback_substitution(refFA, chrom, pos, MH_length, SD_motif_length, max_distance, X_length, W):
    """SD_inverted_substitution. Candidati = tutte le coppie (SDrevc, MHrevc) valide; revcomp(Z) sostituisce X.
    Motivi pesati: SDrevc (distanza = Y), MHrevc (distanza = Z)."""
    SD = SD_motif_length
    MH_seq = F(refFA, chrom, pos - MH_length, pos)
    SD_motif = F(refFA, chrom, pos + X_length, pos + X_length + SD)
    SD_revc = seq2revcomplement(SD_motif)
    MH_revc = seq2revcomplement(MH_seq)
    win1 = F(refFA, chrom, pos + X_length + SD + 1, pos + X_length + SD + 1 + max_distance - 1)
    cands = []
    for k in find_all(win1, SD_revc):
        Y = k + 1
        s = pos + X_length + 2 * SD + Y
        win2 = F(refFA, chrom, s, s + max_distance)
        for Z in find_all(win2, MH_revc):
            if Z < 1:
                continue
            cands.append({"Y": Y, "Z": Z, "motifs": [(Y, SD), (Z, MH_length)]})
    if not cands:
        print("no SD/MH possible in range")
        return "error"
    c, n, p = choose_candidate(cands, W)
    report(c, n, p)
    Y, Z = c["Y"], c["Z"]
    s = pos + X_length + 2 * SD + Y
    print("MH|X|SD|spacer|SDrevc|Zrevc|MHrevc: ")
    print(annot(refFA, chrom, [pos - MH_length, pos, pos + X_length, pos + X_length + SD,
                               pos + X_length + SD + Y, s, s + Z, s + Z + MH_length]))
    Z_seq = F(refFA, chrom, s, s + Z)
    indel_seq = seq2revcomplement(Z_seq)
    X_ANC_seq = F(refFA, chrom, pos, pos + X_length)
    nt_anchor = F(refFA, chrom, pos - 1, pos)
    anc, der = normalize_vcf(X_ANC_seq, indel_seq)
    anc_slide, der_slide, NT_removed_suffix, NT_removed_prefix, total_NT_removed = sliding_vcf(X_ANC_seq, indel_seq)
    if anc_slide == der_slide:
        print("ANC and DER are equal")
        return "error"
    pos_slided = pos + NT_removed_prefix
    variant_id = str(chrom) + "_" + str(pos_slided)
    ancestral_state_vcf = nt_anchor + anc
    derived_state_vcf = nt_anchor + der
    if has_N(indel_seq) or Z == 0 or ancestral_state_vcf in ("N", "-"):
        print("N found in indel")
        return "error"
    distance_motif = len(MH_seq) + Y + X_length
    return (chrom, pos_slided, variant_id, anc_slide, der_slide, distance_motif, n, p)


def parse_args():
    all_args = argparse.ArgumentParser()
    all_args.add_argument("-r", "--ref", required=True,
        help="reference genome in fasta format")
    all_args.add_argument("-o", "--out", required=False,
        help="output file")
    all_args.add_argument("-c", "--chrom", required=False,
        help="chrom of position of DSB; necessary of nsims=1, otherwise ignored")
    all_args.add_argument("-p", "--pos", required=False, type=int,
        help="position of DSB; necessary if nsims=1; otherwise ignored")
    all_args.add_argument("-n", "--nsims", required=False, default=100,
        help="number of random DSB to simulate")
    all_args.add_argument("-m", "--mechanism", required=True, default="deletion_MMEJ",
        help="mechanism of DSB-repair. One of: deletion_NHEJ,deletion_MMEJ,pol_slip_del,pol_slip_ins,"
             "SD_inverted_insertion,SD_inverted_deletion,SD_inverted_substitution,"
             "SD_direct_insertion,SD_direct_deletion,SD_direct_substitution,random_insertion,random_substitution")
    all_args.add_argument("-l", "--MHlength", required=False, type=int, default=2,
        help="length of microhomology motif ( or minimum indel size ) ")
    all_args.add_argument("-s", "--SDlength", required=False, type=int, default=5,
        help="length of SD motif/distance second MH in MMEJtrans")
    all_args.add_argument("-d", "--maxdistance", required=False, type=int, default=10,
        help="maximum distance of MH motif or maximum indel length (= window size)")
    all_args.add_argument("-zx", "--X_length", required=False, type=int, default=10,
        help="maximum distance between MH and SD")
    all_args.add_argument("-vr", "--verbose", action="store_true",
        help="if set, output has extra columns: distance_motif, n_candidates, p_chosen")
    # --- nuovi argomenti: scelta pesata del motivo ---
    all_args.add_argument("--selection", choices=["weighted", "first"], default="weighted",
        help="weighted: estrae il motivo con probabilita' proporzionale alla forza; "
             "first: primo motivo trovato (comportamento vecchio)")
    all_args.add_argument("--lam", type=float, default=1.0,
        help="lambda della forza: exp(-distance/lambda). Usare il valore stimato dall'EM (default 1 = exp(-d))")
    all_args.add_argument("--beta", type=float, default=0.0,
        help="beta della forza: exp(beta*(motif_l-L0)). Default 0")
    all_args.add_argument("--L0", type=float, default=2.0,
        help="lunghezza di riferimento L0 della forza (default 2)")
    all_args.add_argument("--seed", type=int, default=None,
        help="seed random per la riproducibilita'")
    return vars(all_args.parse_args())


def simulate_one(refFA, xmechanism, xchr, xpos, xMHlength, xSDlength, xmaxdistance, xX_length, W):
    if xmechanism == "deletion_MMEJ":
        return fa2deletion_MMEJ(refFA, xchr, xpos, xMHlength, xmaxdistance, W)
    elif xmechanism == "deletion_NHEJ":
        return fa2deletion(refFA, xchr, xpos, xMHlength, xmaxdistance)
    elif xmechanism == "pol_slip_del":
        return pol_slip_del(refFA, xchr, xpos, xMHlength, xmaxdistance)
    elif xmechanism == "pol_slip_ins":
        return pol_slip_ins(refFA, xchr, xpos, xMHlength, xmaxdistance)
    elif xmechanism == "SD_inverted_insertion":
        return fa2insertion_snapback(refFA, xchr, xpos, xMHlength, xSDlength, xmaxdistance, W)
    elif xmechanism == "SD_inverted_deletion":
        return fa2deletion_snapback(refFA, xchr, xpos, xMHlength, xSDlength, xmaxdistance, xX_length, W)
    elif xmechanism == "random_insertion":
        return fa2insertion(refFA, xchr, xpos, xMHlength, xmaxdistance)
    elif xmechanism == "random_substitution":
        return random_substitution(refFA, xchr, xpos, xMHlength, xmaxdistance)
    elif xmechanism == "SD_direct_insertion":
        return fa2SD_direct_insertion(refFA, xchr, xpos, xMHlength, xSDlength, xmaxdistance, W)
    elif xmechanism == "SD_direct_deletion":
        return fa2SD_direct_deletion(refFA, xchr, xpos, xMHlength, xSDlength, xmaxdistance, xX_length, W)
    elif xmechanism == "SD_direct_substitution":
        return fa2SD_direct_substitution(refFA, xchr, xpos, xMHlength, xSDlength, xmaxdistance, xX_length, W)
    elif xmechanism == "SD_inverted_substitution":
        return fa2SD_snapback_substitution(refFA, xchr, xpos, xMHlength, xSDlength, xmaxdistance, xX_length, W)
    raise ValueError("unknown mechanism: " + str(xmechanism))


def main():
    print("Parse arguments")
    args = parse_args()
    if args["seed"] is not None:
        random.seed(args["seed"])
    W = {"lam": args["lam"], "beta": args["beta"], "L0": args["L0"], "selection": args["selection"]}
    if W["lam"] <= 0:
        raise ValueError("--lam must be > 0")

    fastafile = args["ref"]
    outfile = args["out"]
    chrom = args["chrom"]
    pos = args["pos"]
    nsims = int(args["nsims"])
    xmechanism = args["mechanism"]
    xMHlength = args["MHlength"]
    xSDlength = args["SDlength"]
    xmaxdistance = args["maxdistance"]
    xX_length = args["X_length"]
    verbose = args["verbose"]

    refFA = FastaFile(fastafile)

    MIN_CONTIG_LEN = 10000   # soglia (>= 2*2000 + margine per il meccanismo)
    valid_contigs = [(c, l) for c, l in zip(refFA.references, refFA.lengths) if l >= MIN_CONTIG_LEN]
    if not valid_contigs:
        raise ValueError(f"no contig had >= {MIN_CONTIG_LEN} bp in refa")
    valid_names = [c for c, l in valid_contigs]
    valid_lengths = [l for c, l in valid_contigs]
    valid_len_dict = dict(valid_contigs)

    counter_sims = 0
    myvcf = list()
    seen_ids = set()

    while len(myvcf) < nsims:
        if chrom is None:
            xchr = random.choices(valid_names, weights=valid_lengths, k=1)[0]
        else:
            xchr = chrom
        if pos is None:
            xpos = random.randint(2000, valid_len_dict[xchr] - 2000 - 1)
        else:
            xpos = int(pos)
        print(xchr)
        print(xpos)
        res = simulate_one(refFA, xmechanism, xchr, xpos, xMHlength, xSDlength, xmaxdistance, xX_length, W)
        counter_sims = counter_sims + 1

        if res != "error":
            if len(res) == 6:                       # meccanismi senza scelta di motivi
                res = res + (None, None)
            vcf_id = str(res[0]) + "_" + str(res[1])
            if vcf_id in seen_ids:
                print("duplicate simulated vcf, skipping: " + vcf_id)
                continue
            seen_ids.add(vcf_id)
            myvcf.append(res)
            print(res)
        else:
            print("error")
            continue

    print("N of attempted simulations:" + str(counter_sims))
    cols = ['#CHR', 'POS', 'variant_id', 'REF', 'ALT']
    extra = ['distance_motif', 'n_candidates', 'p_chosen']

    if verbose:
        myvcf = pd.DataFrame(myvcf, columns=cols + extra)
    else:
        myvcf = pd.DataFrame([r[:5] for r in myvcf], columns=cols)

    myvcf['ANC'] = 0

    if verbose:
        myvcf = myvcf[cols + ['ANC'] + extra]

    if outfile is not None:
        myvcf.to_csv(outfile, sep="\t", index=None)


if __name__ == "__main__":
    main()

