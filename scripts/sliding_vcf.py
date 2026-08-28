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
