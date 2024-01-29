import os
import sys
from Bio import SeqIO

"""
PF00698.25      Acyl_transf_1   General_PKS
PF02801.26      Ketoacyl-synt_C (ketosynthase)  General_PKS
PF00109.30      ketoacyl-synt (ketosynthase)    General_PKS
PF16073.9       SAT (starter acyltransferase)   NR
TIGR04532.1     PT (product template)   NR
PF21089.1       PKS_DH_N (dehydratase)  R
PF14765.10      PS-DH (dehhdratase)     R
PF08659.14      KR (ketoreductase)      R
PF00378.24      ECH_1 (enoyl reductase) R
PF16113.9       ECH_2 (enoyl reductase) R
PF08240.16      ADH_N (enoyl reductase) R
PF00107.30      ADH_zinc_N (enoyl reductase)    R
PF13602.10      ADH_zinc_N_2 (enoyl reductase)  R
PF02771.20      Acyl-CoA_dh_N (enoyl reductase) R
PF00441.28      Acyl-CoA_dh_1 (enoyl reductase) R
PF08028.15      Acyl-CoA_dh_2 (enoyl reductase) R
PF13561.10      adh_short_C2 (enoyl reductase)  R
PF00975.24      Thiosterase (release domain)    NR
PF14759.10      Reductase_C (release domain)    NR
PF00501.32      AMP-Binding enzyme (adenylation)        NRPS_Hybrid
"""

ENOYL_REDUCTASES = set(['PF00378.24', 'PF16113.9', 'PF08240.16', 'PF00107.30', 'PF13602.10', 'PF02771.20', 'PF00441.28', 'PF08028.15', 'PF13561.10'])
GENERALLY_REDUCING_DOMS = set(['PF13561.10', 'PF08028.15', 'PF00441.28', 'PF02771.20', 'PF13602.10', 'PF00107.30', 'PF08240.16', 'PF16113.9', 'PF00378.24', 'PF08659.14', 'PF14765.10', 'PF21089.1'])
GENERALLY_NONREDUCING_DOMS = set(['PF16073.9', 'TIGR04532.1', 'PF00975.24', 'PF14759.10'])

full_gbk_file = sys.argv[1]
hmm_file = sys.argv[2]
outdir = os.path.abspath(sys.argv[3]) + '/'

cds_faa = outdir + 'proteins.faa'
cds_faa_handle = open(cds_faa, 'w')
with open(full_gbk_file) as orgf:
    for rec in SeqIO.parse(orgf, 'genbank'):
        for feat in rec.features:
            if feat.type == "CDS":
                gene_id = feat.qualifiers["gene"][0]
                sequence = feat.qualifiers["translation"][0]
                cds_faa_handle.write('>' + gene_id + '\n' + str(sequence) + '\n')
cds_faa_handle.close()

hmmer_results = outdir + 'hmmsearch_results.txt'
hmmsearch_cmd = ['hmmsearch', '--cut_ga', '--tblout', hmmer_results, hmm_file, cds_faa]
os.system(' '.join(hmmsearch_cmd))

domains_hit = set([])
with open(hmmer_results) as ohr:
    for line in ohr:
        line = line.strip()
        if line.startswith('#'): continue
        ls = line.split()
        pf = ls[3]
        domains_hit.add(pf)

final_result = outdir + 'domains_found.txt'
final_result_handle = open(final_result, 'w')
final_result_handle.write('\n'.join(sorted(domains_hit)))
final_result_handle.close()
