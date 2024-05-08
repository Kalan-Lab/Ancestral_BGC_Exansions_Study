import os
import sys
import pyhmmer

## Coverage solution taken from apcamargo's answer in a pyhmmer ticket on Github: https://github.com/althonos/pyhmmer/issues/27

db_file = 'dbCAN-HMMdb-V12.txt'
z = 783
cpus = 50

hmm_lengths = {}
try:
    with pyhmmer.plan7.HMMFile(db_file) as hmm_file:
        for hmm in hmm_file:
            hmm_lengths[hmm.name] = len(hmm.consensus)
except:
    raise RuntimeError("Problem getting HMM consensus lengths!")

protein_faa = sys.argv[1]
annotation_result_file = sys.argv[2]

try:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = []
    with pyhmmer.easel.SequenceFile(protein_faa, digital=True, alphabet=alphabet) as seq_file:
        sequences = list(seq_file)

        outf = open(annotation_result_file, 'w')
        with pyhmmer.plan7.HMMFile(db_file) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, sequences, Z=int(z), cpus=cpus):
                for hit in hits:
                    n_aligned_positions = len(hit.best_domain.alignment.hmm_sequence) - hit.best_domain.alignment.hmm_sequence.count(".")
                    hmm_coverage = (n_aligned_positions / hmm_lengths[hit.best_domain.alignment.hmm_name])
                    outf.write('\t'.join([hits.query_name.decode(), 'NA', hit.name.decode(), 'NA', str(hit.evalue), str(hit.score), str(hmm_coverage)]) + '\n')
        outf.close()
except:
    raise RuntimeError('Problem running pyhmmer!')

