import json
from Bio import Entrez
Entrez.email = 'salamzade@wisc.edu'

def accession2taxid(acc: str, db="nucleotide") -> str:
    handle = Entrez.esearch(db=db, term=acc)
    record = Entrez.read(handle)
    gi = record["IdList"][0]
    handle = Entrez.esummary(db=db, id=gi, retmode="json")
    result = json.load(handle)["result"]
    taxid = result[gi]["taxid"]
    return str(taxid)

"""
with open('All_UFCG_Representative_Selections.txt') as ogf:
    for line in ogf:
        tax, gca = line.strip().split('\t')
        ncbi_taxid = accession2taxid(gca, db="assembly")
        print(tax + '\t' + gca + '\t' + ncbi_taxid)
"""

with open('Previous_Overview_File.txt') as opof:
    for i, line in enumerate(opof):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca = ls[3]

        ncbi_taxid = accession2taxid(gca, db="assembly")
        print(gca + '\t' + ncbi_taxid)
