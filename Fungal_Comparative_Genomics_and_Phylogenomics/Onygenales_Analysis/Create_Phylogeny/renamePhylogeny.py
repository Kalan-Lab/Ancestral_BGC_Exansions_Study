import os
import sys
from ete3 import Tree

gca_to_species_name = {}
with open('Ony_Overview_File.txt') as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        species_name = '_'.join(ls[1].split(';')[-1].split())
        gca_to_species_name[gca] = species_name  + '_' + gca

gca_to_species_name['GCA_002240705'] = 'Elaphomyces_granulatus_GCA_002240705'

t = Tree('IQTRee_GHOST_SCC.treefile')
for n in t.traverse('postorder'):
    if n.is_leaf():
        gca = n.name
        sp_name = gca_to_species_name[gca]
        n.name = sp_name

t.write(outfile='Onygenales_with_Outgroups.tre', format=1)
