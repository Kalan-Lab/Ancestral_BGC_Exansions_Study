import os
import sys

genera = ['Actinokineospora', 'Amycolatopsis', 'Kutzneria', 'Micromonospora', 'Streptomyces']
genus_colors = ['#FF66CC', '#FF66CC', '#FF66CC', '#FF66CC', '#9268F6']

print('DATASET_PIECHART')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tGenus')
print('FIELD_LABELS\t' + '\t'.join(genera))
print('FIELD_COLORS\t' + '\t'.join(genus_colors))
print('DATA')
with open('Genome_to_Genus.txt') as oggf:
    for line in oggf:
        line = line.strip()
        ls = line.split('\t')
        row = [ls[0], '1', '1']
        for g in genera:
            if g == ls[1].strip():
                row.append('1')
            else:
                row.append('0')
        print('\t'.join(row))
