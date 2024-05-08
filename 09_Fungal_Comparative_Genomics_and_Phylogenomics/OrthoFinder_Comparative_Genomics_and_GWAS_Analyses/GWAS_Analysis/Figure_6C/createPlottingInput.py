import os
import sys

"""
0       variant
1       af
2       filter-pvalue
3       lrt-pvalue
4       beta
5       key BGC proportion
6       BGC found proportion
7       eggnog annotations
8       Annotation
9       Annotation Category
10      eggnog cogs
11      Notes
"""

print('HOG\tAF\tLRT_Pvalue\tAnnotation_Category')
with open('GWAS_Data.txt') as ogdf:
    for i, line in enumerate(ogdf):
        if i == 0: continue
        line = line.strip('\n')
        ls = line.split('\t')
        hog = ls[0]
        af = ls[1]
        lp = ls[3]
        annot_cat = ls[9]
        print(hog + '\t' + af + '\t' + lp + '\t' + annot_cat)
