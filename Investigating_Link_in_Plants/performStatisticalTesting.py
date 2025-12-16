import os
import sys
from scipy import stats

vascular = []
early = []
with open('Data_from_PlantiSMASH_V2_Preprint_SupTableS4.with_grouping_added.txt') as odf:
    for i, line in enumerate(odf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        if 'Vascular' in ls[-1]:
            vascular.append(int(ls[-2]))
        else:
            early.append(int(ls[-2]))

print(stats.ranksums(vascular, early, alternative='greater'))
