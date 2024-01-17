import os
import sys

for line in open('data.txt'):
    if 'PRJNA' in line:
        print('PRJNA' + line.split('PRJNA')[1].split()[0])
