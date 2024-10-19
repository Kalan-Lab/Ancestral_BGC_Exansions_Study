import os
import sys
from Bio import SeqIO

def parseAntiSMASHGBKForFunctionAndCompleteness(bgc_gbk):
    product = 'NA'
    products = set([])
    with open(bgc_gbk) as obg:
        for rec in SeqIO.parse(obg, 'genbank'):
            for feature in rec.features:
                if feature.type == 'protocluster':
                    try:
                        products.add(feature.qualifiers.get('product')[0])
                    except:
                        pass

    products_filt = products.intersection(set(['thioamide-NRP', 'NRPS', 'NRPS-like', 'NRP-metallophore']))
    product_string = ';'.join(sorted(products_filt))
    nrps_related = False
    if 'thioamide-NRP' in products or 'NRPS' in products or 'NRPS-like' in products or 'NRP-metallophore' in products:
        nrps_related = True
    return([nrps_related, product_string])

overview_file = '../../Overview_File.txt'
as_dir = '../All_antiSMASH_Results/'

gca_to_genus = {}
with open(overview_file) as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        gca_to_genus[ls[0]] = ls[1].split('_')[0]

count = 1
for s in os.listdir(as_dir):
    samp_dir = as_dir + s + '/'
    for f in os.listdir(samp_dir):
        if '.region' in f and f.endswith('.gbk') and not f == 'JADGJL010000003.1.region001.gbk':
            gbk_file = samp_dir + f
            nrps_related, product_string = parseAntiSMASHGBKForFunctionAndCompleteness(gbk_file)

            if nrps_related:
                conden_dom_coords = set([])
                with open(gbk_file) as ogf:
                    for rec in SeqIO.parse(ogf, 'genbank'):
                        for feature in rec.features:
                            if feature.type == 'aSDomain':
                                domain_type = feature.qualifiers["aSDomain"][0]
                                if 'Condensation' in domain_type: 
                                    start = None
                                    end = None
                                    direction = None
                                    all_coords = []
                                    if not 'join' in str(feature.location) and not 'order' in str(feature.location):
                                        start = min([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
                                        end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
                                        direction = str(feature.location).split('(')[1].split(')')[0]
                                        all_coords.append([start, end, direction])
                                    elif 'order' in str(feature.location):
                                        all_starts = []
                                        all_ends = []
                                        all_directions = []
                                        for exon_coord in str(feature.location)[6:-1].split(', '):
                                            start = min([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
                                            end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
                                            direction = exon_coord.split('(')[1].split(')')[0]
                                            all_starts.append(start)
                                            all_ends.append(end)
                                            all_directions.append(direction)
                                            all_coords.append([start, end, direction])
                                            start = min(all_starts)
                                            end = max(all_ends)
                                    else:
                                        all_starts = []
                                        all_ends = []
                                        all_directions = []
                                        for exon_coord in str(feature.location)[5:-1].split(', '):
                                            start = min([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
                                            end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
                                            direction = exon_coord.split('(')[1].split(')')[0]
                                            all_starts.append(start)
                                            all_ends.append(end)
                                            all_directions.append(direction)
                                            all_coords.append([start, end, direction])
                                        start = min(all_starts)
                                        end = max(all_ends)
                                    for pos in range(start-1, end):
                                        conden_dom_coords.add(pos)

                with open(gbk_file) as ogf:
                    for rec in SeqIO.parse(ogf, 'genbank'):
                        for feature in rec.features:
                            if feature.type == 'CDS':
                                sequence = feature.qualifiers["translation"][0]

                                start = None
                                end = None
                                direction = None
                                all_coords = []
                                if not 'join' in str(feature.location) and not 'order' in str(feature.location):
                                    start = min([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
                                    end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
                                    direction = str(feature.location).split('(')[1].split(')')[0]
                                    all_coords.append([start, end, direction])
                                elif 'order' in str(feature.location):
                                    all_starts = []
                                    all_ends = []
                                    all_directions = []
                                    for exon_coord in str(feature.location)[6:-1].split(', '):
                                        start = min([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
                                        end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
                                        direction = exon_coord.split('(')[1].split(')')[0]
                                        all_starts.append(start)
                                        all_ends.append(end)
                                        all_directions.append(direction)
                                        all_coords.append([start, end, direction])
                                        start = min(all_starts)
                                        end = max(all_ends)
                                else:
                                    all_starts = []
                                    all_ends = []
                                    all_directions = []
                                    for exon_coord in str(feature.location)[5:-1].split(', '):
                                        start = min([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
                                        end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
                                        direction = exon_coord.split('(')[1].split(')')[0]
                                        all_starts.append(start)
                                        all_ends.append(end)
                                        all_directions.append(direction)
                                        all_coords.append([start, end, direction])
                                    start = min(all_starts)
                                    end = max(all_ends)
                                conden_related = False
                                for pos in range(start-1, end):
                                    if pos in conden_dom_coords:
                                        conden_related = True
                                if conden_related:
                                    label = 'NA'
                                    try:
                                        label = feature.qualifiers["locus_tag"][0]
                                    except:
                                        label = feature.qualifiers["gene"][0]

                                    print('>' + str(count) + ' ' + gca_to_genus[s] + '|' + s + '|' + label + '\n' + sequence)
                                    count += 1
