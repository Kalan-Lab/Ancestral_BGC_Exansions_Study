import os
import sys
from Bio import SeqIO
from collections import defaultdict
import traceback

def parseCDSCoord(str_gbk_loc):
	try:
		start = None
		end = None
		direction = None
		all_coords = []
		is_multi_part = False
		if not 'join' in str(str_gbk_loc) and not 'order' in str(str_gbk_loc):
			start = min([int(x.strip('>').strip('<')) for x in
						 str(str_gbk_loc)[1:].split(']')[0].split(':')]) + 1
			end = max([int(x.strip('>').strip('<')) for x in
					   str(str_gbk_loc)[1:].split(']')[0].split(':')])
			direction = str(str_gbk_loc).split('(')[1].split(')')[0]
			all_coords.append([start, end, direction])
		elif 'order' in str(str_gbk_loc):
			is_multi_part = True
			all_starts = []
			all_ends = []
			all_directions = []
			for exon_coord in str(str_gbk_loc)[6:-1].split(', '):
				ec_start = min(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
				ec_end = max(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
				ec_direction = exon_coord.split('(')[1].split(')')[0]
				all_starts.append(ec_start)
				all_ends.append(ec_end)
				all_directions.append(ec_direction)
				all_coords.append([ec_start, ec_end, ec_direction])
			assert (len(set(all_directions)) == 1)
			start = min(all_starts)
			end = max(all_ends)
			direction = all_directions[0]
		else:
			is_multi_part = True
			all_starts = []
			all_ends = []
			all_directions = []
			for exon_coord in str(str_gbk_loc)[5:-1].split(', '):
				ec_start = min(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
				ec_end = max(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
				ec_direction = exon_coord.split('(')[1].split(')')[0]
				all_starts.append(ec_start)
				all_ends.append(ec_end)
				all_directions.append(ec_direction)
				all_coords.append([ec_start, ec_end, ec_direction])
			assert (len(set(all_directions)) == 1)
			start = min(all_starts)
			end = max(all_ends)
			direction = all_directions[0]
		return(all_coords, start, end, direction, is_multi_part)
	except Exception as e:
		raise RuntimeError(traceback.format_exc())

def parseAntiSMASHGBKForProductInfo(bgc_gbk):
        products_by_location = defaultdict(set)
        #try:
        with open(bgc_gbk) as obg:
            for rec in SeqIO.parse(obg, 'genbank'):
                for feat in rec.features:
                    if feat.type == 'protocluster':
                        all_coords, start, end, direction, is_multi_part = parseCDSCoord(str(feat.location))
                        try:
                            product = feat.qualifiers.get('product')[0]
                            for pos in range(start, end+1):
                                products_by_location[pos].add(product)
                        except:
                            pass
        #except:
        #    print('error processing genbank')
        #    sys.exit(1)
        return(products_by_location)

as_dir = os.path.abspath(sys.argv[1]) + '/'
for s in os.listdir(as_dir):
    samp_dir = as_dir + s + '/'
    for f in os.listdir(samp_dir):
        if '.region' in f and f.endswith('.gbk'):
            bgc_gbk = samp_dir + f 
            products_by_location = parseAntiSMASHGBKForProductInfo(bgc_gbk)
            
            with open(bgc_gbk) as ogf:
                for rec in SeqIO.parse(ogf, 'genbank'):
                    for feat in rec.features:
                        if feat.type != 'CDS': continue
                        try:
                            lt = feat.qualifiers.get('locus_tag')[0]
                        except:
                            lt = feat.qualifiers.get('gene')[0]
                        rule_based_bgc_cds = False
                        try:
                            if 'rule-based-clusters' in feat.qualifiers.get('gene_functions')[0]:
                                rule_based_bgc_cds = True
                        except:
                            pass
                        if not rule_based_bgc_cds: continue
                        translation = feat.qualifiers.get('translation')[0]
                        all_coords, start, end, direction, is_multi_part = parseCDSCoord(str(feat.location))
                        products = set([])
                        for pos in range(start, end+1):
                            products = products.union(products_by_location[pos])
                        product = '|'.join(products)
                        print('>' + s + '|' + lt + '|' + f.split('.gbk')[0] + '|' + product + '\n' + translation)
