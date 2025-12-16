#!/usr/bin/env python3
"""
Parse GenBank file and annotation TSV to create a simple TSV for gggenes plotting
"""

import pandas as pd
from Bio import SeqIO
import sys

# Configuration
GBK_FILE = "GCA_003751085.1_ASM375108v1_fai-gene-cluster-4.gbk"
TSV_FILE = "Consolidated_Report.tsv"
START_GENE = "AAAC_000067"
END_GENE = "AAAC_000134"
OUTPUT_FILE = "genes_for_plotting.tsv"

def parse_tsv_data(tsv_file):
    """Parse the consolidated report TSV file."""
    df = pd.read_csv(tsv_file, sep='\t')
    
    # Create a dictionary mapping locus tags to their annotations
    annotation_dict = {}
    for _, row in df.iterrows():
        locus_tags = str(row['CDS Locus Tags'])
        if pd.notna(locus_tags) and '|' in locus_tags:
            # Extract the last element which contains the locus tag
            parts = locus_tags.split('|')
            locus_tag = parts[-1] if len(parts) > 0 else None
            
            if locus_tag:
                annotation_dict[locus_tag] = {
                    'manual_annotation': row.get('manual annotation', 'Other'),
                    'conservation_status': row.get('Conservation Status', ''),
                    'og_id': row.get('Ortholog Group (OG) ID', '')
                }
    
    return annotation_dict

def extract_genes_from_gbk(gbk_file, start_gene, end_gene):
    """Extract gene features from GenBank file between start and end genes."""
    genes = []
    record = SeqIO.read(gbk_file, "genbank")
    
    in_range = False
    for feature in record.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
            
            # Check if we've reached the start gene
            if locus_tag == start_gene:
                in_range = True
            
            # If we're in range, add the gene
            if in_range:
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                
                genes.append({
                    'locus_tag': locus_tag,
                    'start': start,
                    'end': end,
                    'strand': '+' if strand == 1 else '-',
                })
            
            # Check if we've reached the end gene
            if locus_tag == end_gene:
                break
    
    return genes

def main():
    """Main function to create the plotting data."""
    print("Parsing TSV file...")
    annotation_dict = parse_tsv_data(TSV_FILE)
    print(f"Found {len(annotation_dict)} annotated genes")
    
    print(f"Extracting genes from {START_GENE} to {END_GENE}...")
    genes = extract_genes_from_gbk(GBK_FILE, START_GENE, END_GENE)
    print(f"Found {len(genes)} genes in range")
    
    # Create DataFrame for plotting
    plot_data = []
    for gene in genes:
        locus_tag = gene['locus_tag']
        info = annotation_dict.get(locus_tag, {
            'manual_annotation': 'Other',
            'conservation_status': '',
            'og_id': ''
        })
        
        manual_annot = info['manual_annotation']
        if pd.isna(manual_annot) or manual_annot == '':
            manual_annot = 'Other'
        
        conservation = info['conservation_status']
        conserved = 'Yes' if (conservation and str(conservation).strip() and 
                             str(conservation).lower() not in ['', 'nan', 'na']) else 'No'
        
        plot_data.append({
            'locus_tag': locus_tag,
            'start': gene['start'],
            'end': gene['end'],
            'strand': gene['strand'],
            'manual_annotation': manual_annot,
            'conserved': conserved,
            'molecule': 'Cluster'  # For gggenes
        })
    
    # Create DataFrame and save
    df = pd.DataFrame(plot_data)
    df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    print(f"\n✓ Data saved to: {OUTPUT_FILE}")
    
    # Print summary
    print("\nGene Summary:")
    print(f"{'Locus Tag':<15} {'Function':<35} {'Strand':<8} {'Conserved':<10}")
    print("-" * 70)
    for _, row in df.iterrows():
        conserved_mark = '✓' if row['conserved'] == 'Yes' else ''
        print(f"{row['locus_tag']:<15} {row['manual_annotation']:<35} {row['strand']:<8} {conserved_mark:<10}")

if __name__ == "__main__":
    main()

