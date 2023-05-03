# Jared Bard
# 4/25/2023
# Almost entirely written by ChatGPT4

import sys

# Define a function to parse the GFF3 file and extract gene information
def parse_gff3(gff3_file):
    genes = {}
    with open(gff3_file, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            # Only process gene features
            if fields[2] != 'gene':
                continue
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            gene_id = fields[8].split(';')[0].split('=')[1]
            genes[gene_id] = (chrom, start, end, strand)
    return genes

# Define a function to check if one gene is entirely inside another gene
def is_inside(gene1, gene2):
    if gene1[0] != gene2[0] or gene1[3] != gene2[3]:
        # Genes are on different chromosomes or strands, so gene1 is not inside gene2
        return False
    if gene1[1] >= gene2[1] and gene1[2] <= gene2[2]:
        # Gene1 is entirely inside gene2
        return True
    else:
        return False

# Define a function to calculate the distance between two genes on the same strand
def calc_distance(gene1, gene2):
    if gene1[0] != gene2[0] or gene1[3] != gene2[3]:
        # Genes are on different chromosomes or strands, so distance is undefined
        return None
    if gene1[2] < gene2[1]:
        # Gene 1 is upstream of gene 2
        return gene2[1] - gene1[2]
    elif gene2[2] < gene1[1]:
        # Gene 2 is upstream of gene 1
        return gene1[1] - gene2[2]
    else:
        # Genes overlap or one is inside the other, so distance is undefined
        return None

# Parse the GFF3 file and extract gene information
genes = parse_gff3(sys.argv[1])
#gff3_file = "/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/UTRs/Saccharomyces_cerevisiae.R64-1-1.109.gff3"
# Loop over the genes and find the closest gene on the same strand
distances = {}
for gene_id in genes:
    gene = genes[gene_id]
    closest_gene_id = None
    closest_distance = None
    for other_gene_id in genes:
        if other_gene_id == gene_id:
            continue
        other_gene = genes[other_gene_id]
        if gene[0] != other_gene[0] or gene[3] != other_gene[3]:
            # Skip genes that are not on the same chromosome and strand
            continue
        if is_inside(other_gene, gene): # modified to check if other gene is inside gene
            # Skip genes that are entirely inside another gene
            continue
        if gene[3] == '+' and other_gene[2] < gene[1]:
            # First gene is on the forward strand and other gene ends before first gene starts
            distance = gene[1] - other_gene[2]
            if closest_distance is None or distance < closest_distance:
                closest_gene_id = other_gene_id
                closest_distance = distance
        elif gene[3] == '-' and other_gene[1] > gene[2]:
            # First gene is on the reverse strand and other gene starts after first gene ends
            distance = other_gene[1] - gene[2]
            if closest_distance is None or distance < closest_distance:
                closest_gene_id = other_gene_id
                closest_distance = distance
    if closest_gene_id is not None:
        distances[gene_id] = (closest_gene_id, closest_distance)

# Print the distances between neighboring genes
print('Gene\tClosest_gene\tDistance')
for gene_id in distances:
    closest_gene_id, distance = distances[gene_id]
    print('{}\t{}\t{}'.format(gene_id.replace("gene:",""),
                              closest_gene_id.replace("gene:",""), distance))

