import sys
import loompy
import pandas as pd
import re
import numpy
import os
from scipy.io import mmread
import argparse

# ------------------------------------------------------------------ #
#find the complete list of gene ids from the gtf file
def get_gtf_gene_ids (gtf_file):
    gene_ids = []
    with open(gtf_file) as f:
        i = 0
        for line in f:
            m = re.search('gene_id "(.+?)"', line)
            if m != None:
                gene_ids.append(m.group(1))
    return(set(gene_ids))

# ------------------------------------------------------------------ #
#given a matrix, row_names, and col_names, find which row names don't exist in gene_ids
#and fill the matrix with rows of zeros for each missing gene
def fill_missing_genes_into_matrix(gene_ids, matrix, row_names, col_names):
    missing = list(set(gene_ids) - set(row_names))

    # append the matrix with rows of zeros
    if(len(missing) > 0):
        new_rows  = numpy.zeros((len(missing), len(col_names)), dtype = int)
        matrix    = numpy.vstack((matrix, new_rows))
        row_names = row_names + missing

    return(matrix, row_names)

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert STAR mtx to loom')
    parser.add_argument('--sample_id',   action="store", dest='sample_id')
    parser.add_argument('--input_file',  action="store", dest="input_file")
    parser.add_argument('--gtf_file',    action="store", dest="gtf_file")
    parser.add_argument('--output_file', action="store", dest="output_file")

    args = parser.parse_args()
    # -------------------------------------------------------------- #
    sample_id    = args.sample_id
    input_file   = args.input_file
    gtf_file     = args.gtf_file
    output_file  = args.output_file

    # -------------------------------------------------------------- #
    print("Reading input files ...")
    matrix = mmread(input_file)

    basepath   = os.path.dirname(input_file)
    path_genes = os.path.join(basepath,'genes.tsv')
    genes      = pd.read_csv(path_genes, sep='\t', header=None)
    genes.columns = ['gene_id','gene_id2']
    genes      = genes.drop(columns = ['gene_id2'])

    path_barcode    = os.path.join(basepath,'barcodes.tsv')
    barcode         = pd.read_csv(path_barcode, sep='\t', header=None)
    barcode.columns = ['cell_id']
    barcode['sample_name'] = sample_id

    # -------------------------------------------------------------- #
    print("Parsing gene ids from gtf file",gtf_file)
    gene_ids = get_gtf_gene_ids(gtf_file)

    # -------------------------------------------------------------- #
    # fills the matrix with zeros for missing genes
    # all input matrices have the same X dimension
    umi_matrix, row_names = fill_missing_genes_into_matrix(gene_ids, matrix.toarray(), genes['gene_id'], barcode['cell_id'])

    # -------------------------------------------------------------- #
    print('Creating loompy file')
    loompy.create(
        output_file,
        matrix.toarray(),
        genes.to_dict("list"),
        barcode.to_dict("list")
    )