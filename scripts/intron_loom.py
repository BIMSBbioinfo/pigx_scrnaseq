import sys
import loompy
import pandas
import re
import numpy 

def load_and_order(path):
	ds = loompy.connect(path)
	genes = ds.ra.Genes
	ordering = genes.argsort()
	ds.permute(ordering, axis = 0)
	return(ds)

def create_intron_loom(path_exon, path_gene):
	ds_exon = load_and_order(path_exon)
	ds_gene = load_and_order(path_gene)

	intron_umi =  ds_gene[:,0:2] - ds_exon[:,0:2]

	for i in range(2, ds_gene.shape[1], 100):
		j = i + 100
		if j >= ds_gene.shape[1]:
			j = ds_gene.shape[1]
		tmp = ds_gene[:,i:j] - ds_exon[:,i:j]
		intron_umi = numpy.concatenate((intron_umi, tmp), axis = 1)

	row_attrs = { "Genes": list(ds_exon.ra.Genes)}
	col_attrs = { "cell_id": list(ds_exon.ca.cell_id) }

	ds_exon.close()
	ds_gene.close()

	return(intron_umi, row_attrs, col_attrs)

if __name__ == '__main__':

	input_exon = sys.argv[1]
	input_gene = sys.argv[2]
	output_filepath = sys.argv[3]

	intron_umi, row_attrs, col_attrs = create_intron_loom(input_exon, input_gene)

	loompy.create(output_filepath, intron_umi, row_attrs, col_attrs)