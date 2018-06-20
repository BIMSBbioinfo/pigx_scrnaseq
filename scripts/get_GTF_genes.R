#!/usr/bin/env Rscript

# Author: EU
# Date: June, 2018
# This script takes as input a GTF file and returns annotation for genes

# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('get_GTF_genes')
# -------------------------------------------------------------------------- #
library(rtracklayer)
library(GenomicRanges)

get_GTF_genes = function(
    gtfFile = NULL,
    output  = NULL

) {
    if (is.null(gtfFile))
        stop("gtfFile not specified")
    if (is.null(output))
        stop("output file not specified")

    gtf   = import.gff2(gtfFile)
    g     = subset(gtf, type == "exon")
    annot = values(gtf)[,c('gene_id','gene_name')]
    annot = unique(annot)
    annot$exon_id         = annot$gene_id
    annot$transcript_id   = annot$gene_id
    annot$transcript_name = annot$gene_id
    annot$type = 'exon'
    annot$exon_number = 1

    gl  = split(g, g$gene_id)
    gl  = unlist(range(gl))
    gl$gene_id = names(gl)
    values(gl) = merge(values(gl), annot, by='gene_id')


    if(!all(names(gl) == gl$gene_id))
        stop('gene_id is not properly ordered')

    export.gff2(gl, output)
}

# -------------------------------------------------------------------------- #
get_GTF_genes(
    gtfFile = argv$input[['infile']],
    output  = argv$output[['outfile']]
)
