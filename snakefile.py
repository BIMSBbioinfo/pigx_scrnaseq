
"""
Snakefile for pigx-scrnaseq pipeline
"""

# ----------------------------------------------------------------------------- #
# libraries and constants
import glob
import os
import re
import subprocess
import yaml
import csv
import inspect
import magic as mg
import sys
import pandas as pd

PATH_SCRIPT = os.path.join(config['locations']['pkglibexecdir'], 'scripts')

# loads function for running Rscripts
include: os.path.join(PATH_SCRIPT, 'Run_Rscript.py')
# functions for java input
include: os.path.join(PATH_SCRIPT, 'Accessory_Functions.py')
include: os.path.join(PATH_SCRIPT, 'validate_input.py')
validate_config(config)
# class for SAMPLE_SHEET
include: os.path.join(PATH_SCRIPT, 'Sample_Sheet_Class.py')

# ----------------------------------------------------------------------------- #
# Software parameters
SOFTWARE = config['tools']
# variables

GENOME_NAME_PRIMARY = config['annotation']['primary']['genome']['name']
REFERENCE_NAMES     = [GENOME_NAME_PRIMARY]
COVARIATES          = config['covariates']

# ----------------------------------------------------------------------------- #
# adapter locations
ADAPTER_PARAMETERS = config['adapter_parameters']

# ----------------------------------------------------------------------------- #
# output count matrices for STAR solo
# Gene and GeneFull is obligatory, the rest might be set as params in conig
# SJ matrix is absent because STAR output is currently wrong
# implemented in two lists because it's the easiest way to pass to functions
STAR_OUTPUT_TYPES_KEYS = ['Gene', 'GeneFull', 'Velocyto', 'Velocyto']
STAR_OUTPUT_TYPES_VALS = ['Counts', 'GeneFull', 'Spliced', 'Unspliced']
# This is a relict variable
STAR_OUTPUT_TYPES = list(set(STAR_OUTPUT_TYPES_KEYS))

# ----------------------------------------------------------------------------- #
# PATHS
OUTPUT_DIR           = config['locations']['output-dir']
PATH_FASTQ           = config['locations']['reads-dir']
TEMPDIR              = config['locations']['tempdir']
PATH_ANNOTATION      = os.path.join(OUTPUT_DIR, 'Annotation')
PATH_FASTQC          = os.path.join(OUTPUT_DIR, 'FASTQC')
PATH_MAPPED          = os.path.join(OUTPUT_DIR, 'Mapped')
PATH_LOG             = os.path.join(OUTPUT_DIR, 'Log')
PATH_SAMPLE_SHEET    = config['locations']['sample-sheet']
PATH_RSCRIPT         = SOFTWARE['Rscript']['executable']

# used in automatic recognition of parameters from the settings files
PARAMS               = config['general']['params']
SOFTWARE             = config['tools']

PATH_ANNOTATION_PRIMARY = os.path.join(PATH_ANNOTATION, GENOME_NAME_PRIMARY)
PATH_REFERENCE_PRIMARY  = config['annotation']['primary']['genome']['fasta']
PATH_GTF_PRIMARY        = config['annotation']['primary']['gtf']

GENOME_NAME_MIX = None
if config['annotation']['secondary']:
    GENOME_NAME_SECONDARY    = config['annotation']['secondary']['genome']['name']
    PATH_REFERENCE_SECONDARY = config['annotation']['secondary']['genome']['fasta']
    PATH_GTF_SECONDARY       = config['annotation']['secondary']['gtf']

    GENOME_NAME_MIX = GENOME_NAME_PRIMARY + '_' + GENOME_NAME_SECONDARY
    PATH_ANNOTATION_MIX = os.path.join(PATH_ANNOTATION, GENOME_NAME_MIX)
    REFERENCE_NAMES = REFERENCE_NAMES + [GENOME_NAME_MIX]

# extracts the mitochondrial chromosome, if defined
MITOCHNODRIAL_CHROMOSOME = 'chrM'
if 'mitochondrial_chromosome' in set(config['annotation']['primary']['genome'].keys()):
    MITOCHNODRIAL_CHROMOSOME = config['annotation']['primary']['genome']['mitochondrial_chromosome']

## Load sample sheet
SAMPLE_SHEET = experiment(config = config)
SAMPLE_SHEET.init_SAMPLE_SHEET(PATH_SAMPLE_SHEET)

SAMPLE_NAMES = SAMPLE_SHEET.fetch_sample_names()

# ----------------------------------------------------------------------------- #
# sets the temporrary directory to default in the working directory if the tempdir does not exist
if TEMPDIR == None:
    if "TMPDIR" in os.environ.keys():
        TEMPDIR = os.environ['TMPDIR']
    else:
        TEMPDIR = '/tmp'

# ----------------------------------------------------------------------------- #
# RULES

# ----------------------------------------------------------------------------- #
# Link primary reference
# TODO : make extension dynamic (to accept both fasta and fasta.gz files)
LINK_REFERENCE_PRIMARY   = os.path.join(PATH_ANNOTATION_PRIMARY,  GENOME_NAME_PRIMARY + '.fasta')
LINK_GTF_PRIMARY         = os.path.join(PATH_ANNOTATION_PRIMARY, GENOME_NAME_PRIMARY + '.gtf')


# ----------------------------------------------------------------------------- #
# Combine primary and secondary reference genomes
COMBINE_REFERENCE = []
GENOME_SECONDARY_IND = not GENOME_NAME_MIX == None
if GENOME_SECONDARY_IND:
    PATH_REFERENCE_MIX = os.path.join(PATH_ANNOTATION_MIX, GENOME_NAME_MIX + '.fasta')
    PATH_GTF_MIX = os.path.join(PATH_ANNOTATION_MIX, GENOME_NAME_MIX + '.gtf')
    COMBINE_REFERENCE = COMBINE_REFERENCE + [PATH_REFERENCE_MIX, PATH_GTF_MIX]

# ----------------------------------------------------------------------------- # FastQC
FASTQC_prep  = list(SAMPLE_SHEET.SAMPLE_SHEET['reads']) + list(SAMPLE_SHEET.SAMPLE_SHEET['barcode'])
FASTQC = [os.path.join(PATH_FASTQC, file + ".fastqc.done") for file in FASTQC_prep]

# ----------------------------------------------------------------------------- #
# Change reference gene_name to gene_id
PATH_GTF_PRIMARY_ID = expand(os.path.join(PATH_ANNOTATION_PRIMARY, "{name}" + "gene_id.gtf"), name=REFERENCE_NAMES)

# ----------------------------------------------------------------------------- #
# STAR INDEX
MAKE_STAR_INDEX = expand(os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX','done.txt'), genome = REFERENCE_NAMES)


# -----------------------------------------------------------------------------
#
# merge reads
# use for loop - expand gives wrong results
MERGE_TECHNICAL_REPLICATES = []
for sample_name in SAMPLE_NAMES:
    MERGE_TECHNICAL_REPLICATES.append(os.path.join(PATH_MAPPED, sample_name, SAMPLE_SHEET.fetch_field(sample_name,'reads_merged')))
    MERGE_TECHNICAL_REPLICATES.append(os.path.join(PATH_MAPPED, sample_name, SAMPLE_SHEET.fetch_field(sample_name,'barcode_merged')))

# -----------------------------------------------------------------------------
#
# filters reads using flexbar
FILTER_READS = []
for sample_name in SAMPLE_NAMES:
    for genome in REFERENCE_NAMES:

        barcode_file = os.path.join(PATH_MAPPED, sample_name, genome, sample_name + '_1.fastq.gz')
        FILTER_READS.append(barcode_file)

        read_file    = os.path.join(PATH_MAPPED, sample_name, genome, sample_name + '_2.fastq.gz')
        FILTER_READS.append(read_file)

# ----------------------------------------------------------------------------- #
# MAPPING
MAP_scRNA = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}", "{name}_Aligned.out.bam"), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# SORT and INDEX BAM
SORT_BAM = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}","{name}.sorted.bam"), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

INDEX_BAM = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}","{name}.sorted.bam.bai"), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# Number of reads per cell calculation
BAM_HISTOGRAM = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_cell_barcode_histogram.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# Number of reads per cell calculation
FIND_CELL_NUMBER_CUTOFF = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadCutoff.yaml'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# UMI matrix in loom format
UMI_LOOM =  expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.matrix.loom'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# Combined UMI matrices in loom format
COMBINED_LOOM_MATRICES = expand(os.path.join(PATH_MAPPED, "{genome}_UMI.loom"), genome = REFERENCE_NAMES)


# ----------------------------------------------------------------------------- #
# Import and preprocess the combined loom files and save as SingleCellExperiment.RDS objects.
SCE_RDS_FILES = expand(os.path.join(PATH_MAPPED, "{genome}.SingleCellExperiment.RDS"), genome = REFERENCE_NAMES)


# ----------------------------------------------------------------------------- #
# Seurat RDS files
SEURAT_RDS_FILES = expand(os.path.join(PATH_MAPPED, "{genome}.Seurat.RDS"), genome = REFERENCE_NAMES)


# ----------------------------------------------------------------------------- #
# READ statistics
READ_STATISTICS = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadStatistics.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)


# ----------------------------------------------------------------------------- #
## Using the preprocessed SingleCellExperiment.RDS file, generates a self-contained HTML report
REPORT_FILES = expand(os.path.join(PATH_MAPPED, "{genome}.scRNA-Seq.report.html"), genome = REFERENCE_NAMES)
# ----------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------- #
# Bam To BigWig
BIGWIG = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}.bw'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
RULE_ALL = []
RULE_ALL = RULE_ALL + [LINK_REFERENCE_PRIMARY, LINK_GTF_PRIMARY]

if len(COMBINE_REFERENCE) > 0:
    RULE_ALL = RULE_ALL + COMBINE_REFERENCE


RULE_ALL = RULE_ALL + MAKE_STAR_INDEX + MERGE_TECHNICAL_REPLICATES + FILTER_READS + BAM_HISTOGRAM + FIND_CELL_NUMBER_CUTOFF + MAP_scRNA + SORT_BAM + INDEX_BAM + UMI_LOOM + COMBINED_LOOM_MATRICES + SCE_RDS_FILES + SEURAT_RDS_FILES + BIGWIG + READ_STATISTICS + REPORT_FILES

# ----------------------------------------------------------------------------- #
rule all:
    input:
        RULE_ALL



# -----------------------------------------------------------------------------
# ----------------------------------------------------------------------------- #
# links the primary annotation to the ./Annotation folder
rule link_primary_annotation:
    input:
        gtf   = PATH_GTF_PRIMARY,
        fasta = PATH_REFERENCE_PRIMARY
    output:
        gtf   = LINK_GTF_PRIMARY,
        fasta = LINK_REFERENCE_PRIMARY
    params:
        threads = config['execution']['rules']['link_primary_annotation']['threads'],
        mem     = config['execution']['rules']['link_primary_annotation']['memory'],
        #ln      = SOFTWARE['ln']['executable']
        ln = 'ln'
    message:
        """
            Linking primary reference files:
                gtf:
                    file: {input.gtf}
                    link: {output.gtf}
                fasta:
                    file: {input.fasta}
                    link: {output.fasta}
        """
    shell:"""
        {params.ln} -s {input.gtf} {output.gtf}
        {params.ln} -s {input.fasta} {output.fasta}
    """



# ----------------------------------------------------------------------------- #
if GENOME_SECONDARY_IND:
    rule combine_reference:
        input:
            primary   =  LINK_REFERENCE_PRIMARY,
            secondary =  PATH_REFERENCE_SECONDARY
        output:
            outfile = PATH_REFERENCE_MIX
        params:
            threads = config['execution']['rules']['combine_reference']['threads'],
            mem     = config['execution']['rules']['combine_reference']['memory'],
            genome_name_primary   = GENOME_NAME_PRIMARY,
            genome_name_secondary = GENOME_NAME_SECONDARY,
            perl = SOFTWARE['perl']['executable'],
            perl_args = SOFTWARE['perl']['executable']['args']
        message:
            """
                Combining fasta files:
                    primary   : {input.primary}
                    secondary : {input.secondary}
                    output: {output.outfile}
            """
        shell:"""
            {params.perl} {params.perl_args} -pe 's|^>|>{params.genome_name_primary}|' < {input.primary} > {output.outfile}
            {params.perl} {params.perl_args} -pe 's|^>|>{params.genome_name_secondary}|' < {input.secondary} >> {output.outfile}
    """

# ----------------------------------------------------------------------------- #
# changes the gene_name field in the GTF file to the gene_id
# this is required for droptools counting
rule change_gtf_id:
    input:
        infile = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.gtf')
    output:
        outfile = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.gene_id.gtf')
    params:
        threads = config['execution']['rules']['change_gtf_id']['threads'],
        mem     = config['execution']['rules']['change_gtf_id']['memory'],
        script  = PATH_SCRIPT,
        Rscript = PATH_RSCRIPT,
        perl    = SOFTWARE['perl']['executable']
    message:
        """
            Changing GTF id:
                input  : {input}
                output : {output}
        """
    run:
        RunRscript(input, output, params, params.script, 'change_gtf_id.R')

        # removes header from the gtf file
        command = " ".join([
            params.perl, "-pi -e 's|^#.+$||'", output.outfile
        ])
        shell(command)


# ----------------------------------------------------------------------------- #
# STAR INDEX
###
rule make_star_reference:
    input:
        fasta = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.fasta'),
        gtf   = rules.change_gtf_id.output.outfile
        # gtf   = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.gtf')
    output:
        outfile = os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX','done.txt')
    params:
        outdir      = os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX'),
        star        = SOFTWARE['STAR']['executable'],
        star_args   = SOFTWARE['STAR']['args'],
        touch       = SOFTWARE['touch']['executable'],
        threads     = config['execution']['rules']['make_star_reference']['threads'],
        mem         = config['execution']['rules']['make_star_reference']['memory'],
        params_STAR = PARAMS['make_star_reference']
    log:
        logfile = os.path.join(PATH_LOG, '{genome}.make_star_reference.log')
    message:"""
        Star reference:
            input:
                fasta : {input.fasta}
                gtf   : {input.gtf}
        """
    run:
        # star genome build command
        command = " ".join(
        [params.star, params.star_args,
         '--runMode genomeGenerate',
        '--genomeDir',        str(params.outdir),
        '--runThreadN',       str(params.threads),
        '--genomeFastaFiles', str(input.fasta),
        '--sjdbGTFfile',      str(input.gtf),
        join_params("STAR", PARAMS, params.params_STAR),
        '2>',log.logfile
        ])

        # touch command for the star index
        command_touch = " ".join([
            params.touch, str(output.outfile),
            '2>',    str(log.logfile)
        ])

        command_final = command + ';' + command_touch
        print_shell(command_final)


# ----------------------------------------------------------------------------- #
# GIVEN PRIMARY AND SECONDARY GTF, COMBINES THEM INTO ONE GTF FILE
if GENOME_SECONDARY_IND:
    rule combine_gtf:
        input:
            primary   =  LINK_GTF_PRIMARY,
            secondary =  PATH_GTF_SECONDARY
        output:
            outfile = PATH_GTF_MIX
        params:
            threads = config['execution']['rules']['combine_gtf']['threads'],
            mem     = config['execution']['rules']['combine_gtf']['memory'],
            genome_name_primary   = GENOME_NAME_PRIMARY,
            genome_name_secondary = GENOME_NAME_SECONDARY,
            perl = SOFTWARE['perl']['executable'],
            perl_args = SOFTWARE['perl']['executable']['args']
        message:
            """
                Combining gtf files:
                    primary   : {input.primary}
                    secondary : {input.secondary}
                    output: {output.outfile}
            """
        shell:"""
            {params.perl} {params.perl_args} -pe 's|^|{params.genome_name_primary}|' < {input.primary} > {output.outfile}
            {params.perl} {params.perl_args} -pe 's|^|{params.genome_name_secondary}|' < {input.secondary} >> {output.outfile}
    """

# ----------------------------------------------------------------------------- #
def fetch_fastq_technical_replicates(wc):
    sample_name = str(wc.name)

    if(str(wc.type) == 'R1'):
        input_fastq = SAMPLE_SHEET.fetch_barcode_path(sample_name)

    if(str(wc.type) == 'R2'):
        input_fastq = SAMPLE_SHEET.fetch_reads_path(sample_name)

    if not input_fastq.__class__ == list:
        input_fastq = list(input_fastq)

    return(input_fastq)


rule merge_technical_replicates:
    input:
        infiles = fetch_fastq_technical_replicates
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", ("{name}" + "_{type,(R1)|(R2)}" + ".fastq.gz"))
    params:
        threads = config['execution']['rules']['merge_technical_replicates']['threads'],
        mem     = config['execution']['rules']['merge_technical_replicates']['memory'],
        tempdir = TEMPDIR,
        cat     = SOFTWARE['cat']['executable'],
        ln      = 'ln'
    log:
        log = os.path.join(PATH_LOG, '{name}.{type,(R1)|(R2)}.merge_technical_replicates.log')
    message:"""
            merge_technical_replicates:
                input:   {input}
                output : {output}
        """
    run:

        # if there is more than one input file, cat them
        # if there is only one, create a symbolic link
        if len(input.infiles) > 1:
            command = ' '.join([
                params.cat,
                " ".join(input.infiles),
                '>'  + str(output.outfile),
                '2>' + str(log.log)
            ])

        else:
            command = ' '.join([
                params.ln,
                str(input.infiles),
                str(output.outfile),
                '2>' + str(log.log)
            ])

        print_shell(command)


# ----------------------------------------------------------------------------- # filters reads based on quality, length and polyA
# uses flexbar to filter reads
def fetch_reads(wc):
    read_hash = {
        'barcode' : os.path.join(PATH_MAPPED, wc.name, SAMPLE_SHEET.fetch_field(wc.name,'barcode_merged')),
        'reads'   : os.path.join(PATH_MAPPED, wc.name, SAMPLE_SHEET.fetch_field(wc.name,'reads_merged'))
    }
    return(read_hash)

rule filter_reads:
    input:
        unpack(fetch_reads)
    output:
        barcode = os.path.join(PATH_MAPPED, "{name}", "{genome}","{name}_1.fastq.gz"),
        reads   = os.path.join(PATH_MAPPED, "{name}", "{genome}","{name}_2.fastq.gz")
    params:
        outpath       = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        name          = '{name}',
        threads       = config['execution']['rules']['filter_reads']['threads'],
        mem           = config['execution']['rules']['filter_reads']['memory'],
        tempdir       = TEMPDIR,
        flexbar       = SOFTWARE['flexbar']['executable'],

        # minimal base quality
        base_quality  = 20,

        # minimal length of the polyA homopolimer
        polya_length  = 10
    log:
        log = os.path.join(PATH_LOG, "{name}.{genome}.filter_reads.log")
    message:"""
        filter reads
                input:   {input}
                output reads   : {output.reads}
                output barcode : {output.barcode}
        """
    run:
        # calculates the minimal adapter size
        adapter_size = get_adapter_size(params.name)

        command = ' '.join([
            params.flexbar,
            '--reads',   str(input.barcode),
            '--reads2',  str(input.reads),
            '--target',  os.path.join(params.outpath, params.name),
            '--threads', str(params.threads),
            '--min-read-length', str(adapter_size),
            # quality trimming
            '--qtrim', 'TAIL',
            '--qtrim-format', 'i1.8',
            '--qtrim-threshold',  str(params.base_quality),
            # homopolyer trimming
            '--htrim-right', 'A',
            '--htrim-min-length', str(params.polya_length),
            '--zip-output', 'GZ'
        ])
        print_shell(command)


# ----------------------------------------------------------------------------- # Barcode histogram
# calculates the number of reads per cell
rule cell_barcode_histogram:
    input:
        infile = rules.filter_reads.output.barcode
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_cell_barcode_histogram.txt')
    params:
        outdir    = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname   = "{name}_{genome}",
        name      = "{name}",
        threads   = config['execution']['rules']['cell_barcode_histogram']['threads'],
        mem       = config['execution']['rules']['cell_barcode_histogram']['memory'],
        jellyfish = SOFTWARE['jellyfish']['executable'],
        perl      = SOFTWARE['perl']['executable'],
        zcat      = SOFTWARE['zcat']['executable'],
        tempdir   = TEMPDIR,
        hash_size = 10000000
    message: """
            cell_barcode_histogram:
                input:  {input.infile}
                output: {output.outfile}
        """
    log:
        logfile = os.path.join(PATH_LOG, '{name}_{genome}_cell_barcode_histogram.log')
    run:
        cb_adapter   = adapter_params(params.name, 'cell_barcode')
        count_file   = os.path.join(params.outdir, params.outname + '.jf')
        tmp_file     = os.path.join(params.outdir, params.outname + '.tmp.fastq')
        # parses the barcodes from the fastq file
        # counts the kmers
        # extracts the cell barcode kmers using a perl oneliner - not optimal but works
        # ------------------------------------------------ #
        perl_extract_cb = ' '.join([
            params.perl,
            '\'-ne if( $. % 4 == 0 | $. % 4 ==2 ){{print substr($_,',
             str(cb_adapter['start'] - 1),
             ',',
             str(cb_adapter['length']),
             ')."\\n";}}else{{print}}\''
        ])
        command_parse = ' '.join([
            params.zcat, input.infile, '|',
            perl_extract_cb, '>',
            tmp_file,
            '2>>', str(log.logfile)
        ])

        # ------------------------------------------------ #
        command_count = ' '.join([
            params.jellyfish, 'count',
            '-t',    str(params.threads),
            '-o',    str(count_file),
            '-m',    str(cb_adapter['length']),
            '-s',    str(params.hash_size),
            tmp_file,
            '2>>', str(log.logfile)
        ])

        # ------------------------------------------------ #
        # outputs the kmer table
        command_dump = ' '.join([
            params.jellyfish, 'dump',
            '--column',
            '--tab',
            '-o',    str(output.outfile),
            count_file,
            '2>>', str(log.logfile)
        ])

        # ------------------------------------------------ #
        # removes the jellyfish database
        command_remove = ' '.join([
            'rm', count_file, tmp_file
        ])

        command_final = ";".join([command_parse, command_count,command_dump,command_remove])
        print_shell(command_final)

# ----------------------------------------------------------------------------- #
# finds the barcode cutoff using inflection method
rule find_absolute_read_cutoff:
    input:
        infile = rules.cell_barcode_histogram.output.outfile
    output:
        outfile_yaml = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadCutoff.yaml'),
        outfile_tab  = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadCutoff.txt')
    params:
        outdir   = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname  = "{name}_{genome}",
        threads  = config['execution']['rules']['find_absolute_read_cutoff']['threads'],
        mem      = config['execution']['rules']['find_absolute_read_cutoff']['memory'],
        cutoff   = config['general']['cell_maximal_number'],
        script   = PATH_SCRIPT,
        Rscript  = PATH_RSCRIPT
    message: """
            find_absolute_read_cutoff:
                input:  {input.infile}
                output: {output.outfile_yaml}
        """
    run:
        RunRscript(input, output, params, params.script, 'Find_Absolute_Read_Cutoff.R')


# ----------------------------------------------------------------------------- #
# Maps single cell data using star and constructs the DGE matrix
rule map_star:
    input:
        barcode   = rules.filter_reads.output.barcode,
        reads     = rules.filter_reads.output.reads,
        genome    = rules.make_star_reference.output,
        whitelist = rules.find_absolute_read_cutoff.output.outfile_tab
    output:
        outfile   = os.path.join(PATH_MAPPED, "{name}", "{genome}","{name}_Aligned.out.bam")
    params:
        name        = "{name}",
        star        = SOFTWARE['STAR']['executable'],
        zcat        = SOFTWARE['zcat']['executable'],
        genome      = os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX'),
        outpath     = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        threads     = config['execution']['rules']['map_star']['threads'],
        mem         = config['execution']['rules']['map_star']['memory'],
        strand      = 'Forward',
        features    = " ".join(STAR_OUTPUT_TYPES),
        tempdir     = TEMPDIR,
        params_STAR = PARAMS['map_star']
    message:"""
        map_star:
            reads:   {input.reads}
            barcode: {input.barcode}
            output:  {output.outfile}
    """
    log:
        logfile = os.path.join(PATH_LOG, "{name}.{genome}.star.log")
    run:

        cb_adapter  = adapter_params(params.name, 'cell_barcode')
        umi_adapter = adapter_params(params.name, 'umi_barcode')
        infiles = " ".join([str(input.reads), str(input.barcode)])

        command = " ".join([
            params.star,
            '--genomeDir',  str(params.genome),
            '--runThreadN', str(params.threads),
            '--outFileNamePrefix', os.path.join(params.outpath, params.name) + '_',
            '--readFilesIn',     infiles,
            '--soloType',        'Droplet',
            '--soloCBwhitelist', str(input.whitelist),
            '--soloCBstart',     str(cb_adapter['start']),
            '--soloCBlen',       str(cb_adapter['length']),
            '--soloUMIstart',    str(umi_adapter['start']),
            '--soloUMIlen',      str(umi_adapter['length']),
            '--soloStrand',      str(params.strand),
            '--soloFeatures',    str(params.features),
            '--outSAMtype', 'BAM Unsorted',
            '--readFilesCommand', str(params.zcat),
            # CR/UR: raw (uncorrected) CellBarcode/UMI
            # CY/UY: quality score for CellBarcode/UMI
            # GX/GN: for gene ID/names
            # sS/sQ: for sequence/quality combined CellBarcode and UMI; sM for barcode match status.
            '--outSAMattributes NH HI nM AS CR UR GX GN sS sQ sM',
            join_params("STAR", PARAMS, params.params_STAR),
            '2>',str(log.logfile)
        ])
        print_shell(command)


# ----------------------------------------------------------------------------- #
rule sort_bam:
    input:
        infile = rules.map_star.output.outfile
    output:
        outfile   = os.path.join(PATH_MAPPED, "{name}", "{genome}","{name}.sorted.bam")
    params:
        samtools   = SOFTWARE['samtools']['executable'],
        threads    = config['execution']['rules']['sort_bam']['threads'],
        mem        = config['execution']['rules']['sort_bam']['memory'],
        tempdir    = TEMPDIR
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.sort_bam.log")
    message:"""
        sort_bam:
            input:  {input.infile}
            output: {output.outfile}
    """

    run:
        command = ' '.join([
            params.samtools, 'sort',
            '-@', str(params.threads),
            '-o', str(output.outfile),
            str(input.infile)
        ])
        print_shell(command)

# ----------------------------------------------------------------------------- #
rule index_bam:
    input:
        infile = rules.sort_bam.output.outfile
    output:
        outfile   = os.path.join(PATH_MAPPED, "{name}", "{genome}","{name}.sorted.bam.bai")
    params:
        samtools   = SOFTWARE['samtools']['executable'],
        threads    = config['execution']['rules']['sort_bam']['threads'],
        mem        = config['execution']['rules']['sort_bam']['memory'],
        tempdir    = TEMPDIR
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.index_bam.log")
    message:"""
        index_bam:
            input: {input.infile}
            output: {output.outfile}
    """

    run:
        command = ' '.join([
            params.samtools, 'index',
            '-b',
            '-@', str(params.threads),
            str(input.infile),
            str(output.outfile)
        ])
        print_shell(command)

# ----------------------------------------------------------------------------- #
# convert UMI matrix from txt format into one loom format
rule convert_matrix_from_mtx_to_loom:
    input:
        bamfile       = rules.sort_bam.output.outfile,
    output:
        outfile       = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.matrix.loom')
    params:
        name    = '{name}',
        python  = SOFTWARE['python']['executable'],
        threads = config['execution']['rules']['convert_matrix_from_mtx_to_loom']['threads'],
        mem     = config['execution']['rules']['convert_matrix_from_mtx_to_loom']['memory'],

        # input file
        indir   = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_Solo.out'),

        # STARsolo output types
        star_output_types_keys = " ".join(STAR_OUTPUT_TYPES_KEYS),
        star_output_types_vals = " ".join(STAR_OUTPUT_TYPES_VALS),

        # input gtf
        gtf               = lambda wildcards: os.path.join(PATH_ANNOTATION, wildcards.genome, '.'.join([wildcards.genome, 'gtf'])),
        script            = PATH_SCRIPT,
        sample_sheet_file = PATH_SAMPLE_SHEET
    log:
        logfile = os.path.join(PATH_LOG, "{name}.{genome}.convert_matrix_from_mtx_to_loom.log")
    message: """
            convert_matrix_from_mtx_to_loom:
                input:  {params.indir}
                output: {output.outfile}
        """
    run:
        command = ' '.join([
            params.python, os.path.join(params.script, 'convert_matrix_from_mtx_to_loom.py'),
            '--sample_id',              params.name,
            '--input_dir',              params.indir,
            '--gtf_file',               params.gtf,
            '--star_output_types_keys', params.star_output_types_keys,
            '--star_output_types_vals', params.star_output_types_vals,
            '--output_file',            output.outfile,
            '--sample_sheet_file',      params.sample_sheet_file,
            '--path_script',            params.script,
            '&>', str(log.logfile)
        ])
        print_shell(command)


# ----------------------------------------------------------------------------- #
## combines multiple loom files into one loom file
def fetch_loom_files(wc):
    loom_files  = expand(os.path.join(PATH_MAPPED, "{name}", wc.genome, '_'.join(["{name}", wc.genome, 'UMI.matrix.loom'])), name = SAMPLE_NAMES)
    return(loom_files)

rule combine_loom_files:
    input:
        infile   = fetch_loom_files
    output:
        outfile  = os.path.join(PATH_MAPPED, "{genome}_UMI.loom")
    params:
         python = SOFTWARE['python']['executable'],
         threads    = config['execution']['rules']['combine_loom_files']['threads'],
         mem        = config['execution']['rules']['combine_loom_files']['memory'],
         script = PATH_SCRIPT
    log:
        logfile = os.path.join(PATH_LOG, "{genome}.combine_loom_files.log")
    message: """
            combine_loom_files:
                input:  {input.infile}
                output: {output.outfile}
        """

    run:
        command = ' '.join([
            params.python, os.path.join(params.script, 'combine_loom_matrices.py'),
            '--input_files', " ".join(input.infile),
            '--output_file', output.outfile,
            '&>', str(log.logfile)
        ])
        print_shell(command)




# ----------------------------------------------------------------------------- #
## Imports and preprocesses the combined loom files and saves as SingleCellExperiment.RDS objects.
#  THIS function needs to be adapted if the pipeline will be extended to more than one genome
def fetch_gtf_path(wc):
    PATH_GTF = ''
    if config['annotation']['primary']['genome']['name'] == wc.genome:
        PATH_GTF = PATH_GTF_PRIMARY

    if 'secondary' in set(config['annotation'].keys()):
        if not config['annotation']['secondary'] == None:
            if config['annotation']['secondary']['genome']['name'] == wc.genome:
                PATH_GTF = PATH_GTF_SECONDARY

    if len(PATH_GTF) == 0:
        sys.exit('genome name is not properly defined')

    return(PATH_GTF)

rule convert_loom_to_singleCellExperiment:
    input:
        infile        = os.path.join(PATH_MAPPED, "{genome}_UMI.loom")
    output:
        outfile       = os.path.join(PATH_MAPPED, "{genome}.SingleCellExperiment.RDS")
    log:
        logfile = os.path.join(PATH_LOG, "{genome}.convert_loom_to_singleCellExperiment.log")
    params:
        gtf_file          = fetch_gtf_path,
        script            = PATH_SCRIPT,
        Rscript           = PATH_RSCRIPT,
        threads  = config['execution']['rules']['convert_loom_to_singleCellExperiment']['threads'],
        mem      = config['execution']['rules']['convert_loom_to_singleCellExperiment']['memory'],
        star_output_types_vals = " ".join(STAR_OUTPUT_TYPES_VALS)
    message: """
            convert_loom_to_singleCellExperiment:
                input:  {input.infile}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, params.script, 'convert_loom_to_singleCellExperiment.R')


# ----------------------------------------------------------------------------- #
rule convert_loom_to_seurat:
    input:
        infile        = os.path.join(PATH_MAPPED, "{genome}_UMI.loom")
    output:
        outfile       = os.path.join(PATH_MAPPED, "{genome}.Seurat.RDS")
    log:
        log = os.path.join(PATH_LOG, "{genome}.convert_loom_to_seurat.log")
    params:
        script            = PATH_SCRIPT,
        Rscript           = PATH_RSCRIPT,
        genome_version    = "{genome}",
        threads  = config['execution']['rules']['convert_loom_to_seurat']['threads'],
        mem      = config['execution']['rules']['convert_loom_to_seurat']['memory'],
        star_output_types_vals = " ".join(STAR_OUTPUT_TYPES_VALS)

    message: """
            convert_loom_to_seurat:
                input:  {input.infile}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, params.script, 'convert_loom_to_Seurat.R')


# ----------------------------------------------------------------------------- #
rule bam_to_BigWig:
    input:
        bamfile  = rules.sort_bam.output.outfile,
        index    = rules.index_bam.output
    output:
        bwfile   = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}.bw')
    params:
        threads = config['execution']['rules']['bam_to_BigWig']['threads'],
        mem     = config['execution']['rules']['bam_to_BigWig']['memory'],
        script  = PATH_SCRIPT,
        Rscript = PATH_RSCRIPT
    message: """
            bam_to_BigWig:
                input:  {input.bamfile}
                output: {output.bwfile}
            """
    run:
        RunRscript(input, output, params, params.script, 'BamToBigWig.R')

# ----------------------------------------------------------------------------- #
rule fastqc:
    input:
        infile = os.path.join(PATH_FASTQ, "{name}")
    output:
        outfile =  os.path.join(PATH_FASTQC, "{name}.fastqc.done")
    params:
        outpath = os.path.join(PATH_FASTQC),
        threads = config['execution']['rules']['fastqc']['threads'],
        mem     = config['execution']['rules']['fastqc']['memory'],
        java    = SOFTWARE['java']['executable'],
        fastqc  = SOFTWARE['fastqc']['executable'],
        touch   = SOFTWARE['touch']['executable']
    log:
        log = os.path.join(PATH_LOG, "{name}.fastqc.log")
    message: """
            fastqc:
                input: {input.infile}
                output: {output.outfile}
            """
    run:
        command = ' '.join([
            params.fastqc,
            '-j', params.java,
            '-t', params.threads,
            '-o', params.outpath,
            " ".join(input.infile),
            '&>', str(log.logfile)
        ])
        command_touch = ' '.join([
            params.touch, output.outfile
        ])
        command_final = ";".join([command, command_touch])
        print_shell(command_final)

# ----------------------------------------------------------------------------- #
rule extract_read_statistics:
    input:
        bamfile = rules.sort_bam.output.outfile,
        index   = rules.index_bam.output
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadStatistics.txt')
    params:
        sample            = "{name}",
        threads           = config['execution']['rules']['extract_read_statistics']['threads'],
        mem               = config['execution']['rules']['extract_read_statistics']['memory'],
        script            = PATH_SCRIPT,
        star_output_types = " ".join(STAR_OUTPUT_TYPES),
        Rscript           = PATH_RSCRIPT,
        mito_chr          = MITOCHNODRIAL_CHROMOSOME
    message: """
            extract_read_statistics:
                input:  {input.bamfile}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, params.script, 'Extract_Read_Statistics.R')


# ----------------------------------------------------------------------------- #
## Using the preprocessed SingleCellExperiment.RDS file, generates a self-contained HTML report
rule report:
    input:
        sceRds_file   = os.path.abspath(os.path.join(PATH_MAPPED, "{genome}.SingleCellExperiment.RDS")),
        read_stats    = READ_STATISTICS
    output:
        outfile       = os.path.join(PATH_MAPPED, "{genome}.scRNA-Seq.report.html")
    params:
        report_file   = os.path.join(PATH_SCRIPT, "scrnaReport.Rmd"),
        workdir       = os.path.abspath(PATH_MAPPED),
        Rscript       = PATH_RSCRIPT,
        script        = PATH_SCRIPT,
        path_mapped   = PATH_MAPPED,
    log:
        log = os.path.join(PATH_LOG, "{genome}.scRNA-Seq.report.log")
    message: """
            Generate an HTML report:
                input:  {input.sceRds_file}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, params.script, 'renderReport.R')
