# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Read_Statistics')


# -------------------------------------------------------------------------- #
.extractReadStatistics = function(
  bamfile,
  yieldSize = 1e6,
  norm.factor = 1e6,
  read_categories = c('INTRONIC','INTERGENIC','CODING','UTR')
  ){

    bam = scanBam(BamFile(bamfile), yieldSize=yieldSize,
                          param=ScanBamParam(
                            what='mapq',
                            tag=c('XF', 'GE', 'XC')
                          ))[[1]]

    default = setNames(rep(0, length(read_categories)),read_categories)
    non.ge.types = table(bam$tag$XF[which(bam$mapq == 255 & is.na(bam$tag$GE))])
    if(length(non.ge.types == 0)){
        non.ge.types = default

    }else{
        default[names(non.ge.types)] = non.ge.types
        non.ge.types = default
    }


    df = data.frame('total_reads'    = round(length(bam$mapq)/norm.factor, 4),
                    'uniq_mapped'    = round(sum(bam$mapq == 255, na.rm = T)/norm.factor, 1),
                    'intronic'       = round(non.ge.types['INTRONIC']/norm.factor, 4),
                    'intergenic'     = round(non.ge.types['INTERGENIC']/norm.factor, 4),
                    'coding'         = round(non.ge.types['CODING']/norm.factor, 4),
                    'UTR'            = round(non.ge.types['UTR']/norm.factor, 4),
                    'total_barcodes' = round(length(unique(bam$tag$XC[which(bam$mapq == 255)]))/norm.factor, 4))
    return (df)
  }

# -------------------------------------------------------------------------- #
Extract_Read_Statistics = function(
    bamfile = NULL,
    outfile = NULL,
    outname = 'Sample'

){
  if(is.null(bamfile))
    stop('bamfile not specified')

  if(is.null(outfile))
      stop('outfile not specified')

  suppressPackageStartupMessages(library(Rsamtools))

    message('Extracting read statistics ...')
      stat = suppressWarnings(.extractReadStatistics(bamfile))
      rownames(stat) = outname

    message('Writing read statistics ...')
      write.table(stat, outfile,
        row.names=TRUE, col.names=TRUE,
        sep='\t',quote=FALSE)
}

# -------------------------------------------------------------------------- #
Extract_Read_Statistics(
      bamfile = argv$input[['bamfile']],
      outfile = argv$output[['outfile']],
      outname = argv$params[['outname']]
  )
