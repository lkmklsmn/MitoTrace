#' @title  MitoTrace - an R package for the investigation of mitochondrial heteroplasmies.
#' 
#' @description  The MitoTrace function calculates mitochondrial heteroplasmies in (single-cell) RNA sequencing data based on the reads pileups of the mitochondrial genome.
#' 
#' @usage MitoTrace(bam_list = bams, ref_fasta = fasta_loc, name = "MT", max_depth = "", 
#'        min_base_quality = "", min_mapq = "", min_nucleotide_depth = "", min_minor_allele_depth = "")
#' 
#' @param bams_list Vector of absolute path(s) pointing to BAM alignment file(s).
#' @param ref_fasta Absolute path to the mitochondrial reference genome in FASTA format.
#' @param name Name of mitochondrial genome as specified in the BAM files. Sequence names can be check with the checkSequenceNames() function.
#' @param max_depth The maximum depth of reads considered at any position.
#' @param min_base_quality The minimum read base quality below which the base is ignored when summarizing pileup information.
#' @param min_mapq The minimum mapping quality below which the entire reads is ignored.
#' @param min_nucleotide_depth integer(1); minimum count of each nucleotide at a given position required for said nucleotide to appear in the result.
#' @param min_minor_allele_depth integer(1);  minimum count of all nucleotides other than the major allele at agiven position.
#' 
#' @details result <- MitoTrace(bam_list = bams, ref_fasta = fasta_loc, chr_name = "MT", max_depth = "1e6", min_base_quality=25, min_mapq=30, min_nucleotide_depth=0, min_minor_allele_depth=0)
#' @details See packageDescription("MitoTrace") for more details.
#' 
#' @note This package could not only apply for the analysis of single-cell data but also for bulk seuqencing data.
#' 
#' @return Read counts matrix and coverage matrix.
#' 
#' @author Mingqiang WANG <Mingqiang.Wang@uth.tmc.edu>, Simon Lukas <lkmklsmn@gmail.com>
#' 
#' @references The current source code of MitoTrace is from https://github.com/lkmklsmn/MitoTrace.


# check BAM chromosome names
checkSequenceNames <- function(bam_list){
  require(Rsamtools)
  tmp <- scanBamHeader(bam)
  targets <- names(tmp[[1]][[1]])
  text <- names(tmp[[1]][[2]])
  targets[which(text == "@SQ")]
}

# define the MitoTrace main function
MitoTrace <- function(bam_list = bams, fasta = fasta_loc, chr_name = "chrM",
                      max_depth= 1e6, min_base_quality=25, min_mapq=30, min_nucleotide_depth=0, min_minor_allele_depth=0){
  require(Rsamtools)
  require(Matrix)
  require(seqinr)
  
  bases <- c("A", "C", "G", "T")
  
  # Load in reference FASTA 
  reffasta <- seqinr::read.fasta(fasta)
  mitoChr = attr(reffasta, "name")
  maxpos = length(reffasta[[1]])
  reference_genome <- data.frame(
    postion = 1:maxpos,
    base = toupper(unname(reffasta)[[1]])[1:maxpos]
  )
  
  which <- GRanges(seqnames = chr_name, ranges = IRanges(1, maxpos))
  
  if(length(bam_list) == 1) combinedBam <- TRUE
  if(length(bam_list) > 1) combinedBam <- FALSE
  
  if(combinedBam){
    print("Extracting barcodes from single BAM file and running pileup for each barcode separately")
    
    params <- ScanBamParam(tag = "CB", which = gr)
    barcodes <- scanBam(bam, param = params)
    
    good_barcodes <- names(which(table(barcodes[[1]][[1]][[1]]) > 1000))
    
    # Run pileup command
    total_mpileups <- lapply(good_barcodes, function(x){
      
      filter <- list(x)
      names(filter) <- "CB"
      
      pileup_bam <- pileup(bam,
                           scanBamParam=ScanBamParam(tagFilter = l, which=gr),
                           pileupParam=PileupParam(distinguish_strands= FALSE, 
                                                   max_depth= max_depth, 
                                                   min_base_quality=min_base_quality, 
                                                   min_mapq=min_mapq, 
                                                   min_nucleotide_depth=min_nucleotide_depth, 
                                                   min_minor_allele_depth=min_minor_allele_depth
                           ))
      bases <- c("A", "T", "C", "G")
      base_counts <- lapply(bases, function(base){
        mutation <- subset(pileup_bam, nucleotide == base)
        data.frame(mutation$pos, mutation$count)
      })
      
      names(base_counts) <- bases
      base_counts
      
    })
    
    names(total_mpileups) <- good_barcodes
  }
  else{
    print("Running pileup on each BAM file separately")
    
    # Define list of BAM files
    bam_name_list_array <- BamFileList(bam_list)
    
    # Run pileup command
    total_mpileups <- lapply(bam_name_list_array, function(x){
      
      pileup_bam <- pileup(x,
                           scanBamParam=ScanBamParam(which=which),
                           pileupParam=PileupParam(distinguish_strands= FALSE, 
                                                   max_depth= max_depth, 
                                                   min_base_quality=min_base_quality, 
                                                   min_mapq=min_mapq, 
                                                   min_nucleotide_depth=min_nucleotide_depth, 
                                                   min_minor_allele_depth=min_minor_allele_depth
                           ))
      bases <- c("A", "T", "C", "G")
      base_counts <- lapply(bases, function(base){
        mutation <- subset(pileup_bam, nucleotide == base)
        data.frame(mutation$pos, mutation$count)
      })
      
      names(base_counts) <- bases
      base_counts
      
    })
    
    nom <- unlist(lapply(bam_list, basename))
    names(total_mpileups) <- nom
  }
  
  # Count number of bases at each position
  res_counts <- lapply(bases, function(base){
    
    allpos <- 1:maxpos
    
    counts_base_allcells <- do.call(cbind, lapply(total_mpileups, function(x){
      count_base_cell <- x[[base]]
      count_base_cell$mutation.count[match(allpos, count_base_cell$mutation.pos)]
    }))
    
    counts_base_allcells[which(is.na(counts_base_allcells))] <- 0
    rownames(counts_base_allcells) <- paste0(allpos, reference_genome$base, ">", base)
    colnames(counts_base_allcells) <- names(total_mpileups)
    as(counts_base_allcells, "sparseMatrix")
  })
  names(res_counts) <- bases
  
  # Calculate total read coverage at each position
  coverage <- matrix(0, nrow(res_counts[[1]]), ncol(res_counts[[1]]))
  rownames(coverage) <- as.character(1:maxpos)
  colnames(coverage) <- colnames(res_counts[[1]])
  lapply(1:4, function(x) coverage <<- coverage + data.matrix(res_counts[[x]]))
  
  # Remove reference allele sites
  res_counts2 <- lapply(bases, function(base){
    tmp <- res_counts[[base]]
    tmp[which(reference_genome$base != base),]
  })
  names(res_counts2) <- bases
  
  # Merge into one table and order
  allpos <- unique(unlist(lapply(res_counts2, rownames)))
  pos_tmp <- allpos
  subs <- c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")
  lapply(subs, function(x) pos_tmp <<- gsub(x, "", fixed = T, pos_tmp))
  pos_tmp <- as.numeric(pos_tmp)
  allpos <- allpos[order(pos_tmp)]
  matr <- matrix(0, length(allpos), ncol(res_counts2[["A"]]))
  rownames(matr) <- allpos
  colnames(matr) <- colnames(res_counts2[["A"]])
  lapply(res_counts, function(x){
    ok <- intersect(rownames(matr), rownames(x))
    matr[ok,] <<- data.matrix(x[ok,])
  })
  
  counts <- matr
  
  counts <- as(counts, "sparseMatrix")
  coverage <- as(coverage, "sparseMatrix")
  
  return(list(read_counts = counts, coverage = coverage))
}

