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
MitoTrace <- function(bam_list = bams, 
                      fasta = fasta_loc, 
                      chr_name = "MT",
                      max_depth= 1e6, 
                      min_base_quality= 25, 
                      min_mapq = 30, 
                      min_nucleotide_depth = 0, 
                      min_minor_allele_depth = 0){
  require(Rsamtools)
  require(Matrix)
  require(seqinr)
  
  if(length(bam_list) == 1) {
    singlefile <- T
    }
  else{ 
    singlefile=F
    }
  
  if(singlefile) bam_list <- rep(bam_list, 2)
  
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
  # Define list of BAM files
  bam_name_list_array <- BamFileList(bam_list)

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
  pos_tmp <- unlist(lapply(pos_tmp, function(x) substr(x, 1, nchar(x) - 3)))
  pos_tmp <- as.numeric(pos_tmp)
  allpos <- allpos[order(pos_tmp)]
  matr <- matrix(0, length(allpos), ncol(res_counts2[["A"]]))
  rownames(matr) <- allpos
  colnames(matr) <- colnames(res_counts2[["A"]])
  lapply(res_counts2, function(x){
    ok <- intersect(rownames(matr), rownames(x))
    matr[ok,] <<- data.matrix(x[ok,])
  })
  
  counts <- matr
  
  counts <- as(counts, "sparseMatrix")
  coverage <- as(coverage, "sparseMatrix")
  
  if(singlefile){
    nom <- colnames(counts)[1]
    counts <- as.data.frame(counts[,1])
    coverage <- as.data.frame(coverage[,1])
    colnames(counts) <- colnames(coverage) <- nom
  }
    
  return(list(read_counts = counts, coverage = coverage))
}




# define the MitoTrace plot coverage depth function
MitoDepth <- function(mae = mae_res, species = "human", mt_ann = mt_ann){

  # set the gene label location
    if(species == "human"){
    location_human <- c(288.55, 577, 1124.55, 1602, 2450.05, 3230, 3784.55, 
                3900, 4300, 4700, 4990.55, 5150, 5550, 5950, 6350, 6750, 
                7150.55, 7300, 7718, 7927.55, 8295, 8469.05, 8867.05, 9598.55,
                9991, 10231.55, 10405, 10618.05, 11448.55, 11750, 12200, 12650, 
                13242.55, 14411.05, 14674, 15317.05, 15700, 16200, 16584.55)
    }
  
    # set color
    require("RColorBrewer")
    col= c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 11, name = "Spectral"), brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Accent"))

    # plot the coverave depth # single sample
    ymax <- max(log(mae[[2]][,1]))
    plot(row.names(mae[[2]]), log(mae[[2]][,1]), 
         type = "l", 
         col="brown1", 
         xlab="MT genome postion", 
         ylab="log2(Coverage depth)",
         main="scRNA-seq coverage depth cross MT genome", 
         ylim=c(-5,ymax))
    
    # plot the MT gene bar 
    for (i in 1:39) {
      if(grepl("TR", mt_ann$gene[i])){
        rect(mt_ann[i,2], -2.5, mt_ann[i,3], -1, col=col[i], border = "white")
        text(location[i], -0.5, mt_ann$gene[i], cex = 0.9)
      }
      else if(grepl("RN", mt_ann$gene[i])){
        rect(mt_ann[i,2], -2.5, mt_ann[i,3], -1, col=col[i], border = "white")
        text(location[i], -0.5, mt_ann$gene[i], cex = 0.9)
      }
      else{
        rect(mt_ann[i,2], -4, mt_ann[i,3], -2.5, col=col[i], border = "white")
        text(location[i], -4.5, mt_ann$gene[i], cex = 0.9)
      }
    }
    
  }


