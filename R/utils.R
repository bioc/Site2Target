

#' @title Take Genomic Ranges from a table file
#' @description Read a table file and derive genomic ranges from user provided
#' @description column names.
#' @param fileName A table delimited file
#' @param chrColName Chromosomes column name (default: "Chr")
#' @param startColName Start column name (default: "start")
#' @param endColName End column name (default: "end")
#' @return granges format of given coordinates
#' @examples
#'
#' geneFile=system.file("extdata", "gene_expression.tsv", package="Site2Target")
#' grs <- Table2Granges(fileName=geneFile,
#'                       chrColName="chr",
#'                        startColName="start",
#'                        endColName="end")
#' grs
#'
#' @export

Table2Granges  <-
  function(fileName, chrColName="chr", startColName="start", endColName="end")
    {

  if (!(file.exists(fileName))){
    stop(fileName, " file does not exist")
  }

  Table <- utils::read.table(fileName, header=TRUE, stringsAsFactors=FALSE)
  colnamesTable <- colnames(Table)

  chrInd <- which((colnamesTable==chrColName))
  if (length(chrInd)==0)
  {
    stop("Provided chrColName column does not exist in the table")
  }

  startInd <- which((colnamesTable==startColName))
  if (length(startInd)==0)
  {
    stop("Provided startColName column does not exist in the table")
  }

  endInd <- which((colnamesTable==endColName))
  if (length(endInd)==0)
  {
    stop("Provided endColName column does not exist in the table")
  }

  granges <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(Table[,chrInd]),
    ranges=IRanges::IRanges(Table[,startInd], Table[,endInd]))
  rm(Table)
  gc()
  return(granges)
}



#' @title Return center of the given granges files
#' @description Get a granges and find the center of it
#' @param gr granges coordinate
#' @return granges format of the center
#' @examples
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#' TFCoordsCenters <- getCenterOfPeaks(TFCoords)
#' TFCoordsCenters
#'
#' @export

getCenterOfPeaks <- function(gr)
{
  chrs <- as.character(GenomeInfoDb::seqnames(gr))
  starts<- as.numeric(BiocGenerics::start(gr))
  ends<- as.numeric(BiocGenerics::end(gr))
  gc()
  centers <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(chrs),
    ranges =
      IRanges::IRanges(round((starts+ends)/2), end=round((starts+ends)/2))
  )
  return(centers)
}



#' @title Return the distance between paired peaks and genes
#' @description Get a granges of genes and peaks and return their distances
#' @param geneCoordinates granges coordinates of genes
#' @param peakCoordinates granges coordinates of peaks
#' @return the respective distances of paired genes and peaks

site2GeneDistance <- function(geneCoordinates, peakCoordinates)
{
  peakNumber <- length(peakCoordinates)
  geneNumber <- length(geneCoordinates)
  if(peakNumber!=geneNumber)
  {
    stop("site2GeneDistance requires paired genes and peaks")
  }

  peak <- getCenterOfPeaks(peakCoordinates)
  locPeaks <- BiocGenerics::start(peakCoordinates)

  startGenes <- BiocGenerics::start(geneCoordinates)
  endGenes <- BiocGenerics::end(geneCoordinates)
  geneLenghts <- endGenes - startGenes

  if(min(geneLenghts) < 0)
  {
    stop("site2GeneDistance requires correct gene coordinates")
  }


  distStarts <- abs(startGenes - locPeaks)
  distEnds <- abs(endGenes - locPeaks)

  # get the minimum distance of peaks to each ends of their paired genes
  D <- pmin(distStarts, distEnds)

  # put zero for peaks inside genes
  diffs <- abs(distEnds - distStarts) - geneLenghts
  insideGenesInds <- which((diffs==0)==FALSE)
  D[insideGenesInds] <- 0

  return(D)

}


#' @title Conver granes to strings of cooridnates
#' @description Get genomic coordinates granges and convert them to strings
#' @param gr granges coordinates
#' @return string of coordinates
#' @examples
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#' strCoords <- granges2String(TFCoords)
#' head(strCoords)
#'
#' @export

granges2String <- function(gr)
{
  chrs <- GenomeInfoDb::seqnames(gr)
  starts <- BiocGenerics::start(gr)
  ends <- BiocGenerics::end(gr)
  str <- paste0(chrs, ":", starts, "-", ends)
  return(str)
}



#' @title Conver strings to granges of cooridnates
#' @description Get genomic coordinates as trings and convert them to grangess
#' @param strCoordinates string of coordinates
#' @return Genomic coordinates in granges format
#' @examples
#' string2Granges(c("chr1:1112-1231", "ch2:3131-3221"))
#'
#' @export

string2Granges <- function(strCoordinates) {
  len <- length(strCoordinates)
  tmp <- unlist(strsplit(strCoordinates, ":"))
  chrs <- tmp[(c(1:len)*2-1)]
  Ranges <- tmp[(c(1:len)*2)]
  tmp <- unlist(strsplit(Ranges, "-"))
  start <-  tmp[(c(1:len)*2-1)]
  end <- tmp[(c(1:len)*2)]
  granges <-
    GenomicRanges::GRanges(
      seqnames=S4Vectors::Rle(chrs),
      ranges = IRanges::IRanges(as.numeric(start), as.numeric(end))
      )
  return(granges)

}

