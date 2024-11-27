

#' @title Get names of genes or peaks related to a query coordinates
#' @description Get names and coordinates of genes or peaks. It also get the
#' @description coordinates of query regions and returns the realted genes or
#' @description peak names.
#' @param names Names of genes or peaks
#' @param coordinates Coordinates of genes or peaks in granges format
#' @param queryCoordinates Coordinates of the query regions in granges format
#' @return Names of genes or peaks in queried reiongs

getNameFromCoordinates <- function(names, coordinates, queryCoordinates)
{
  queryCoordinates <- getCenterOfPeaks(queryCoordinates)
  grCoordinate <- string2Granges(coordinates)
  map <- GenomicRanges::findOverlaps(queryCoordinates, grCoordinate)
  mapQuery <- S4Vectors::queryHits(map)
  mapSubject <- S4Vectors::subjectHits(map)
  queryLen <- length(queryCoordinates)
  queryNames <- rep("", queryLen)
  for(i in seq_len(queryLen))
  {
    tmpInd <- which((mapQuery==i)==TRUE)
    if(length(tmpInd)>0)
    {
      indInTable <- mapSubject[tmpInd[1]]
      queryNames[i] <- names[indInTable]
    }
  }
  return(queryNames)
}


#' @title Remove reserved characters from a string
#' @description Remove reserved characters (such as *, +, -, etc) from a string
#' @param name A string of characters
#' @return A string without reserved characters
#' @examples
#' removeReserveCharacter("A&%B^f6")

removeReserveCharacter <- function(name)
{
  name <- gsub(" ", "", name, fixed = TRUE)
  name <- gsub(",", "", name, fixed = TRUE)
  name <- gsub("<", "", name, fixed = TRUE)
  name <- gsub(">", "", name, fixed = TRUE)
  name <- gsub("/", "", name, fixed = TRUE)
  name <- gsub("*", "", name, fixed = TRUE)
  name <- gsub("+", "", name, fixed = TRUE)
  name <- gsub("-", "", name, fixed = TRUE)
  name <- gsub("(", "", name, fixed = TRUE)
  name <- gsub(")", "", name, fixed = TRUE)
  name <- gsub("{", "", name, fixed = TRUE)
  name <- gsub("}", "", name, fixed = TRUE)
  name <- gsub("[", "", name, fixed = TRUE)
  name <- gsub("]", "", name, fixed = TRUE)
  name <- gsub("%", "", name, fixed = TRUE)
  name <- gsub("^", "", name, fixed = TRUE)
  name <- gsub("&", "", name, fixed = TRUE)
  name <- gsub("#", "", name, fixed = TRUE)
  name <- gsub("!", "", name, fixed = TRUE)
  name <- gsub("?", "", name, fixed = TRUE)
  name <- gsub(".", "", name, fixed = TRUE)
  name <- gsub("|", "", name, fixed = TRUE)
  name <- gsub("~", "", name, fixed = TRUE)
  name <- gsub("$", "", name, fixed = TRUE)
  name <- gsub("@", "", name, fixed = TRUE)
  name <- gsub(":", "", name, fixed = TRUE)
  name <- gsub(";", "", name, fixed = TRUE)
  return(name)
}



#' @title Generate genewise association between genes and peaks
#' @description Get genomic coordinates of a set of genes and a set of peaks
#' @description associate them by a fixed distance (default 50K nt). It also
#' @description associate genes and peaks for provided DNA-DNA interaction from
#' @description a dataset like HiC. This function can also associate genes and
#' @description user provided regions (ex. TADs, subTADs, etc). It generates
#' @description three tables: Gene table, peak table, and Gene-Peak association
#' @description table.
#' @param associationBy Can be "distance", "regions", or "DNAinteractions"
#' @param geneCoordinates Gene coordinates in granges format
#' @param geneNames Gene names can be provided by the user
#' @param peakCoordinates Peak coordinates in granges format
#' @param peakNames Peak names can be provided by the user
#' @param distance The maximum distance to associate peaks to genes. default 50K
#' @param givenRegions  granges coordinates of given regions (ex. TAD or loops)
#' @param strand1 granges of DNA strand1 linked to DNA strand2
#' @param strand2 granges of DNA strand2 linked to DNA strand1
#' @param outFile The name of the output folder (default "genewiseAssociation")
#' @return A vector of portions of linked genes and linked peaks
#' @examples
#' geneFile=system.file("extdata", "gene_expression.tsv", package="Site2Target")
#' geneCoords <- Table2Granges(geneFile)
#' geneTable <- read.table(geneFile, header=TRUE)
#'
#' geneDEIndices <- which((abs(geneTable$logFC)>1)==TRUE)
#' indicesLen <- length(geneDEIndices)
#' if(indicesLen >0)
#' {
#'     geneTable <- geneTable[geneDEIndices,]
#'     geneCoords <- geneCoords[geneDEIndices]
#' }
#' geneDENames <- geneTable$name
#' geneDElogFC <- geneTable$logFC
#' geneCoordsDE <- geneCoords
#'
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#' tfTable <- read.table(tfFile, header=TRUE)
#'
#' stats <-
#' genewiseAssociation(associationBy="distance",
#'                     geneCoordinates=geneCoordsDE,
#'                     geneNames=geneDENames,
#'                     peakCoordinates=TFCoords,
#'                     distance=50000,
#'                     outFile="Gene_TF_50K")
#' stats
#' @export



genewiseAssociation <-
  function(associationBy="distance", geneCoordinates=NULL, geneNames=NULL,
           peakCoordinates=NULL, peakNames = NULL, distance=50000,
           givenRegions=NULL, strand1=NULL, strand2=NULL,
           outFile="genewiseAssociation")
{
    geneNumber <- length(geneCoordinates)

    if(geneNumber==0)
    {
      stop("geneCoordinates must be provided in granges format")
    }

    peakNumber <- length(peakCoordinates)
    if(peakNumber==0)
    {
      stop("peakCoordinates must be provided in granges format")
    }

  peakCoordinateStr <- granges2String(peakCoordinates)
  if(length(peakNames)==0)
  {
    peakNames <- peakCoordinateStr
  } else {
    if(length(peakNames)!=peakNumber)
    {
      stop("peakNames should have same length as peakCoordinates")
    }
  }

  # Get center of TF sites and genes
  peakCoordinates <- getCenterOfPeaks(peakCoordinates)


  geneCoordinateStr <- granges2String(geneCoordinates)
  if(length(geneNames)==0)
  {
    geneNames  <-  geneCoordinateStr
  } else {
    if(length(geneNames)!=geneNumber)
    {
      stop("geneNames should have same length as gene Coordinates")
    }
  }


  if(associationBy=="distance")
  {
    extendRegions <-
      GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(GenomeInfoDb::seqnames(peakCoordinates)),
        ranges = IRanges::IRanges(BiocGenerics::start(peakCoordinates),
                                  end = BiocGenerics::end(peakCoordinates)
                                  ) + distance
      )
  } else if (associationBy=="regions")
  {
    givenRegionNumber <- length(givenRegions)
    if(givenRegionNumber==0)
    {
        stop("For extending sites in regions, the regions must be provided")
    }

    extendRegions <-
      extendSitesInGivenRegions(sites=peakCoordinates,
                                distance=distance,
                                givenRegions=givenRegions
      )

  } else if (associationBy=="DNAinteractions")
  {
    extendRegions <-
      GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(GenomeInfoDb::seqnames(peakCoordinates)),
        ranges = IRanges::IRanges(BiocGenerics::start(peakCoordinates),
                         end = BiocGenerics::end(peakCoordinates)
                         ) + distance
      )
  } else
  {
    stop("Peak to gene is associated either by distance or regions")
  }

  map <- GenomicRanges::findOverlaps(geneCoordinates, extendRegions)
  geneMap <- S4Vectors::queryHits(map)
  peakMap <- S4Vectors::subjectHits(map)
  geneCoverage <- length(unique(geneMap))/length(geneCoordinates)
  peakCoverage <- length(unique(peakMap))/length(peakCoordinates)

  tmpGeneNumb <- length(geneMap)
  if(tmpGeneNumb>0)
  {
      tmpGeneCoords <- geneCoordinates[geneMap]
      tmpPeakNumb <- length(peakMap)
      if(tmpPeakNumb>0)
      {
          tmpPeakCoords <- peakCoordinates[peakMap]

      }

      D <- site2GeneDistance(geneCoordinates=tmpGeneCoords,
                             peakCoordinates=tmpPeakCoords)

      # write link table
      df <- data.frame(
          geneNames=geneNames[geneMap],
          peak=peakNames[peakMap],
          distance=D
      )

  }






  if (associationBy=="DNAinteractions")
  {
    strandNumber <- length(strand1)

    if(strandNumber==0)
    {
      stop("Strand1 is empty")
    }
    strand2len <- length(strand2)

    if(strand2len!=strandNumber)
    {
      stop("Strand1 and Strand2 are not from the same length")
    }


    # Remove interactions lower than distance  ############# <----- This can become a function
    strand1Center <- getCenterOfPeaks(strand1)
    center1 <- BiocGenerics::start(strand1Center)
    rm(strand1Center)
    gc()

    strand2Center <- getCenterOfPeaks(strand2)
    center2 <- BiocGenerics::start(strand2Center)
    rm(strand2Center)
    gc()

    D <- abs(center1 - center2)

    distantInteractomInds <-
      which((D > (distance-1))==TRUE)

    strand1 <- strand1[distantInteractomInds]
    strand2 <- strand2[distantInteractomInds]

    ## Now add interactions to association table  ### This part can be a new function

    peakCoord <- getCenterOfPeaks(peakCoordinates)
    # Get TFsites in DNA inetarction strand2
    mapPeak <- GenomicRanges::findOverlaps(peakCoord, strand2)
    mapPeakInds <- S4Vectors::queryHits(mapPeak)
    mapPeakStrandInds <- S4Vectors::subjectHits(mapPeak)

    geneCoord <- getCenterOfPeaks(geneCoordinates)
    # Get GeneCorrdinates in DNA interactions strand1
    mapGene <- GenomicRanges::findOverlaps(geneCoord, strand1)
    mapGeneInds <- S4Vectors::queryHits(mapGene)
    mapGeneStrandInds <- S4Vectors::subjectHits(mapGene)

    # Take the gene-peak links which are mapped to a given link
    commonStrands <- intersect(mapPeakStrandInds, mapGeneStrandInds)
    commonStrandsNumber <- length(commonStrands)



    geneNamesDistal <- NULL
    peakDistal <- NULL
    distanceDistal <- NULL

    for(i in seq_len(commonStrandsNumber))
    {
      currentStrand <- commonStrands[i]

      strandIndsPeak <- which((mapPeakStrandInds==currentStrand)==TRUE)
      currentPeakInds <- mapPeakInds[strandIndsPeak]
      currentPeak <- peakNames[currentPeakInds]
      currentPeakCoord <- peakCoordinates[currentPeakInds]
      currentPeakNumber <- length(currentPeak)

      strandIndsGene <- which((mapGeneStrandInds==currentStrand)==TRUE)
      currentGeneInds <- mapGeneInds[strandIndsGene]
      currentGene <- geneNames[currentGeneInds]
      currentGeneCoord <- geneCoordinates[currentGeneInds]
      currentGeneNumber <- length(currentGene)

      # Take pairwise distances
      tmpDistance <- NULL
      for(j in seq_len(currentGeneNumber))
      {
          if(currentPeakNumber>0)
          {
              tmpGeneCoordinate <- rep(currentGeneCoord[j], currentPeakNumber)
              tmpDistance <- c(tmpDistance,
                               site2GeneDistance(
                                   geneCoordinates=tmpGeneCoordinate,
                                   peakCoordinates=currentPeakCoord
                               )
              )
          }
      }

      # Make repeated vector of Genes and peaks
      currentGene <- rep(currentGene, currentPeakNumber)

      # Make repeated vector of Genes and peaks
      tmpcurrentPeak <- NULL
      for(j in seq_len(currentGeneNumber))
      {
        tmpcurrentPeak <- c(tmpcurrentPeak, currentPeak)
      }
      currentPeak <- tmpcurrentPeak



      geneNamesDistal <- c(geneNamesDistal, currentGene)
      peakDistal <- c(peakDistal, currentPeak)
      distanceDistal <- c(distanceDistal , tmpDistance)


    }

    distalInteractionInds <- which((distanceDistal>distance)==TRUE)

    dfDistal <- data.frame(
      geneNames=geneNamesDistal[distalInteractionInds],
      peak=peakDistal[distalInteractionInds],
      distance=distanceDistal[distalInteractionInds]
    )

    # now combined two data frame
    geneNamesCombined <- c(df$geneNames, dfDistal$geneNames)
    peakCombined <- c(df$peak, dfDistal$peak)
    distanceCombined <- c(df$distance, dfDistal$distance)

    df <- data.frame(
      geneNames=geneNamesCombined,
      peak=peakCombined,
      distance=distanceCombined
    )

    geneCoverage <- length(unique(df$geneNames))/length(geneCoordinates)
    peakCoverage <- length(unique(df$peak))/length(peakCoordinates)


  }

  # sort table based on gene names
  tmpInds <-  order(df$geneNames)

  df <- df[tmpInds, ]

  if(!(dir.exists(outFile)))
  {
    dir.create(outFile)
  }




  utils::write.table(df, file=file.path(outFile, "link.tsv"), row.names=FALSE,
              col.names=TRUE,  quote=FALSE, sep="\t")


  #write gene table
  df <- data.frame(
    geneNames=geneNames,
    geneCoordinate=geneCoordinateStr
  )
  utils::write.table(df, file=file.path(outFile, "gene.tsv"), row.names=FALSE,
              col.names=TRUE,  quote=FALSE, sep="\t")



  #write peak table
  df <- data.frame(
    peakNames=peakNames,
    peakCoordinate=peakCoordinateStr
  )
  utils::write.table(df, file=file.path(outFile, "peak.tsv"), row.names=FALSE,
              col.names=TRUE,  quote=FALSE, sep="\t")


  stat <- data.frame(geneCoverage=geneCoverage,
                     peakCoverage=peakCoverage)

  return(stat)

}















#' @title Add column to gene-wise association
#' @description Add a column of values based on the type either genes or peaks.
#' @description The Input is either coordinates or names of genes or peaks plus
#' @description a column of relevant values. This function add these values as
#' @description a column to gene or peak table as well as the interaction table.
#' @param type type of columns to be added. Either "gene" or "peak"
#' @param name Names of genes or peaks
#' @param coordinates Coordinates of genes or peaks in granges format
#' @param columnName Column name that should be added to the tables
#' @param column Column values that should be added to the tables
#' @param inFile The name of the input folder (default "genewiseAssociation")
#' @param outFile The name of the output folder (default "genewiseAssociation")
#' @examples
#' geneFile=system.file("extdata", "gene_expression.tsv", package="Site2Target")
#' geneCoords <- Table2Granges(geneFile)
#' geneTable <- read.table(geneFile, header=TRUE)
#'
#' geneDEIndices <- which((abs(geneTable$logFC)>1)==TRUE)
#' indicesLen <- length(geneDEIndices)
#' if(indicesLen >0)
#' {
#'     geneTable <- geneTable[geneDEIndices,]
#'     geneCoords <- geneCoords[geneDEIndices]
#' }
#' geneDENames <- geneTable$name
#' geneDElogFC <- geneTable$logFC
#' geneCoordsDE <- geneCoords
#'
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#' tfTable <- read.table(tfFile, header=TRUE)
#' tfIntensities <- tfTable$intensities
#'
#' stats <-
#' genewiseAssociation(associationBy="distance",
#'                     geneCoordinates=geneCoordsDE,
#'                     geneNames=geneDENames,
#'                     peakCoordinates=TFCoords,
#'                     distance=50000,
#'                     outFile="Gene_TF_50K")
#' stats
#'
#' # add expression log fold changes to the table
#' addColumn2geneWiseAssociation(type="gene", name=geneDENames,
#'    columnName="Expr_logFC", column=geneDElogFC, inFile="Gene_TF_50K",
#'    outFile="Gene_TF_50K")
#'
#' # add peak intensitites to the table
#' addColumn2geneWiseAssociation(type="peak", coordinates=TFCoords,
#'    columnName="Binding_Intensities", column=tfIntensities,
#'    inFile="Gene_TF_50K", outFile="Gene_TF_50K")
#' @seealso
#' \code{\link{genewiseAssociation}}
#' @export

addColumn2geneWiseAssociation <-
  function(type="", name=NULL, coordinates=NULL, columnName=NA,
           column, inFile="geneWiseAssociation", outFile="geneWiseAssociation")
{

  if(is.na(columnName))
  {
    stop("Column name should be provided")
  }else{
    columnName <- removeReserveCharacter(columnName)

  }


  if (!dir.exists(inFile)){
    stop("The user provided directory does not exist")
  }

  # read interaction table
  if (!file.exists(file.path(inFile,"link.tsv"))){
    stop("The gene-peaak link file does not exist in the directory")
  }
  interactionTable <-
    utils::read.table(file.path(inFile,"link.tsv"), header=TRUE, sep = "\t")

  # read gene table
  if (!file.exists(file.path(inFile,"gene.tsv"))){
    stop("The gene information file does not exist in the directory")
  }
  geneTable <-
    utils::read.table(file.path(inFile,"gene.tsv"), header=TRUE, sep = "\t")

  if (!file.exists(file.path(inFile,"peak.tsv"))){
    stop("The peaak file does not exist in the directory")
  }
  peakTable <-
    utils::read.table(file.path(inFile,"peak.tsv"), header=TRUE, sep = "\t")

  lenName <- length(name)
  lenCoord <- length(coordinates)
  lenColumn <- length(column)

  if((lenName==0)&&(lenCoord==0))
  {
    stop("Either the names or coordinates must be provided by user")
  }else if((lenName>0)&&(lenCoord>0)) {
    stop("Either the names or coordinates must be provided by user")
  }else if((lenName>0)&&(lenCoord==0)) {
    dataType <- "name"
    if(lenName!=lenColumn)
    {
      stop("Provided names and columns must have the same size")
    }
  } else {
    dataType <- "coordinates"
    if(lenCoord!=lenColumn)
    {
      stop("Provided coordinates and columns must have the same size")
    }

  }

  if(type=="gene")
  {

    geneNumber <- dim(geneTable)[1]
    newCol <- rep("", geneNumber)

    if( dataType== "coordinates")
    {
      # find name of the genes and put them in name
      name <- getNameFromCoordinates(
        names=geneTable$geneNames,
        coordinates=geneTable$geneCoordinate,
        queryCoordinates=coordinates
      )
    }

    tmpInds <- match(name, geneTable$geneNames)
    mapNumber <- length(tmpInds)
    for(i in seq_len(mapNumber))
    {
      if(!is.na(tmpInds[i]))
      {
        if(newCol[tmpInds[i]]=="")
        {
          newCol[tmpInds[i]] <- as.character(column[i])
        } else {
          newCol[tmpInds[i]] <-  paste0(newCol[tmpInds[i]], " - ",
                                        as.character(column[i]))
        }
      }
    }

    # add new column to gene table
    eval(parse(text=
                 paste0("geneTable$", columnName, "<-newCol")))

    # add new column to interaction table
    map <- match(interactionTable$geneNames, geneTable$geneNames)

    ineratctionNumber <- dim(interactionTable)[1]
    newCol <- rep("", ineratctionNumber)

    tmpVec <-
      eval(parse(text=
                   paste0("geneTable$", columnName)))

    newCol <- tmpVec[map]

    eval(parse(text=
                 paste0("interactionTable$", columnName, "<-newCol")))


  } else if(type=="peak")
  {

    peakNumber <- dim(peakTable)[1]
    newCol <- rep("", peakNumber)

    if( dataType== "coordinates")
    {
      # find name of the peak and put them in name
      name <- getNameFromCoordinates(
        names=peakTable$peakName,
        coordinates=peakTable$peakCoordinate,
        queryCoordinates=coordinates
      )
    }

    tmpInds <- match(name, peakTable$peakName)
    mapNumber <- length(tmpInds)
    for(i in seq_len(mapNumber))
    {
      if(!is.na(tmpInds[i]))
      {
        if(newCol[tmpInds[i]]=="")
        {
          newCol[tmpInds[i]] <- as.character(column[i])
        } else {
          newCol[tmpInds[i]] <-  paste0(newCol[tmpInds[i]], " - ",
                                        as.character(column[i]))
        }
      }
    }

    # add new column to peak table
    eval(parse(text=
                 paste0("peakTable$", columnName, "<-newCol")))

    # add new column to interaction table
    map <- match(interactionTable$peak, peakTable$peakName)

    ineratctionNumber <- dim(interactionTable)[1]
    newCol <- rep("",ineratctionNumber)


    tmpVec <-
      eval(parse(text=
                   paste0("peakTable$", columnName)))

    newCol <- tmpVec[map]


    eval(parse(text=
                 paste0("interactionTable$", columnName, "<-newCol")))

  }else {
    stop("type should be provided as gene or peak")
  }


  if(!(dir.exists(outFile)))
  {
    dir.create(outFile)
  }


  # write link table
  utils::write.table(interactionTable, file=file.path(outFile, "link.tsv"),
              row.names=FALSE, col.names=TRUE,  quote=FALSE, sep="\t")


  #write gene table
  utils::write.table(geneTable, file=file.path(outFile, "gene.tsv"),
                     row.names=FALSE, col.names=TRUE,  quote=FALSE, sep="\t")


  #write peak table
  utils::write.table(peakTable, file=file.path(outFile, "peak.tsv"),
                     row.names=FALSE, col.names=TRUE,  quote=FALSE, sep="\t")


}






# Add relation into interaction table

#' @title Add a relation column to gene-peak interaction table
#' @description Get coordinates of interactions (ex. HiC interactions) and a
#' @description column of interaction values (ex. HiC intensities ) and add them
#' @description as a column to gene-peak interaction table.
#' @param strand1 granges of DNA strand1 linked to DNA strand2
#' @param strand2 granges of DNA strand2 linked to DNA strand1
#' @param columnName Column name that should be added to the interaction table
#' @param column Column values that should be added to the interaction table
#' @param inFile The name of the input folder (default "genewiseAssociation")
#' @param outFile The name of the output folder (default "genewiseAssociation")
#' @examples
#' geneFile=system.file("extdata", "gene_expression.tsv", package="Site2Target")
#' geneCoords <- Table2Granges(geneFile)
#' geneTable <- read.table(geneFile, header=TRUE)
#'
#' geneDEIndices <- which((abs(geneTable$logFC)>1)==TRUE)
#' indicesLen <- length(geneDEIndices)
#' if(indicesLen >0)
#' {
#'     geneTable <- geneTable[geneDEIndices,]
#'     geneCoords <- geneCoords[geneDEIndices]
#' }
#' geneDENames <- geneTable$name
#' geneDElogFC <- geneTable$logFC
#' geneCoordsDE <- geneCoords
#'
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#' tfTable <- read.table(tfFile, header=TRUE)
#'
#' stats <-
#' genewiseAssociation(associationBy="distance",
#'                     geneCoordinates=geneCoordsDE,
#'                     geneNames=geneDENames,
#'                     peakCoordinates=TFCoords,
#'                     distance=50000,
#'                     outFile="Gene_TF_50K")
#' stats
#'
#' HiCFile =system.file("extdata", "HiC_intensities.tsv", package="Site2Target")
#' HiCstr1 <- Table2Granges(HiCFile, chrColName="Strand1_chr",
#'                      startColName="Strand1_start", endColName="Strand1_end")
#' HiCstr2 <- Table2Granges(HiCFile, chrColName="Strand2_chr",
#'                      startColName="Strand2_start", endColName="Strand2_end")
#' HiCTable <- read.table(HiCFile, header=TRUE)
#' HiCintensities <- HiCTable$intensities
#'
#' addRelation2geneWiseAssociation(strand1=HiCstr1, strand2=HiCstr2,
#'      columnName="HiC_Intensities", column=HiCintensities,
#'      inFile="Gene_TF_50K", outFile="Gene_TF_50K")
#' @seealso
#' \code{\link{genewiseAssociation}}
#' @export


addRelation2geneWiseAssociation <-
  function(strand1=NULL, strand2=NULL, columnName, column,
           inFile="geneWiseAssociation", outFile="geneWiseAssociation")
{
  if(is.na(columnName))
  {
    stop("Column nam should be provided")
  }else{
    columnName <- removeReserveCharacter(columnName)

  }


  if (!dir.exists(inFile)){
    stop("The user provided directory does not exist")
  }

  # read interaction table
  if (!file.exists(file.path(inFile,"link.tsv"))){
    stop("The gene-peaak link file does not exist in the directory")
  }
  interactionTable <-
    utils::read.table(file.path(inFile,"link.tsv"), header=TRUE, sep = "\t")


  # read gene table
  if (!file.exists(file.path(inFile,"gene.tsv"))){
    stop("The gene information file does not exist in the directory")
  }
  geneTable <-
    utils::read.table(file.path(inFile,"gene.tsv"), header=TRUE, sep = "\t")

  geneCoord <- string2Granges(geneTable$geneCoordinate)
  geneCoord <- getCenterOfPeaks(geneCoord)

  if (!file.exists(file.path(inFile,"peak.tsv"))){
    stop("The peaak file does not exist in the directory")
  }
  peakTable <-
    utils::read.table(file.path(inFile,"peak.tsv"), header=TRUE, sep = "\t")
  peakCoord <- string2Granges(peakTable$peakCoordinate)
  peakCoord <- getCenterOfPeaks(peakCoord)


  lenColumn <- length(column)

  strandNumber <- length(strand1)

  if(strandNumber==0)
  {
    stop("strand1 and strand2 must be provided")
  }

  if(length(strand2)!=strandNumber)
  {
    stop("Strand1 and Strand2 should be from the same length")
  }

  if(lenColumn!=strandNumber)
  {
    stop("Column and Strands should be from the same length")
  }


  if(strandNumber > 0)
  {
    ineratctionNumber <- dim(interactionTable)[1]
    newCol <- rep("", ineratctionNumber)


    # Get TFsites in HiC strand1
    mapPeak <- GenomicRanges::findOverlaps(peakCoord, strand1)
    mapPeakInds <- S4Vectors::queryHits(mapPeak)
    mapPeakStrandInds <- S4Vectors::subjectHits(mapPeak)

    # Get GeneCorrdinates in HiC strand2
    mapGene <- GenomicRanges::findOverlaps(geneCoord, strand2)
    mapGeneInds <- S4Vectors::queryHits(mapGene)
    mapGeneStrandInds <- S4Vectors::subjectHits(mapGene)

    # Take the gene-peak links which are mapped to a given link
    commonStrands <- intersect(mapPeakStrandInds, mapGeneStrandInds)
    commonStrandsNumber <- length(commonStrands)

    # Add the given value to these peaks
    for(i in seq_len(commonStrandsNumber))
    {
      currentStrand <- commonStrands[i]

      strandIndsPeak <- which((mapPeakStrandInds==currentStrand)==TRUE)
      currentPeakInds <- mapPeakInds[strandIndsPeak]
      currentPeak <- peakTable$peakName[currentPeakInds]

      strandIndsGene <- which((mapGeneStrandInds==currentStrand)==TRUE)
      currentGeneInds <- mapGeneInds[strandIndsGene]
      currentGene <- geneTable$geneNames[currentGeneInds]

      ## Find the indices on the interactome table

      # Peak indices in the table
      tmpIndsPeaks <- which((interactionTable$peak==currentPeak[1]))
      currentPeaklen <- length(currentPeak)
      for(j in seq_len(currentPeaklen))
      {
        tmpIndsPeaks <-
          union(tmpIndsPeaks,
                which((interactionTable$peak==currentPeak[j]))
          )

      }

      # Gene indices in the table
      tmpIndsGenes <- which((interactionTable$geneNames==currentGene[1]))
      currentGeneslen <- length(currentGene)
      for(j in seq_len(currentGeneslen))
      {
        tmpIndsGenes <-
          union(tmpIndsGenes,
                which((interactionTable$geneNames==currentGene[j]))
          )

      }

      tmpInds <- intersect(tmpIndsGenes, tmpIndsPeaks)
      tmpIndsLen <- length(tmpInds)

      if(length(tmpIndsLen)>0)
      {
        newCol[tmpInds] <-  column[currentStrand]
      }


    }


  } else {
    stop("No valid link exist to be added to the interaction table")
  }


  # add new column to interaction table
  eval(parse(text=
               paste0("interactionTable$", columnName, "<-newCol")))

  if(!(dir.exists(outFile)))
  {
    dir.create(outFile)
  }


  # write link table
  utils::write.table(interactionTable, file=file.path(outFile, "link.tsv"),
              row.names=FALSE, col.names=TRUE,  quote=FALSE, sep="\t")


}


