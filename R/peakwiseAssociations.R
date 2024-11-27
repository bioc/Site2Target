



#' @title generate number of sites per gene given distances
#' @description Get genes and sites coordinates, and associate them by given
#' @description distance.
#' @param geneCoordinates granges coordinates of genes
#' @param sites granges coordinates of sites
#' @param distance the maximum distance to associate sites to genes. default 50K
#' @return A vector sites number matched to each gene
#' @examples
#' geneFile=system.file("extdata", "gene_expression.tsv", package="Site2Target")
#' geneCoords <- Table2Granges(geneFile)
#'
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#'
#' targetNum <- getTargetGenesNumber( geneCoords, TFCoords)

getTargetGenesNumber <- function(geneCoordinates=NA, sites=NA, distance=50000)
{

  geneNumber <- length(geneCoordinates)
  if(geneNumber<2)
  {
    stop("At least two genes corrdinats must be provided")
  }

  siteNumber <- length(sites)
  if( siteNumber<2)
  {
    stop("At least two sites must be provided")
  }

  extendRegions <-
    GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(GenomeInfoDb::seqnames(sites)),
      ranges = IRanges::IRanges(BiocGenerics::start(sites),
                                end = BiocGenerics::end(sites)
                                ) + distance
      )

  targets <-
    S4Vectors::queryHits(
      GenomicRanges::findOverlaps(geneCoordinates, extendRegions)
      )

  tmpDF <- as.data.frame(table(targets), stringsAsFactors = FALSE)

  targetNumber <- rep(0, geneNumber)

  for(i in seq_len(geneNumber))
  {
    targetNumber[as.numeric(tmpDF$targets[i])] <- as.numeric(tmpDF$Freq[i])
  }

  return(targetNumber)
}




#' @title Extend sites given regions boundaries
#' @description Get sites and given regions (ex. TADs or loops) coordinates.
#' @description It extends sites in a give region using a distance function
#' @param givenRegions  granges coordinates of given regions (ex. TAD or loops)
#' @param sites granges coordinates of sites
#' @param distance the maximum distance to associate sites to regions
#' @return A granges of the extended sites in given regions

#' @examples
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#'
#' TADsFile =system.file("extdata", "TADs.tsv",package="Site2Target")
#' TADs <- Table2Granges(TADsFile)
#'
#' extendSitesInGivenRegions(TADs, TFCoords)


extendSitesInGivenRegions <- function(givenRegions, sites, distance=100000)
{
  # Entend +/- distance around sites
  extendRegions <-
    GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(GenomeInfoDb::seqnames(sites)),
      ranges = IRanges::IRanges(BiocGenerics::start(sites),
                                end = BiocGenerics::end(sites)
                                ) + distance
      )

  # Find map of sites and given regions
  siteRegionOverlap <- GenomicRanges::findOverlaps(sites, givenRegions)
  sitesInRegions <- S4Vectors::queryHits(siteRegionOverlap)

  # Remove sites out of given regions
  extendRegions <- extendRegions[sitesInRegions]

  # Get given regions coordinate of sites
  regionsOfSites <- givenRegions[S4Vectors::subjectHits(siteRegionOverlap)]


  # keep the extension of sites inside the given regions
  Starts <-
    pmax(BiocGenerics::start(extendRegions),BiocGenerics::start(regionsOfSites))

  Ends <-
    pmin(BiocGenerics::end(extendRegions),BiocGenerics::end(regionsOfSites))

  extendRegionsGiven <-
    GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(GenomeInfoDb::seqnames(extendRegions)),
            ranges = IRanges::IRanges(Starts, end=Ends)
      )

  return(extendRegionsGiven)

}




#' @title Fit Negative binomial distribution to target genes
#' @description Get genes and sites coordinates, and associate them by given
#' @description distance or given regions (ex. TADs or loops). It tests the
#' @description distribution of sites around genes either by poisson or
#' @description negative binomial test.
#' @param geneCoordinates granges coordinates of genes
#' @param sites granges coordinates of sites
#' @param distance the maximum distance to associate sites to genes. default 50K
#' @param dist either "negative binomial" or "poisson"
#' @param associationBy either "distance" or "regions"
#' @param givenRegions user provided granges regions like TADs or loops
#' @return A vector of pvalue distribution for target genes
#' @examples
#'
#' geneFile=system.file("extdata", "gene_expression.tsv", package="Site2Target")
#' geneCoords <- Table2Granges(geneFile)
#'
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#'
#' pvals <- getTargetGenesPvals( geneCoordinates=geneCoords, sites=TFCoords)
#'
#' @export

getTargetGenesPvals <-
  function(associationBy="distance", dist= "negative binomial",
           geneCoordinates=NA, sites=NA, distance=50000, givenRegions=NA)
{

    geneNumber <- length(geneCoordinates)
    if(geneNumber<2)
    {
      stop("At least two genes corrdinats must be provided")
    }

    siteNumber <- length(sites)
    if( siteNumber<2)
    {
      stop("At least two sites must be provided")
    }

  sites <- getCenterOfPeaks(sites)
  if(associationBy=="distance")
  {
    targetNumber <-
      getTargetGenesNumber(geneCoordinates=geneCoordinates,
                           sites=sites,
                           distance=distance
                           )

  } else if (associationBy=="regions")
  {
    givenRegionNumber <- length(givenRegions)
    if(givenRegionNumber<2)
    {
      if(is.na(givenRegions))
      {
        stop("For extending sites in regions, the regions must be provided")
      }
    }
    extendRegions <-
      extendSitesInGivenRegions(
        sites=sites,
        distance=distance,
        givenRegions=givenRegions
      )


    targetNumber <-
      getTargetGenesNumber(geneCoordinates=geneCoordinates,
                           sites=extendRegions,
                           distance=0
      )


  } else
  {
    stop("Peak to gene is associated either by distance or regions")
  }




  # Don't consider high binding sites to find distribution
  eps <- 1
  log2ScaleCount <- log2(targetNumber + eps)

  upperbound <-
    2^(ceiling(stats::quantile(log2ScaleCount ,0.75)+
                 3*stats::IQR(log2ScaleCount)))

  if(upperbound<4)
  {
    warning("Insufficeint interactions to model")
    acceptedInds <- c(1:geneNumber)
  }
  else
  {
    acceptedInds <- which((targetNumber<upperbound)==TRUE)
  }



  if(dist=="negative binomial")
  {
    flag <- TRUE
    try({
      distNB <- MASS::fitdistr(targetNumber[acceptedInds],
                               densfun = "negative binomial")
      pvals <- stats::pnbinom(targetNumber,
                              size=(distNB$estimate)[1],
                              mu=(distNB$estimate)[2],
                              lower.tail = FALSE
                              )
      flag <- FALSE
    })

    if(flag)
    {
      stop("negative binomial distribution could not be fitted try poisson")
    }


  } else if(dist=="poisson"){
    distP <- MASS::fitdistr(targetNumber[acceptedInds], densfun = "poisson")
    pvals <- stats::ppois(targetNumber,
                          lambda = as.numeric(distP[1]),
                          lower.tail = FALSE
                          )

  } else
  {
    stop("The distribution should be either negative binomial or poisson")
  }
  return(pvals)
}






#' @title Fit Negative binomial distribution to target genes
#' @description Get genes and sites coordinates, and associate them by given
#' @description distance and user provided DNA interaction (ex. HiC). It tests
#' @description the distribution of sites around genes either by poisson or
#' @description negative binomial test.
#' @param geneCoordinates granges coordinates of genes
#' @param sites granges coordinates of sites
#' @param distance the maximum distance to associate sites to genes. default 50K
#' @param dist either "negative binomial" or "poisson"
#' @param strand1 granges of DNA strand1 linked to DNA strand2
#' @param strand2 granges of DNA strand2 linked to DNA strand1
#' @return A vector of pvalue distribution for target genes
#' @examples
#'
#' geneFile=system.file("extdata", "gene_expression.tsv", package="Site2Target")
#' geneCoords <- Table2Granges(geneFile)
#'
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#'
#' HiCFile =system.file("extdata", "HiC_intensities.tsv", package="Site2Target")
#' HiCstr1 <- Table2Granges(HiCFile, chrColName="Strand1_chr",
#'                      startColName="Strand1_start", endColName="Strand1_end")
#' HiCstr2 <- Table2Granges(HiCFile, chrColName="Strand2_chr",
#'                      startColName="Strand2_start", endColName="Strand2_end")
#'
#' pvals <- getTargetGenesPvalsWithDNAInteractions(
#'                geneCoordinates=geneCoords, sites=TFCoords, strand1=HiCstr1,
#'                strand2=HiCstr2)
#'
#' @export

getTargetGenesPvalsWithDNAInteractions <-
  function(dist= "negative binomial", geneCoordinates=NA, sites=NA,
           strand1=NA, strand2=NA, distance=50000)
{

    geneNumber <- length(geneCoordinates)
    if(geneNumber<2)
    {
      stop("At least two genes corrdinats must be provided")
    }

    siteNumber <- length(sites)
    if( siteNumber<2)
    {
      stop("At least two sites must be provided")
    }

    LenStrand1 <- length(strand1)
    LenStrand2 <- length(strand2)

    if(LenStrand1<2)
    {
      stop("At least two DNA-DNA interactions must be provided")
    }


    if(LenStrand1!=LenStrand2)
    {
      stop("The length of Gstrand and Sstrand must be equal")
    }

  # Get center of TF sites and genes
  sites <- getCenterOfPeaks(sites)

  targetNumber <- rep(0, geneNumber)

  # First add interactions by distance

  if(distance > -1)
  {
    # Remove interactions lower than distance
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

    InteractionNumber <- length(distantInteractomInds)
    if(InteractionNumber > 0)
    {
      strand1 <- strand1[distantInteractomInds]
      strand2 <- strand2[distantInteractomInds]
    }

    targetNumber <-
      getTargetGenesNumber(geneCoordinates=geneCoordinates,
                           sites=sites,
                           distance=distance
      )

  }


  # Add the user give interaction

  geneCoordinates <- getCenterOfPeaks(geneCoordinates)

  InteractionNumber <- LenStrand1

    # Get sites in strand1
    mapSite <- GenomicRanges::findOverlaps(sites, strand1)
    mapSiteInds <- S4Vectors::queryHits(mapSite)
    mapSiteStrandInds <- S4Vectors::subjectHits(mapSite)

    # Get geneCoordinates in strand2
    mapGene <- GenomicRanges::findOverlaps(geneCoordinates,strand2)
    mapGeneInds <- S4Vectors::queryHits(mapGene)
    mapGeneStrandInds <- S4Vectors::subjectHits(mapGene)


    commonStrands <- intersect(mapSiteStrandInds, mapGeneStrandInds)

    targetNumberDNAIntacts <- rep(0, geneNumber)


    for(i in seq_len(geneNumber))
    {
      # detect gene DNA interactome
      tmpGeneInds <- which((mapGeneInds==i)==TRUE)
      tmpGeneIndsNum <- length(tmpGeneInds)

      if(tmpGeneIndsNum > 0)
      {
        # Take the HiC indices in gene mapping
        tmpStrandInds <- mapGeneStrandInds[tmpGeneInds]

        intersect(tmpStrandInds, commonStrands)

        # Take the HiC indices in TF mapping
        tmpSiteInds <- match(tmpStrandInds, mapSiteStrandInds)
        tmpSiteInds <- tmpSiteInds[!is.na(tmpSiteInds)]

        targetNumberDNAIntacts[i] <-
        targetNumberDNAIntacts[i] + length(tmpSiteInds)
      }

    }
    targetNumber <- targetNumber + targetNumberDNAIntacts



  # Don't consider high binding sites to find distribution
  eps <- 1
  log2ScaleCount <- log2(targetNumber + eps)

  upperbound <-
    2^(ceiling(stats::quantile(log2ScaleCount ,0.75)+
                 3*stats::IQR(log2ScaleCount)))

  if(upperbound<4)
  {
    warning("Insufficeint interactions to model")
    acceptedInds <- c(1:geneNumber)
  }
  else
  {
    acceptedInds <- which((targetNumber<upperbound)==TRUE)
  }


  if(dist=="negative binomial")
  {
    flag <- TRUE
    try({
      distNB <- MASS::fitdistr(targetNumber[acceptedInds],
                               densfun = "negative binomial")
      pvals <- stats::pnbinom(targetNumber,
                              size=(distNB$estimate)[1],
                              mu=(distNB$estimate)[2],
                              lower.tail = FALSE
      )
      flag <- FALSE
    })

    if(flag)
    {
      stop("negative binomial distribution could not be fitted try poisson")
    }


  } else if(dist=="poisson"){
    distP <- MASS::fitdistr(targetNumber[acceptedInds], densfun = "poisson")
    pvals <- stats::ppois(targetNumber,
                          lambda = as.numeric(distP[1]),
                          lower.tail = FALSE
    )

  } else
  {
    stop("The distribution should be either negative binomial or poisson")
  }
  return(pvals)
}



#' @title Fit log-normal distribution to target genes
#' @description Get genes and sites coordinates, and associate them by given
#' @description distance or given regions (ex. TADs or loops). It tests the
#' @description distribution of log-intensities of sites around genes by
#' @description log-normal test. This function consider both binding sites and
#' @description intensities.
#' @param geneCoordinates granges coordinates of genes
#' @param sites granges coordinates of sites
#' @param intensities intensity values associated to sites
#' @param distance the maximum distance to associate sites to genes. default 50K
#' @param associationBy either "distance" or "regions"
#' @param givenRegions user provided granges regions like TADs or loops
#' @return A vector of pvalue distribution for target genes
#' @examples
#'
#' geneFile=system.file("extdata", "gene_expression.tsv", package="Site2Target")
#' geneCoords <- Table2Granges(geneFile)
#'
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#' tfTable <- read.table(tfFile, header=TRUE)
#' tfIntensities <- tfTable$intensities
#'
#' pvals <- getTargetGenesPvalsWithIntensities(geneCoordinates=geneCoords,
#'                       sites=TFCoords, intensities=tfIntensities)
#'
#' @export

getTargetGenesPvalsWithIntensities <-
  function(associationBy="distance", intensities,
             geneCoordinates=NA, sites=NA, distance=50000, givenRegions=NA)

{

    geneNumber <- length(geneCoordinates)
    if(geneNumber<10)
    {
      stop("At least ten genes corrdinats must be provided")
    }

    siteNumber <- length(sites)
    if( siteNumber<10)
    {
      stop("At least ten sites must be provided")
    }


  sites <- getCenterOfPeaks(sites)
  if(associationBy=="distance")
  {
    extendRegions <-
      GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(GenomeInfoDb::seqnames(sites)),
        ranges = IRanges::IRanges(BiocGenerics::start(sites),
                                  end = BiocGenerics::end(sites)
        ) + distance
      )

  } else if (associationBy=="regions")
  {

    givenRegionNumber <- length(givenRegions)
    if(givenRegionNumber<2)
    {
      if(is.na(givenRegions))
      {
        stop("For extending sites in regions, the regions must be provided")
      }
    }

    extendRegions <-
      extendSitesInGivenRegions(sites=sites,
                                distance=distance,
                                givenRegions=givenRegions
                                )
  } else
  {
    stop("Peak to gene is associated either by distance or regions")
  }

  overlaps <- GenomicRanges::findOverlaps(geneCoordinates, extendRegions)
  targetNumber <- rep(0, geneNumber)
  overlapsNumber <- length(overlaps)

  if(overlapsNumber > 10)
  {
    for(i in seq_len(overlapsNumber))
    {
      targetNumber[S4Vectors::queryHits(overlaps[i])] <-
        targetNumber[S4Vectors::queryHits(overlaps[i])] +
        intensities[S4Vectors::subjectHits(overlaps[i])]

    }
  }
  else
  {
    stop("Genes and sites are far from each other")
  }

  # Don't consider high binding sites to find distribution
  eps <- 1
  log2ScaleCount <- log2(targetNumber + eps)

  nonZeroInds <- which((log2ScaleCount>0)==TRUE)
  lowerbound <-
    ceiling(stats::quantile(log2ScaleCount[nonZeroInds] ,0.25)-
              1.5*stats::IQR(log2ScaleCount[nonZeroInds]))

  upperbound <-
    ceiling(stats::quantile(log2ScaleCount[nonZeroInds] ,0.75)+
              1.5*stats::IQR(log2ScaleCount[nonZeroInds]))

  acceptedInds <- intersect(which((log2ScaleCount<upperbound)==TRUE),
                            which((log2ScaleCount>lowerbound)==TRUE) )

  # try to fit log normal distribution without using outliers
  flag <- TRUE
  try({
    distN <- MASS::fitdistr(log2ScaleCount[acceptedInds], densfun = "normal")
    pvals <- stats::pnorm(log2ScaleCount,
                          mean=(distN$estimate)[1],
                          sd=(distN$estimate)[2],
                          lower.tail = FALSE
                          )
    flag <- FALSE
  })

  if(flag)
  {
    warning("Low number of sites and genes")

    # Try to fit the distribution with considering outliers due to low samples
    try({
      distN <- MASS::fitdistr(log2ScaleCount, densfun = "normal")
      pvals <- stats::pnorm(log2ScaleCount,
                            mean=(distN$estimate)[1],
                            sd=(distN$estimate)[2],
                            lower.tail = FALSE
      )
      flag <- FALSE
    })

  }

  if(flag)
  {
    stop("Cannot fit the log-normal distirbution. Use more sites and genes")
  }


  return(pvals)

}




#' @title Fit log-normal distribution to target genes
#' @description Get genes and sites coordinates, and associate them by given
#' @description distance and user provided DNA interaction (ex. HiC). It tests
#' @description the distribution of log-intensities of sites around genes by
#' @description log-normal test. This function consider both binding sites and
#' @description intensities.
#' @param geneCoordinates granges coordinates of genes
#' @param sites granges coordinates of sites
#' @param intensities intensity values associated to sites
#' @param distance the maximum distance to associate sites to genes. default 50K
#' @param strand1 granges of DNA strand1 linked to DNA strand2
#' @param strand2 granges of DNA strand2 linked to DNA strand1
#' @return A vector of pvalue distribution for target genes
#' @examples
#'
#' geneFile=system.file("extdata", "gene_expression.tsv", package="Site2Target")
#' geneCoords <- Table2Granges(geneFile)
#'
#' tfFile =system.file("extdata", "MEIS_binding.tsv", package="Site2Target")
#' TFCoords <- Table2Granges(tfFile)
#' tfTable <- read.table(tfFile, header=TRUE)
#' tfIntensities <- tfTable$intensities
#'
#' HiCFile =system.file("extdata", "HiC_intensities.tsv", package="Site2Target")
#' HiCstr1 <- Table2Granges(HiCFile, chrColName="Strand1_chr",
#'                      startColName="Strand1_start", endColName="Strand1_end")
#' HiCstr2 <- Table2Granges(HiCFile, chrColName="Strand2_chr",
#'                      startColName="Strand2_start", endColName="Strand2_end")
#'
#' pvals <- getTargetGenesPvalsWithIntensitiesAndDNAInteractions(
#'                        geneCoordinates=geneCoords, sites=TFCoords,
#'                        intensities=tfIntensities, strand1=HiCstr1,
#'                         strand2=HiCstr2)
#'
#' @export

getTargetGenesPvalsWithIntensitiesAndDNAInteractions <-
  function(geneCoordinates, sites, intensities, strand1, strand2,
           distance=50000)
{
    geneNumber <- length(geneCoordinates)
    if(geneNumber<10)
    {
      stop("At least ten genes corrdinats must be provided")
    }

    siteNumber <- length(sites)
    if( siteNumber<10)
    {
      stop("At least ten sites must be provided")
    }

    LenStrand1 <- length(strand1)
    LenStrand2 <- length(strand2)

    if(LenStrand1<2)
    {
      stop("At least two DNA-DNA interactions must be provided")
    }


    if(LenStrand1!=LenStrand2)
    {
      stop("The length of Gstrand and Sstrand must be equal")
    }

    sites <- getCenterOfPeaks(sites)
    targetNumber <- rep(0, geneNumber)

  if(distance > -1)
  {
    # Remove interactions lower than distance
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

    InteractionNumber <- length(distantInteractomInds)
    if(InteractionNumber > 0)
    {
      strand1 <- strand1[distantInteractomInds]
      strand2 <- strand2[distantInteractomInds]
    }


    # Add nearby links: up to user provided distance
    extendRegions <-
      GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(GenomeInfoDb::seqnames(sites)),
        ranges = IRanges::IRanges(BiocGenerics::start(sites),
                                  end = BiocGenerics::end(sites)
        ) + distance
      )

    overlaps <- GenomicRanges::findOverlaps(geneCoordinates, extendRegions)

    overlapsNumber <- length(overlaps)

    if(overlapsNumber > 10)
    {
      for(i in seq_len(overlapsNumber))
      {
        targetNumber[S4Vectors::queryHits(overlaps[i])] <-
          targetNumber[S4Vectors::queryHits(overlaps[i])] +
          intensities[S4Vectors::subjectHits(overlaps[i])]

      }
    }
    else
    {
      stop("Genes and sites are far from each other")
    }


  }

  geneCoordinates <- getCenterOfPeaks(geneCoordinates)

  # Add interaction links
  InteractionNumber <- length(strand1)

    # Get sites in strand1
    mapSite <- GenomicRanges::findOverlaps(sites, strand1)
    mapSiteInds <- S4Vectors::queryHits(mapSite)
    mapSiteStrandInds <- S4Vectors::subjectHits(mapSite)

    # Get geneCoordinates in strand2
    mapGene <- GenomicRanges::findOverlaps(geneCoordinates,strand2)
    mapGeneInds <- S4Vectors::queryHits(mapGene)
    mapGeneStrandInds <- S4Vectors::subjectHits(mapGene)

    commonStrands <- intersect(mapSiteStrandInds, mapGeneStrandInds)
    targetNumberDNAIntacts <- rep(0, geneNumber)

    for(i in seq_len(geneNumber))
    {
      # detect gene DNA interactome
      tmpGeneInds <- which((mapGeneInds==i)==TRUE)
      tmpGeneIndsNumber <- length(tmpGeneInds)
      if(tmpGeneIndsNumber > 0)
      {
        # Take the HiC indices in gene mapping
        tmpStrandInds <- mapGeneStrandInds[tmpGeneInds]

        intersect(tmpStrandInds, commonStrands)

        # Take the HiC indices in TF mapping
        tmpSiteInds <- match(tmpStrandInds, mapSiteStrandInds)
        tmpSiteInds <- tmpSiteInds[!is.na(tmpSiteInds)]
        tmpSiteNumber <- length(tmpSiteInds)
        for(siteCounter in seq_len(tmpSiteNumber))
        {
          targetNumberDNAIntacts[i] <- targetNumberDNAIntacts[i] +
            intensities[mapSiteInds[tmpSiteInds[siteCounter]]]
        }

      }

    }
  targetNumber <- targetNumber + targetNumberDNAIntacts


  eps <- 1
  log2ScaleCount <- log2(targetNumber + eps)

  nonZeroInds <- which((log2ScaleCount>0)==TRUE)
  lowerbound <-
    ceiling(stats::quantile(log2ScaleCount[nonZeroInds] ,0.25)-
              1.5*stats::IQR(log2ScaleCount[nonZeroInds]))

  upperbound <-
    ceiling(stats::quantile(log2ScaleCount[nonZeroInds] ,0.75)+
              1.5*stats::IQR(log2ScaleCount[nonZeroInds]))

  acceptedInds <- intersect(which((log2ScaleCount<upperbound)==TRUE),
                            which((log2ScaleCount>lowerbound)==TRUE) )

  # try to fit log normal distribution without using outliers
  flag <- TRUE
  try({
    distN <- MASS::fitdistr(log2ScaleCount[acceptedInds], densfun = "normal")
    pvals <- stats::pnorm(log2ScaleCount,
                          mean=(distN$estimate)[1],
                          sd=(distN$estimate)[2],
                          lower.tail = FALSE
    )
    flag <- FALSE
  })

  if(flag)
  {
    warning("Low number of sites and genes")

    # Try to fit the distribution with considering outliers due to low samples
    try({
      distN <- MASS::fitdistr(log2ScaleCount, densfun = "normal")
      pvals <- stats::pnorm(log2ScaleCount,
                            mean=(distN$estimate)[1],
                            sd=(distN$estimate)[2],
                            lower.tail = FALSE
      )
      flag <- FALSE
    })

  }

  if(flag)
  {
    stop("Cannot fit the log-normal distirbution. Use more sites and genes")
  }


  return(pvals)

  }




