---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# Site2Target

The goal of Site2Target is to to associate sets of sites/peaks to target genes. It provides peakwise-associations to associate target genes for a given set of peaks. It also provides genewise-associations which start from genes (ex. differential expressed genes) and associate peaks/sites to these genes.

## Installation

To install this package, start R (version "4.4") and enter:

``` 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Site2Target")

```

## Example

Here is an example of a peak-wise and a gene-wise association of differential genes WT vs KO of a transcription factor and binding sites of this transcription factor:

```
library(Site2Target)

## peak-wise association example 
# Read gene expression coordinates 
geneFile=system.file("extdata", "gene_expression.tsv",
                     package="Site2Target")
geneCoords <- Table2Granges(geneFile)

# Read gene expression table
geneTable <- read.table(geneFile, header=TRUE)

# Read peak coordinates of MEIS binding sites
tfFile =system.file("extdata", "MEIS_binding.tsv",
                    package="Site2Target")
TFCoords <- Table2Granges(tfFile)
tfTable <- read.table(tfFile, header=TRUE)


# Predict targets of MEIS using peakwise-association
pvals <- getTargetGenesPvals( geneCoordinates=geneCoords,
                              sites=TFCoords, distance = 50000)

topTargetNum <- 5
topTargetIndex <- order(pvals)[1:topTargetNum]

# Make a data frame of peak targets pvalues and expression logFCs

dfTopTarget <- 
  data.frame(name=geneTable$name[topTargetIndex],
             pvalue=pvals[topTargetIndex],
             exprLogC=geneTable$logFC[topTargetIndex]
             )
dfTopTarget



## gene-wise association example

# Take differential genes iformation
geneDEIndices <- which((abs(geneTable$logFC)>1)==TRUE)
indicesLen <- length(geneDEIndices)
if(indicesLen >0)
{
    geneTable <- geneTable[geneDEIndices,]
    geneCoords <- geneCoords[geneDEIndices]
}
geneDENames <- geneTable$name
geneDElogFC <- geneTable$logFC
geneCoordsDE <- geneCoords

# Associate peaks located up to 50k bp to differential genes
stats <-
genewiseAssociation(associationBy="distance",
                    geneCoordinates=geneCoordsDE,
                    geneNames=geneDENames,
                    peakCoordinates=TFCoords,
                    distance=50000,
                    outFile="Gene_TF_50K")
stats


```
# Site2Target
