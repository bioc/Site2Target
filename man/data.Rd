% Generated manually by Peyman Zarrineh
\name{data}
\alias{data}
\title{MEIS cardiomyocytes datasets used in the package}
\format{
Gene expression WT vs MEIS KO in chr21.
MEIS binding sites in chr21.
TADs, and HiC interactions in chr21.
\describe{
\item{gene_expression.tsv}{Gene expression}

\item{MEIS_binding.tsv}{MEIS binding sites}

\item{TADs.tsv}{TADs}

\item{HiC_intensities.tsv}{HiC interactions}
}
}
\description{
Human cardiomyocytes datasets are reduced in size by only using chr21.
log fold changes of Gene expression WT vs MEIS KO from RNA-seq experiments, 
and binding sites of MEIS derived from a ChIP-seq experiment are the main
experimental datasets representing relevant gene and peak information.
HiC interactions and topologically associating domains (TADs) are derived
from a HiC experiments are auxiliary datasets related to DNA-DNA interactions.
}
\value{
Just description of data
}
\examples{

## Gene expression table

# Read gene coordinates 
geneFile=system.file("extdata", "gene_expression.tsv",
                     package="Site2Target")
geneCoords <- Table2Granges(geneFile)

# Read gene table
geneTable <- read.table(geneFile, header=TRUE)



## TF binding table

# Read peak coordinates
tfFile =system.file("extdata", "MEIS_binding.tsv",
                    package="Site2Target")
TFCoords <- Table2Granges(tfFile)

# Read MEIS binding intensities
tfTable <- read.table(tfFile, header=TRUE)


## DNA-DNA interactions

# Read TAD regions
TADsFile =system.file("extdata", "TADs.tsv",
                               package="Site2Target")
TADs <- Table2Granges(TADsFile)


# Read HiC interactions
HiCFile =system.file("extdata", "HiC_intensities.tsv",
                               package="Site2Target")
HiCstr1 <- Table2Granges(HiCFile, chrColName="Strand1_chr",
                      startColName="Strand1_start", endColName="Strand1_end")
HiCstr2 <- Table2Granges(HiCFile, chrColName="Strand2_chr",
                     startColName="Strand2_start", endColName="Strand2_end")

HiCTable <- read.table(HiCFile, header=TRUE)


}
