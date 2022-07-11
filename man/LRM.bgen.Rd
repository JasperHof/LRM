\name{LRM.bgen}
\alias{LRM.bgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear regression on martingale residuals for .bgen files.
}
\description{
A fast and accurate method for a recurrent event GWAS
}
\usage{
LRM.bgen(
  bgenfile,
  gIDs,
  obj.null,
  output.file,
  chr = NULL,
  missing.cutoff = 0.05,
  min.maf = 0.05,
  p.cutoff = 0.001,
  memory = 512,
  maxchunksize = 5e4
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgenfile}{
    character, the name of the .bgen file, without the '.bgen' extension
}
  \item{gIDs}{
    character vector of genotype IDs, in the same order as the .bgen file
}
  \item{obj.null}{
    an R object returned from 'Null_model()'
}
  \item{output.file}{
    character, the name of the output file for GWAS summary statistics
  }
  \item{chr}{
  chromosome number
  }
  \item{missing.cutoff}{
  cut-off rate for number of missing genotypes
  }
  \item{min.maf}{
  minimum minor allele frequency
  }
  \item{p.cutoff}{
  cut-off P value for implementing a saddle-point approximation
  }
  \item{memory}{
  memory of the device in MB, used to compute chunk size
  }
  \item{maxchunksize}{
  maximum size of chunks of SNPs to be simultaneously analyzed in LRM
  }
}
\details{
To run LRM for .bgen files, the user should first fit a null model based on the phenotype data. This null model computes martingale residuals, and subsequently the distribution of martingale residuals is empirically estimated using a the saddle-point approximation.

After computing the null model, LRM.bgen() computes association statistics between SNPs in the .bed file and martingale residuals.
}

\value{
The output file contains a matrix with the following columns:
\item{rsid}{
  The rs-id of the SNP
}
\item{SNP}{
  Identifier of SNP in 'CHR_BP_A1_A2' format
}
\item{Chr}{
  Chromosome
}
\item{BP}{
  Position on the chromosome, in base pair
}
\item{A1}{
  A1 allele
}
\item{A2}{
  A2 allele
}
\item{Beta}{
  Slope coefficient of the linear regression
}
\item{pMG}{
  P value from linear regression on martingale residuals without SPA correction
}
\item{pSPA}{
  P value from linear regression, corrected using a saddle-point approximation (recommended)
}
\item{SE}{
  Standard error of beta in linear regression
}
\item{Z}{
  Z-statistic in linear regression
}
\item{MAF}{
  Minor allele frequency of SNP
}
\item{Missing}{
  Fraction of missing SNPs
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jasper Hof
}
\note{
Application is similar to the LRM.bed() function for .bed files, however the sample IDs of the .bgen file should be exactly the same as the 'gIDs' argument.
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

