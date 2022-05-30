\name{LRM.bgen}
\alias{LRM.bgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear regression on martingale residuals for .bgen files.
}
\description{
Carry out a linear regression on martingale residuals for .bgen files
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
  p.cutoff = 0.001
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgenfile}{
    character, the name of the .bgen file, without '.bgen'.
}
  \item{gIDs}{
    character vector of genotype IDs, in the same order as the .bgen file.
}
  \item{obj.null}{
    an R object returned from 'Null_model()'
}
  \item{output.file}{
    character, the name of the output file for GWAS summary statistics
  }
  \item{chr}{
  Chromosome number
  }
  \item{missing.cutoff}{
  Cut-off rate for number of missing genotypes
  }
  \item{min.maf}{
  Minimum minor allele frequency
  }
  \item{p.cutoff}{
  Cut-off P value for implementing a saddle-point approximation.
  }
}
\details{
To run LRM for .bgen files, the user should first fit a null model based on the phenotype data. This null model is used to compute martingale residuals, and their distribution is empirically estimated using a the saddle-point approximation.

Subsequently, this function can be used to compute association statistics between SNPs in the .gen file and martingale residuals.
}

\value{
The output file contains a matrix with the following columns
\item{SNP}{
  Name of SNP
}
\item{Beta}{
  Slope coefficient of the linear regression
}
\item{pMG}{
  P value from linear regression on martingale residuals
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

