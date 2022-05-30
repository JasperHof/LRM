\name{LRM.bed}
\alias{LRM.bed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear regression on martingale residuals for .bed files.
}
\description{
Carry out a linear regression on martingale residuals for .bed files
}
\usage{
LRM.bed(
  bedfile,
  gIDs,
  obj.null,
  output.file,
  chr = NULL,
  missing.cutoff = 0.05,
  min.maf = 0.05,
  p.cutoff = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bedfile}{
    character, the name of the .bed file, without '.bed'.
}
  \item{gIDs}{
    character vector of genotype IDs, in the same order as the .bgen file.
}
  \item{obj.null}{
    an R object returned from 'Null_model()'.
}
  \item{output.file}{
    character, the name of the output file for GWAS summary statistics.
  }
  \item{chr}{
  Chromosome number.
  }
  \item{missing.cutoff}{
  Cut-off rate for number of missing genotypes.
  }
  \item{min.maf}{
  Minimum minor allele frequency.
  }
  \item{p.cutoff}{
  Cut-off P value for implementing a saddle-point approximation.
  }
}
\details{
To run LRM for .bed files, the user should first fit a null model based on the phenotype data. This null model is used to compute martingale residuals, and their distribution is empirically estimated using a the saddle-point approximation.

Subsequently, this function can be used to compute association statistics between SNPs in the .bed file and martingale residuals.
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Load required packages
library(coxme)
library(seqminer)

# Load genotype and phenotype data
famfile = system.file('extdata', 'bedfile.fam', package = 'LRM')
fam = read.table(famfile)

bedfile = gsub('.fam', '', famfile)

datafile = system.file('extdata', 'Recurrencedata', package = 'LRM')
data <- read.csv2(datafile, stringsAsFactors = FALSE)
data = as.data.frame(data)

# Fit a non-stratified calendar model as null model
fitme = coxme(Surv(tstart, tstop, Status) ~ X1 + X2 + X3 + (1|subject), data = data, ties = 'breslow')

# Compute the martingale residuals and saddle-point approximation
obj.null = Null_model(fitme, data, IDs = unique(data$subject))

# Perform association tests between .bed file and martingale residuals. Output is written to 'bedoutput'
LRM.bed(bedfile, fam$v2 , obj.null, 'bedoutput')
}
