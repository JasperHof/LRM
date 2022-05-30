\name{LRM.bed}
\alias{LRM.bed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear regression on martingale residuals for .bed files.
}
\description{
Carry out a linear regression on martingale residuals for .bgen files
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
  Name of the .bed file, excluding '.bed' file name extension
  }
  \item{gIDs}{
  Names of the genotype IDs in the .bed file
  }
  \item{obj.null}{
    an R object returned from 'Null_model()'
  }
  \item{output.file}{
  Name of the output file for GWAS summary statistics
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
To run LRM for .bed files, the user should first fit a null model based on the phenotype data. This null model is used to compute martingale residuals, and their distribution in the saddle-point approximation.

Subsequently, this function can be used to compute association statistics between SNPs in the .bed file and martingale residuals.
}
\value{
The output contains name of SNP, P values and other summary statistics of the linear regression
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
# Load genotype and phenotype data
bedfile = system.file('extdata', 'bedfile', package = 'LRM')
datafile = system.file('Recurrencedata', 'bedfile', package = 'LRM')
data <- read.csv2(datafile, stringsAsFactors = F)
data = as.data.frame(data)

# Fit a non-stratified calendar model as null model
fitme = coxme(Surv(tstart, tstop, Status) ~ X1 + X2 + X3 + (1|subject), data = data, ties = 'breslow')

# Compute the martingale residuals and saddle-point approximation
obj.null = Null_model(fitme, data, IDs = unique(data$subject))

# Perform association tests between .bed file and martingale residuals. Output is written to 'bedoutput'
LRM.bed('bedfile', fam$id , obj.null, 'bedoutput')
}
