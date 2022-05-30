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
    character, the name of the .bgen file, without '.bgen'
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
%%  ~~ If necessary, more details than the description above ~~
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

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
