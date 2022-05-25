\name{LRM.bgen}
\alias{LRM.bgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Linear regression on martingale residuals for .bgen files
}
\description{
 Carry out a linear regression on martingale residuals for .bgen files
}
\usage{
LRM.bgen(bgenfile, gIDs, obj.null, chr = NULL, missing.cutoff = 0.05, min.maf = 0.05, p.cutoff = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgenfile}{
    The name of the .bgen file, without '.bgen'
}
  \item{gIDs}{
    Name of the genotype IDs, in the same order as the .bgen file
}
  \item{obj.null}{
    Null model object
}
  \item{chr}{
    Number of the chromosome
}
  \item{missing.cutoff}{
    Threshold for rate of missing SNPs
}
  \item{min.maf}{
    MAF threshold
}
  \item{p.cutoff}{
    p value threshold for applying saddle-point approximation
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output contains name of SNP, P values and other outcomes of the linear regression
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Jasper Hof
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
