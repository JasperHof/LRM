\name{LRM}
\alias{LRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Linear regression on martingale residuals
}
\description{
Carry out the linear regression on martingale resduals based on an input null model and genotype data.
}
\usage{
LRM(obj.null, Geno.mtx, missing.cutoff = 0.05, min.maf = 0.05, p.cutoff = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.null}{
    Null model
}
  \item{Geno.mtx}{
    Genotype matrix, with subjects in row names and SNP names in column names.
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
#
# Bla bla
#
#
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
