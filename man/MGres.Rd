\name{MGres}
\alias{MGres}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Compute martingale residuals from a coxme null model
}
\description{
   The cumulative baseline hazards are computed from the null model, based on Breslows formula. Subsequently, martingale residuals are computed.
}
\usage{
MGres(
  fitme,
  data
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitme}{
    A 'coxme' object that is fitted for the null model
}
  \item{data}{
    The data frame that is used for fitting the coxme object. This data frame should include the samples under column name 'subjects'.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    A numeric vector containing outcome martingale residuals for all subjects
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
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
