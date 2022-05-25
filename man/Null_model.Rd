\name{Null_model}
\alias{Null_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute martingale residuals and implement saddle-point approximation
}
\description{
  Computes martingale residuals based on a specified coxme null model, and computes a saddle-point approximation for the distribution of regression coefficients in LRM
}
\usage{
Null_model(fitme, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitme}{
    A 'coxme' object that is fitted for the null model
}
  \item{data}{
    The data frame that is used for fitting the coxme object. This data frame should include the samples under column name 'subjects'. }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    Object of class 'NULL_model' which includes martingale residuals, the derivatives of the CGF and the IDs.
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
