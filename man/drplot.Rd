\name{drplot}
\alias{drplot}
\title{Dose-Response Curve Plotting}
\description{
  Function for plotting a dose-response curve of an ``mtd'' object or an
  ``sim'' object.
}

\usage{
drplot(x,n=length(x$level),include=n) 
}

\arguments{
  \item{x}{The object to be plotted.  Either an ``mtd'' or ``sim'' object.}
  \item{n}{Sample size of the trial in the object.}
  \item{include}{A vector of patient numbers indicating which curve(s) to be
    plotted.  Its length cannot be greater than \code{n}.  The default
    is set as ``n''.  That is, plot the curve based on all patients in
    the trial.}
}

\details{
  Can also use the generic function \code{print} for plotting an ``mtd''
  or ``sim'' object.
}

\examples{
PI <- c(0.10,0.20,0.40,0.50,0.60,0.65)
prior <- c(0.05,0.10,0.20,0.35,0.50,0.70)
target <- 0.2
# Generate a trial of size 24
foo <- titesim1(PI,prior,target,24,obswin=6,rate=4,accrual="poisson")
drplot(foo)  # Plot the dose-toxicity curve based on all 24 patients
}

\keyword{datasets}
