\name{crmsens}
\alias{crmsens}
\alias{print.homesets}
\title{Model Sensitivity in the CRM}
\description{
  Returns an object of class ``homesets'' that summarizes the sensitvity
  (robustness) of the dose-toxicity model used in the CRM.}
}

\usage{
crmsens(prior, target, eps=1e-8, maxit=100, incr=0.1, LB=0, UB=Inf, 
detail=FALSE)
}

\arguments{
  \item{prior}{A vector of initial estimates of toxicity probabilities
    associated the doses.}
  \item{target}{The target DLT rate.}
  \item{detail}{If TRUE, the details of the home sets will also be
    shown.  Default is set as FALSE.}
  \item{eps}{A control parameter for absolute tolerance.  Defaults to
    `1e-8'.}
  \item{maxit}{Maximum number of iterations allowed.  Defaults to 100.}
  \item{incr}{Step size.  Defaults to 0.1.}
  \item{LB}{Lower bound of the search space.  Defaults to 0.}
  \item{UB}{Upper bound of the search space.  Defaults to `Inf'.}    
}

\value{
  An object of class ``homesets'' is returned, describing the model
  sensitivity by indifference intervals and a summary statement.

  \item{homeset}{The home sets of the model parameter $\exp(\beta)$.}
  \item{iint}{The indifference intervals of the dose-toxicity model.}
  \item{prior}{Initial estimates of toxicity probabilities.}
  \item{target}{The target probability of toxicity at the MTD.}
}
\details{
  Dose-toxicity relationship is assumed as an empiric power model
  $a_i^{\exp(\beta)}$ where $a_i$ is the initial estimate of toxicity.
  Different choice of the initial probabilities results in various level
  of robustness of the CRM, which can be summarized by the indifference
  intervals (\code{iint}).  See Cheung and Chappell (2002) in the
  reference section.
}

\references{
  Cheung, Y. K. and Chappell, R. (2002). A simple technique to evaluate
  model sensitivity in the continual reassessment method.  Biometrics
  58:671-674.
}

\examples{
# Create a simple data set
prior <- c(0.05,0.10,0.20,0.35,0.50,0.70)
target <- 0.2
foo <- crmsens(prior,target)
print(foo)
}

\keyword{datasets}

