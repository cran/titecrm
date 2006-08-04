\name{cohere}
\alias{cohere}
\title{Coherence of Two-stage CRM}
\description{
  Returns a message on the coherence status of a two-stage CRM
  design.}

\usage{
cohere(prior, target, n, x0, method="bayes", scale=sqrt(1.34), detail=TRUE)
}

\arguments{
  \item{prior}{A vector of initial estimates of toxicity probabilities
    associated the doses.}
  \item{target}{The target DLT rate.}
  \item{n}{Sample size of the trial.}
  \item{x0}{A vector of treatment sequence according to the initial
    design.  Must be of length \code{n}.}
  \item{method}{A character string to specify the method for parameter
    estimation.  The default method ``bayes'' estimates the model
    parameter by the posterior mean.  Estimation using ``mle'' is to be
    made available.}  
  \item{scale}{Standard deviation of the normal prior of the model
    parameter.  Default is sqrt(1.34).}
  \item{detail}{If TRUE, details about incoherent escalations
    will be displayed.}
}
\value{
  A string character giving a message on the coherence status of a
  two-stage CRM design.
}

\details{
  Dose-toxicity relationship is assumed as an empiric power model
  $a_i^{\exp(\beta)}$ where $a_i$ is the initial estimate of toxicity
  probability of dose level i and the model parameter
  $\beta$ has a normal
  prior with mean 0 and scale to be provided by users.

  A two-stage CRM is defined by two components.  First, an initial
  design sequence (\code{x0}) for dose assignments before any DLT is seen.
  When a DLT is seen, dose assignment decisions are switched to the
  CRM.  See Cheung (2005) in the reference section.
}

\references{
  Cheung, Y. K. and Chappell, R. (2000). Sequential designs for phase I
  clinical trials with late-onset toxicities.  Biometrics 56:1177-1182.

  Cheung, Y. K. (2005). Coherence principles in dose-finding studies.
  Biometrika 92:863-873.
}

\examples{
# Create a simple data set
prior <- c(0.05,0.10,0.20,0.35,0.50,0.70)
target <- 0.2
n <- 24
x0 <- c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,9))
# The above design is coherent
foo <- cohere(prior,target,n,x0)
print(foo)

# The design is incoherent if a larger target DLT rate is used.
target2 <- 0.3
foo2 <- cohere(prior,target2,n,x0)
}

\keyword{datasets}
