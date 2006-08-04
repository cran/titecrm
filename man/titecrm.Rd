\name{titecrm}
\alias{titecrm}
\alias{print.mtd}
\alias{plot.mtd}
\alias{crmh}
\alias{crmht}
\alias{vcrm}
\alias{myjitter}
\title{TITE-CRM}
\description{
  Returns an object of class \code{mtd} that summarizes the dose
  assignments and recommends a dose for the next patient
  in a phase I trial using TITE-CRM.
}

\usage{
titecrm(prior, target, tox, level, n=length(level), weights=NULL, 
followup=NULL, obswin=NULL, scheme="linear", dosename=NULL, include=1:n, 
pid=1:n, method="bayes", scale=sqrt(1.34), model.detail=TRUE, 
patient.detail=TRUE)
}

\arguments{
  \item{prior}{A vector of initial estimates of toxicity probabilities
    associated the doses.}
  \item{target}{The target DLT rate.}
  \item{tox}{A vector of patient outcomes; 1 indicates a toxicity, 0
    otherwise.}
  \item{level}{A vector of dose levels assigned to patients.  The length
    of \code{level} must be equal to that of \code{tox}.}
  \item{weights}{A vector of weights assigned to observations.  A
    weight must be between 0 and 1.  If given, the arguments
    \code{followup}, \code{obswin}, and \code{scheme} will be ignored.
    If not supplied, users must provide \code{followup} and
    \code{obswin}.  The length of \code{weights} must be equal to that
    of \code{tox}.}
  \item{n}{The number of enrollments.}
  \item{followup}{A vector of follow-up times of patients.  If not
    supplied, users must provide \code{weights}.}
  \item{obswin}{The observation window with respect to which the MTD is
    defined.  If not supplied, users must provide \code{weights}.}
  \item{scheme}{A character string to specify the method for assigning
    weights.  Default is ``linear''.  Adaptive weight using Kaplan-Meier
    ``KM'' is to be made available.}
  \item{dosename}{A vector containing the names of the regimens/doses
    used.  The length of \code{dosename} must be equal to that of
    \code{prior}.}
  \item{include}{A subset of patients included in the dose calculation.}
  \item{pid}{Patient ID provided in the study.  Its length must be equal
    to that of \code{level}.}
  \item{method}{A character string to specify the method for parameter
    estimation.  The default method ``bayes'' estimates the model
    parameter by the posterior mean.  Estimation using ``mle'' is to be
    made available.}  
  \item{scale}{Standard deviation of the normal prior of the model
    parameter.  Default is sqrt(1.34).}
  \item{model.detail}{If TRUE, the model content of an ``mtd'' object will
    be displayed in detail.}
  \item{patient.detail}{If TRUE, patient summary will be given in
    detail.}
}
\value{
  An object of class ``mtd'' is returned, consisting of the summary
  of dose assignments thus far and the recommendation of dose for the
  next patient.

  \item{prior}{Initial estimates of toxicity probabilities.}
  \item{ptox}{Updated estimates of toxicity probabilities.}
  \item{target}{The target probability of toxicity at the MTD.}
  \item{recommend}{The recommended dose level for the next patient.}
  \item{scale}{The standard deviation of the normal prior.}
  \item{estimate}{Estimate of the model parameter.}
  \item{level}{Dose levels assigned to patients.}
  \item{tox}{Patients' toxicity indications.}
  \item{followup}{Follow-up times of patients.}
  \item{obswin}{Observation window with respect to which the MTD is defined.}
  \item{weights}{Weights assigned to patients.}
  \item{scheme}{Weighting scheme.}
}
\details{
  Dose-toxicity relationship is assumed as an empiric power model
  $a_i^{\exp(\beta)}$ where $a_i$ is the initial estimate of toxicity
  probability of dose level i and the model parameter $\beta$ has a normal
  prior with mean 0 and scale to be provided by users.
}

\references{
  Cheung, Y. K. and Chappell, R. (2000). Sequential designs for phase I
  clinical trials with late-onset toxicities.  Biometrics 56:1177-1182.
}

\examples{
# Create a simple data set
prior <- c(0.05,0.10,0.20,0.35,0.50,0.70)
target <- 0.2
level <- c(3,4,4,3,3,4,3,2,2,2)
y <- c(0,0,1,0,0,1,1,0,0,0)
u <- c(1,1,0.8,1,1,0.6,0.45,0.25,1/6,1/12)
tau <- 1
foo <- titecrm(prior,target,y,level,followup=u,obswin=1)
rec <- foo$recommend  # recommend a dose level for next patient
}

\keyword{datasets}
