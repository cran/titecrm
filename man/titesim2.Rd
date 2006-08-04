\name{titesim2}
\alias{titesim2}
\title{Two-stage TITE-CRM Simulator}
\description{
  Returns an object of class ``sim'' that generates and summarizes the dose
  assignments of a simulated trial by a two-stage TITE-CRM.
}

\usage{
titesim2(PI, prior, target, n, x0, obswin=1, tgrp=obswin, rate=1, 
accrual="fixed", surv="uniform", scheme="linear", method="bayes",
scale=sqrt(1.34), seed=1099)
}

\arguments{
  \item{PI}{A vector of the true toxicity probabilites associated with
    the doses.}
  \item{prior}{A vector of initial estimates of toxicity probabilities
    associated with the doses.  Must be of same length as \code{PI}.}
  \item{target}{The target DLT rate.}
  \item{n}{Sample size of the trial.}
  \item{x0}{A vector of treatment sequence according to the initial
    design.  Must be of length \code{n}.}
  \item{obswin}{The observation window with respect to which the MTD is
    defined.}
  \item{tgrp}{The minimum waiting time between two dose cohorts at the
    initial stage.  Default is set as \code{obswin}, i.e., complete
    follow-up in all current patients is required before escalation to
    the next dose group.}
  \item{rate}{Patient arrival rate: Expected number of arrivals per
    observation window.  Example: obswin=6 and rate=3 means expecting 3
    patients arrive in 6 time units.}
  \item{accrual}{Patient arrival scheme.  Default is ``fixed'' whereby
    inter-patient arrival is fixed.  Alternatively, use ``poisson'' to
    simulate patient arrivals by the Poisson process.}
  \item{surv}{Distribution for time-to-toxicity.  Default is ``uniform''
    where toxicity, if occurs, occurs uniformly on the interval
    [0,\code{obswin}].  Other survival distributions including exponential and
    Weibull are to be made available.}
  \item{scheme}{A character string to specify the method for assigning
    weights.  Default is ``linear''.  Adaptive weight using Kaplan-Meier
    ``KM'' is to be made available.}
  \item{method}{A character string to specify the method for parameter
    estimation.  The default method ``bayes'' estimates the model
    parameter by the posterior mean.  Estimation using ``mle'' is to be
    made available.}  
  \item{scale}{Standard deviation of the normal prior of the model
    parameter.  Default is sqrt(1.34).}
  \item{seed}{Seed of the random number generator.}
}
\value{
  An object of class ``mtd'' is returned, consisting of the summary
  of dose assignments and the final dose recommendation in a simulated trial.
  
  \item{PI}{True toxicity probabilites.}
  \item{prior}{Initial estimates of toxicity probabilities.}
  \item{target}{The target probability of toxicity at the MTD.}
  \item{recommend}{The recommended dose level for the next patient.}
  \item{scale}{The standard deviation of the normal prior.}
  \item{estimate}{Estimate of the model parameter.}
  \item{level}{Dose levels assigned to patients.}
  \item{tox}{Patients' toxicity indications.}
  \item{arrival}{Patients' arrival times.}
  \item{ttox.pt}{Patients' times-to-toxicity.}
  \item{ttox.cal}{Patients' times-to-toxicity on study time.}
  \item{obswin}{Observation window with respect to which the MTD is defined.}
  \item{weights}{Weights assigned to patients.}
  \item{scheme}{Weighting scheme.}
  \item{accrual}{Patient's arrival pattern.}
  \item{rate}{Rate of patient's arrival.}
  \item{surv}{Distribution of time-to-toxicity.}
}
\details{
  Dose-toxicity relationship is assumed as an empiric power model
  $a_i^{\exp(\beta)}$ where $a_i$ is the initial estimate of toxicity
  probability of dose level i and the model parameter $\beta$ has a normal
  prior with mean 0 and scale to be provided by users.

  The simulation of a trial run by a two-stage TITE-CRM.  Users need to
  provide an initial sequence before switching the TITE-CRM.
}


\references{
  Cheung, Y. K. and Chappell, R. (2000). Sequential designs for phase I
  clinical trials with late-onset toxicities.  Biometrics 56:1177-1182.

  Cheung, Y. K. (2005). Coherence principles in dose-finding studies.
  Biometrika 92:863-873.
}

\examples{
PI <- c(0.10,0.20,0.40,0.50,0.60,0.65)
prior <- c(0.05,0.10,0.20,0.35,0.50,0.70)
target <- 0.2
x0 <- c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,9))
# Generate a trial of size 24
foo <- titesim2(PI,prior,target,24,x0, obswin=6,rate=4,accrual="poisson")
rec <- foo$recommend  # recommend a dose level for next patient
plot(foo)  # summarize trial graphically
}

\keyword{datasets}
