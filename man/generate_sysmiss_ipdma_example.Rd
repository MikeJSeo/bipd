% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.example.data.R
\name{generate_sysmiss_ipdma_example}
\alias{generate_sysmiss_ipdma_example}
\title{Generate a simulated IPD-MA data with systematically missing covariates}
\usage{
generate_sysmiss_ipdma_example(
  Nstudies = 10,
  Ncov = 5,
  sys_missing_prob = 0.1,
  magnitude = 0.3,
  heterogeneity = 0.1,
  interaction = TRUE
)
}
\arguments{
\item{Nstudies}{number of studies. Default is 10.}

\item{Ncov}{number of covariates in total. Options are 5 or 10 studies. Default is set to 5.}

\item{sys_missing_prob}{probability of systematically missing studies for each covariates. Default is set to 0.3.}

\item{magnitude}{magnitude of the regression estimates (the mean). Default is set to 0.2.}

\item{heterogeneity}{heterogeneity of regression estimates across studies. Default is set to 0.1.}

\item{interaction}{whether to include treatment indicator and treatment}
}
\value{
returns simulated IPD-MA data with systematically missing covariates
}
\description{
Generate a simulated IPD-MA data with systematically missing covariates
}
\examples{
simulated_dataset <- generate_sysmiss_ipdma_example(Nstudies = 10, Ncov = 5, sys_missing_prob = 0.3, 
magnitude = 0.2, heterogeneity = 0.1)
head(simulated_dataset)
}
