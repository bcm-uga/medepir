% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simu.R
\name{create_exp_grp}
\alias{create_exp_grp}
\title{Function to create the experimental group}
\usage{
create_exp_grp(n, plates = NA, sex = F, age = F)
}
\arguments{
\item{n}{The number of patients}

\item{plates}{NA by default, if you want a plate effect put a vector of plates names.}

\item{sex}{False by default, True if you want a sex parameter.}

\item{age}{False by defaul, True if you want an age parameter.}
}
\value{
This function return a dataframe with different random covariables.
}
\description{
This function make a random experimental group.
}
