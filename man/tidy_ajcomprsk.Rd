% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy.R
\name{tidy_ajcomprsk}
\alias{tidy_ajcomprsk}
\title{Tidy competing risk mstate survival fit.}
\usage{
tidy_ajcomprsk(
  formula,
  data,
  weights,
  subset,
  na.action,
  main.event,
  competing.event,
  count,
  id,
  timefix,
  tidy_types,
  ...
)
}
\arguments{
\item{formula}{A standard model formula, with survival on the left and covariates on the right}

\item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model}

\item{weights, subset, na.action, , count, id, timefix}{parameter passed to finegray model. See \code{\link[survival:finegray]{finegray}},}

\item{main.event, competing.event}{main and competing event to plot}

\item{tidy_types}{[c(main='risk', competing='survival')] summarisiation for main and competing events}

\item{...}{other parameters passed to finegray (excluding etype)}
}
\value{
a tibble
}
\description{
This complements gg_ajurvplot2, provides the intermediate data for advanced ggplot manipulation
}
