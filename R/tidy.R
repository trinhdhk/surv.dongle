#' Tidy competing risk mstate survival fit.
#' @description
#' This complements gg_ajurvplot2, provides the intermediate data for advanced ggplot manipulation
#' @param formula A standard model formula, with survival on the left and covariates on the right
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model
#' @param weights,subset,na.action,,count,id,timefix parameter passed to finegray model. See \code{\link[survival:finegray]{finegray}},
#' @param main.event,competing.event main and competing event to plot
#' @param tidy_types [c(main='risk', competing='survival')] summarisiation for main and competing events
#' @param ... other parameters passed to finegray (excluding etype)
#' @return a tibble
#' @export
tidy_ajcomprsk <- function(formula, data, weights, subset, na.action, main.event, competing.event, count, id, timefix, tidy_types, ...){
  dot <- list(...)
  fgargs <- match.call()[-1]
  fgargs <- as.list(fgargs[setdiff(names(fgargs), names(dot))])
  fgargs$tidy_types <- NULL
  etypes <- with(fgargs, list(main.event, competing.event))
  fgargs.each <- lapply(etypes, function(etype) {
    new.arg <- fgargs
    new.arg$etype <- etype
    new.arg$main.event <- new.arg$competing.event <- NULL
    new.arg
  })
  if(is.null(names(tidy_types))) names(tidy_types) <- c('main', 'competing')
  fg1 <- do.call(survival::finegray, fgargs.each[[1]])
  fg2 <- do.call(survival::finegray, fgargs.each[[2]])
  # browser()
  fml <- force(update(formula, survival::Surv(fgstart,fgstop,fgstatus)~.))

  sf1 <- eval(substitute(survival::survfit(fml, data = fg1,  weights=fgwt), list(fml=fml)))
  sf2 <- eval(substitute(survival::survfit(fml, data = fg2,  weights=fgwt), list(fml=fml)))

  tidy1 <- tidy_survfit2(sf1, type=tidy_types[['main']]) |> mutate(Event = main.event)
  tidy2 <- tidy_survfit2(sf2, type=tidy_types[['competing']]) |> mutate(Event = competing.event)

  rbind(tidy1, tidy2) |> tibble::as_tibble()
}


#' Tidy survfit with strata save
#' @description
#' Same as \link[ggsurvfit:tidy_survfit]{tidy_survfit} but with strata saved.
#' @param x,times,type same as \link[ggsurvfit:tidy_survfit]{tidy_survfit}
#' @export
tidy_survfit2 <- function(x, times=NULL, type='survival'){
  ggsurvfit::tidy_survfit(x, times=times, type) |>
    untidy_strata()
}

#' Untidy the tidy_survfit
untidy_strata <- function(x, strata.col='strata', strata.lab.col='strata_label') {
  if (!'strata' %in% names(x)) return(x)
  # tidy of survfit has strata_lab to be strata
  if (all(x[[strata.lab.col]]=='strata'))
    return(._untidy_survfit_(x, strata.col='strata'))
  ._untidy_survfit2_(x, strata.col='strata', strata.lab.col='strata_label')

}

._untidy_survfit_ <-
  function(x, strata.col='strata'){
    x2 <- x |>
      dplyr::mutate(strata.split = strsplit(as.character(x[[strata.col]]), ', ')) |>
      tidyr::unnest(cols=strata.split) |>
      dplyr::mutate(
        strata.split = unlist(strata.split),
        strata.names = strsplit(strata.split, '=') %>% purrr::transpose() %>% `[[`(1) |> unlist(),
        strata.values = strsplit(strata.split, '=') %>% purrr::transpose() %>% `[[`(2)  |> unlist(),
        strata.split=NULL
      ) |>
      tidyr::pivot_wider(names_from = strata.names, values_from = strata.values)

    x2
  }

._untidy_survfit2_ <-
  function(x, strata.col='strata', strata.lab.col='strata_label'){
    x2 <- x |>
      dplyr::mutate(
        strata.split = strsplit(as.character(x[[strata.col]]), ', '),
        strata.lab.split = strsplit(as.character(x[[strata.lab.col]]), ', ')
        ) |>
      tidyr::unnest(cols=c(strata.split, strata.lab.split)) |>
      tidyr::pivot_wider(names_from = strata.lab.split, values_from = strata.split)
  }
