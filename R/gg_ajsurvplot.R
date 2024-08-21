#' Plot the Aalen-Johansen curve with ggplot for competing risks
#' @description
#' Extend \link[survminer:ggsurvplot]{ggsurvplot} functionality to plot competing risk curve
#' @param formula A standard model formula, with survival on the left and covariates on the right
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model
#' @param event event of interest, passed to finegray as etype.
#' @param weights,subset,na.action, count,id,timefix parameter passed to finegray model. See \code{\link[survival:finegray]{finegray}},
#' @param backend [\code{ggsurvfit}] either 'survminer' or 'ggsurvfit'
#' @param ... optional parameters passed to \code{ggsurvplot} or ].
#' @export
gg_ajsurvplot <- function(formula, data, weights, subset, na.action, etype, count, id, timefix, backend = c('ggsurvfit', 'survminer'), ...){
  dot <- list(...)
  fgargs <- match.call()[-1]
  fgargs <- as.list(fgargs[setdiff(names(fgargs), names(dot))])
  fgargs$etype <- event
  fgargs$event <- NULL
  fg <- do.call(survival::finegray, fgargs)
  backend <- match.arg(backend)

  fml <- force(update(formula, Surv(fgstart,fgstop,fgstatus)~.))
  sf <- eval(substitute(survival::survfit(fml, data = fg,  weights=fgwt), list(fml=fml)))
  if (backend == 'survminer')
    return(survminer::ggsurvplot(
      sf,
      data=fg,
      ...
  ))

  ggsurvfit::ggsurvfit(
    sf,
    data=fg,
    ...
  )
}

#' ggsurvfit with strata saved
#' @description
#' Same as \link[ggsurvfit:ggsurvfit]{ggsurvfit} but with strata saved, so you can facet by covariables.
#'
#' @param x,type,linetype_aes,theme,... parameters passed directly \link[ggsurvfit:ggsurvfit]{ggsurvfit}.
#' @export
ggsurvfit2 <- rlang::new_function(
  args = rlang::fn_fmls(ggsurvfit::ggsurvfit),
  body = {
    # require(cli)
    require(ggsurvfit)
    tt <- ggsurvfit::ggsurvfit |>
    rlang::fn_body() |>
    as.character() |>
    gsub('tidy_survfit', 'surv.dongle::tidy_survfit2', x=_)|>
    c('}') |> paste(collapse = '\n') |>
    # parse(text=_) |>
    rlang::parse_expr()
    environment(tt) <- asNamespace('ggsurvfit')
    tt
    }
)
environment(ggsurvfit2) <- asNamespace('ggsurvfit')


#' Plot 2 Aalen-Johansen curves with ggplot for competing risks
#' @description
#' Extend \link[ggsurvfit:ggsurvfit]{ggsurvfit} functionality to plot competing risk curve, with main risk and competing risk in one plot
#'
#' @param formula A standard model formula, with survival on the left and covariates on the right
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model
#' @param weights,subset,na.action,,count,id,timefix parameter passed to finegray model. See \code{\link[survival:finegray]{finegray}},
#' @param main.event,competing.event main and competing event to plot
#' @param facet.by formula. default to all strata
#' @param ci [TRUE] Plot confidence band?
#' @param theme [FALSE] plot in monochrome? either FALSE (default) or a string of color, TRUE is equivalent to "black"
#' @param style ['alternate'] plot in alternative or overlaid style?
#' @param ... other parameters passed to finegray (excluding etype)
#' @import ggplot2
#' @export
gg_ajsurvplot2 <- function(formula, data, weights, subset, na.action, main.event, competing.event, facet.by = ~strata, count, id, timefix, ci=TRUE, monochrome = FALSE, style=c('alternate', 'overlaid'), ...){

  dot <- list(...)
  gargs <- match.call()[-1]
  gargs <- as.list(gargs[setdiff(names(gargs), names(dot))])
  gargs$facet.by <- gargs$ci <- gargs$monochrome <- NULL
  style <- match.arg(style)
  gargs$tidy_type <- c(main='risk', competing='survival')
  if (style == 'overlaid') gargs$tidy_type[['competing']] <- 'risk'

  dt <- do.call(tidy_competingevent, gargs)

  # if any variable in data is a factor, dt should be so
  for (v in formula.tools::rhs.vars(formula)){
    if (is.factor(data[[v]])) dt[[v]] <- factor(dt[[v]], levels=levels(data[[v]]))
  }

  strata_levels_order <-
    lapply(formula.tools::rhs.vars(formula)[!formula.tools::rhs.vars(formula) %in%
                                              formula.tools::get.vars(facet.by)],
           function(v){
             if (is.factor(data[[v]])) return(as.numeric(factor(dt[[v]], levels=levels(data[[v]]))))
             as.numeric(factor(dt[[v]]))
           }) |> do.call(order, args=_)

  facet <- facet.vars <- NULL
  if (!is.null(facet.by)){
    facet <-
      if (length(formula.tools::lhs.vars(facet.by))) facet_grid(facet.by) else facet_wrap(facet.by)
    facet.vars <- formula.tools::get.vars(facet.by)
  }

  for (v in facet.vars) dt$strata <-
    sapply(strsplit(as.character(dt$strata), ', '),
           function(x) x[!grepl(paste0(v,'='), x, fixed=TRUE)]|>paste(collapse=', '))

  if (all(dt$strata == '')) plt <- ggplot(dt,aes(x=time, y=estimate, ymin=conf.low, ymax=conf.high, group=Event))
  else if (!isFALSE(monochrome)) {

    dt$strata <- factor(dt$strata, levels = unique(dt$strata)[strata_levels_order])
    plt <- ggplot(dt,aes(x=time, y=estimate, ymin=conf.low, ymax=conf.high, linetype=strata, color=Event))
    if (ci) plt <- plt +
      ggsurvfit::stat_stepribbon(alpha=.5, linewidth=.2, fill='transparent')
    plt <-  plt +
      geom_step(linewidth=1) +
      scale_color_manual(values=rep(if (isTRUE(monochrome)) 'black' else as.character(monochrome),2), guide=NULL)
  }

  else {

    dt$strata <- factor(dt$strata, levels = unique(dt$strata)[strata_levels_order])
    plt <- ggplot(dt,aes(x=time, y=estimate, ymin=conf.low, ymax=conf.high, linetype=Event, color=strata, fill=strata))
    if (ci) plt <- plt +  ggsurvfit::stat_stepribbon(alpha=.3, color='transparent')
    plt <- plt +  geom_step(linewidth=1) +
      scale_linetype_manual(values=rep('solid',2), guide=NULL)
  }

  plt + ggsurvfit::theme_ggsurvfit_default() + facet

}
