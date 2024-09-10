library(tidybayes)

args_null=function (par, argl, default) 
{
  if (is.null(argl[[par]])) 
    return(default)
  return(argl[[par]])
}


summary.pibblefit=function (object, pars = NULL, use_names = TRUE, as_factor = FALSE, 
                            gather_prob = FALSE, ...) 
{
  if (is.null(pars)) {
    pars <- c()
    if (!is.null(object$Eta)) 
      pars <- c(pars, "Eta")
    if (!is.null(object$Lambda)) 
      pars <- c(pars, "Lambda")
    if (!is.null(object$Sigma)) 
      pars <- c(pars, "Sigma")
    pars <- pars[pars %in% names(object)]
  }
  if (summary_check_precomputed(object, pars)) 
    return(object$summary[pars])
  mtidy <- dplyr::filter(pibble_tidy_samples(object, use_names, 
                                             as_factor), .data$Parameter %in% pars)
  suppressWarnings({
    vars <- c()
    if ("Eta" %in% pars) 
      vars <- c(vars, "coord", "sample")
    if ("Lambda" %in% pars) 
      vars <- c(vars, "coord", "covariate")
    if (("Sigma" %in% pars) & (object$coord_system != "proportions")) {
      vars <- c(vars, "coord", "coord2")
    }
    vars <- unique(vars)
    vars <- rlang::syms(vars)
    mtidy <- dplyr::group_by(mtidy, .data$Parameter, !!!vars)
    if (!gather_prob) {
      mtidy <- mtidy %>% summarise_posterior(.data$val, 
                                             ...) %>% dplyr::ungroup() %>% split(.$Parameter) %>% 
        purrr::map(~dplyr::select_if(.x, ~!all(is.na(.x))))
    }
    else if (gather_prob) {
      mtidy <- mtidy %>% dplyr::select(-.data$iter) %>% 
        tidybayes::mean_qi(.data$val, .width = c(0.5, 
                                                 0.8, 0.95, 0.99)) %>% dplyr::ungroup() %>% 
        split(.$Parameter) %>% purrr::map(~dplyr::select_if(.x, 
                                                            ~!all(is.na(.x))))
    }
  })
  return(mtidy)
}

summary_check_precomputed=function (m, pars) 
{
  if (!is.null(m$summary)) {
    if (all(!is.null(m$summary[pars]))) 
      return(TRUE)
  }
  return(FALSE)
}
