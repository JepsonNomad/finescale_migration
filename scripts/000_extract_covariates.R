
extract_covariates <- function(x, ...) {
  UseMethod("extract_covariates", x)
}

#' @export
#' @rdname extract_covariates
extract_covariates.track_xy <- function(x, covariates, ...) {
  extract_covar_base(x, covariates, ...)
}

#' @export
#' @rdname extract_covariates
extract_covariates.random_points <- function(x, covariates, ...) {
  extract_covar_base(x, covariates, ...)
}

#' @export
#' @rdname extract_covariates
extract_covariates.steps_xy <- function(x, covariates, where = "end", ...) {
  if (class(covariates) == "SpatRaster") {
    if (where == "both") {
      x_start <- terra::extract(covariates, 
                                as.matrix(x[, c("x1_", "y1_")]),
                                ...)
      names(x_start) <- paste0(names(x_start), "_start")
      x_end <- terra::extract(covariates, as.matrix(x[, c("x2_", "y2_")]),
                               ...)
      names(x_end) <- paste0(names(x_end), "_end")
      x_all <- cbind(x_start, x_end)
      x[names(x_all)] <- as.data.frame(x_all)
    } else {
      x[names(covariates)] <- if (where == "end") {
        as.data.frame(terra::extract(covariates, 
                                     as.matrix(x[, c("x2_", "y2_")]), ...))
      } else if (where == "start") {
        as.data.frame(terra::extract(covariates, 
                                     as.matrix(x[, c("x1_", "y1_")]), ...))
      }
    }
    x
  } else {
    stop("no terra")
  }
}

extract_covar_base <- function(x, covars, ...) {
  if (class(covars) == "SpatRaster") {
    x[names(covars)] <- as.data.frame(terra::extract(covars, 
                                                     x[, c("x_", "y_")], ...))
    x
  } else {
    stop("no terra")
  }
}





extract_covariates_var_time <- function(x, ...) {
  UseMethod("extract_covariates_var_time", x)
}

#' @export
#' @rdname extract_covariates
extract_covariates_var_time.track_xyt <- function(
  x, covariates, when = "any", max_time,
  name_covar = "time_var_covar", ...) {
  x[name_covar] <- extract_covar_var_time_base(
    cbind(x$x_, x$y_),
    x$t_, covariates, when, max_time, ...)
  x
}

#' @export
#' @rdname extract_covariates
extract_covariates_var_time.steps_xyt <- function(
  x, covariates, when = "any", max_time, name_covar = "time_var_covar",
  where = "end", ...) {
  
  if (where == "start") {
    x[name_covar] <- extract_covar_var_time_base(
      cbind(x$x1_, x$y1_),
      x$t1_, covariates, when, max_time, ...)
  } else if (where == "end") {
    x[name_covar] <- extract_covar_var_time_base(
      cbind(x$x2_, x$y2_),
      x$t2_, covariates, when, max_time, ...)
  } else if (where == "both") {
    x[paste0(name_covar, "_start")] <- extract_covar_var_time_base(
      cbind(x$x1_, x$y1_),
      x$t1_, covariates, when, max_time, ...)
    x[paste0(name_covar, "_end")] <- extract_covar_var_time_base(
      cbind(x$x2_, x$y2_),
      x$t2_, covariates, when, max_time, ...)
  }
  x
}


extract_covariates_var_time <- function(x, ...) {
  UseMethod("extract_covariates_var_time", x)
}

#' @export
#' @rdname extract_covariates
extract_covariates_var_time.track_xyt <- function(
  x, covariates, when = "any", max_time,
  name_covar = "time_var_covar", ...) {
  x[name_covar] <- extract_covar_var_time_base(
    cbind(x$x_, x$y_),
    x$t_, covariates, when, max_time, ...)
  x
}

#' @export
#' @rdname extract_covariates
extract_covariates_var_time.steps_xyt <- function(
  x, covariates, when = "any", max_time, name_covar = "time_var_covar",
  where = "end", ...) {
  
  if (where == "start") {
    x[name_covar] <- extract_covar_var_time_base(
      cbind(x$x1_, x$y1_),
      x$t1_, covariates, when, max_time, ...)
  } else if (where == "end") {
    x[name_covar] <- extract_covar_var_time_base(
      cbind(x$x2_, x$y2_),
      x$t2_, covariates, when, max_time, ...)
  } else if (where == "both") {
    x[paste0(name_covar, "_start")] <- extract_covar_var_time_base(
      cbind(x$x1_, x$y1_),
      x$t1_, covariates, when, max_time, ...)
    x[paste0(name_covar, "_end")] <- extract_covar_var_time_base(
      cbind(x$x2_, x$y2_),
      x$t2_, covariates, when, max_time, ...)
  }
  x
}

extract_covar_var_time_base <- function(
  xy, t, covariates, when = "any",
  max_diff, ...) {
  
  if (is.null(terra::time(covariates))) {
    stop("Covariates do not have a Z column.")
  }
  
  if (!is(max_diff, "Period")) {
    stop("`max_diff` is not of class `Period`.")
  }
  max_diff <- lubridate::period_to_seconds(max_diff)
  t_covar <- as.numeric(as.POSIXct(terra::time(covariates)))
  t_obs <- as.numeric(as.POSIXct(t))
  
  # Fun to find closest point
  which_rast <- function(t_diffs, where, max_diff) {
    wr <- if (when == "after") {
      which.min(t_diffs[t_diffs >= 0])
    } else if (when == "before") {
      which.min(abs(t_diffs[t_diffs <= 0])) + sum(t_diffs > 0)
    } else if (when == "any") {
      which.min(abs(t_diffs))
    }
    if (length(wr) == 0) {
      NA
    } else if (max_diff < abs(t_diffs[wr])) {
      NA
    } else {
      wr
    }
  }
  
  wr <- sapply(t_obs, function(x) which_rast(x - t_covar, when, max_diff))
  # ev <- raster::extract(covariates, cbind(xy), ...)
  ev <- terra::extract(covariates, cbind(xy))
  cov_val <- ev[cbind(seq_along(wr), wr)]
  return(cov_val)
}