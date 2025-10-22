#' Plot Figure 1 with effect of CHC on mortality.
#'
#' @param tau_sdid   Numeric vector. The vector contains treatment effect.
#' @param tau_b_sdid Numeric array. The bootstrapped tau array utilized for CI.
#' @param save_pdf   Bool. If TRUE, saves pdf beside generating a plot.
#' @param path       Character.
#' @export
plot_chc_effect <- function(
    tau_sdid, tau_b_sdid,
    save_pdf = FALSE, path = ""
) {
    orig_values <- c(
      -3.295, -4.272, -5.222, -6.346, -8.771,
      -9.265, -10.226, -10.071, -11.533, -14.266
    )
    tau_upper <- apply(tau_b_sdid, 1, FUN = function(x) quantile(x, probs = 0.975))
    tau_lower <- apply(tau_b_sdid, 1, FUN = function(x) quantile(x, probs = 0.025))

    plot(0:(length(tau_sdid)-1), tau_sdid, type = 'b', xlim = c(0, 10), ylim = c(-30, 5), xlab = "Lag", ylab = "Estimated effect")
    lines(0:(length(orig_values)-1), orig_values, col = 'blue', lwd = 0.5, lty = 2)
    lines(0:(length(tau_upper)-1), tau_upper, col = 'grey', lwd = 0.5, lty = 2)
    lines(0:(length(tau_lower)-1), tau_lower, col = 'grey', lwd = 0.5, lty = 2)
    legend(0, -25, legend=c("Synth DiD", "Original Results"),
           col=c("black", "blue"), lty=1, cex=0.8)
    lines(-1:(length(tau_upper) - 1), rep(0, length(tau_upper) + 1), col = 'red', lwd = 0.5, lty = 3)

    if (save_pdf) {
      if (length(path) > 0 & path[length(path)] != "/") {
        full_path = paste0(path, "/", "CHC CI.pdf")
      } else {
        full_path = paste0(path, "CHC CI.pdf")
      }

      pdf(, width = 16, height = 9)
      plot(0:(length(tau_sdid)-1), tau_sdid, type = 'b', xlim = c(0, 10), ylim = c(-30, 5), xlab = "Lag", ylab = "Estimated effect")
      lines(0:(length(orig_values)-1), orig_values, col = 'blue', lwd = 0.5, lty = 2)
      lines(0:(length(tau_upper)-1), tau_upper, col = 'grey', lwd = 0.5, lty = 2)
      lines(0:(length(tau_lower)-1), tau_lower, col = 'grey', lwd = 0.5, lty = 2)
      legend(0, -25, legend=c("Synth DiD", "Original Results"),
          col=c("black", "blue"), lty=1, cex=0.8)
      lines(-1:(length(tau_upper) - 1), rep(0, length(tau_upper) + 1), col = 'red', lwd = 0.5, lty = 3)
      dev.off()
    }
}


#' Generate t-statistic plot that compares the emperical t-statistic distribution with the standard normal distribution.
#' @description
#' There is should be no treatment effect since the matrix completion was used on the untreated part of our dataset.
#'
#' @param t_stat_array Numeric array. The array containing t-statistic for all lags.
#' @param lags         Integer or integer vector. Lag number: starts with 0 which stands for the first element of the estimate vector.
#' @param type         Character. The type of the estimator: "sdid" or "did".
#' @param sim_name     Character. The simulation name: "CHC" or "CPS".
#' @param signal_share Numeric. The ratio of the interactive FE variance to the noise variace for 'CHC' simulation.
#' @param breaks       Integer. Number of breaks used in histagram.
#' @param save_pdf     Bool. If TRUE, saves pdf beside generating a plot.
#' @param path         Character.
#' @export
plot_t_stat <- function(
    t_stat_array,
    lags = 0,
    type = "sdid",
    sim_name = "CHC",
    signal_share = 0,
    breaks = 25,
    width = 9,
    height = 9,
    save_pdf = FALSE,
    path = NULL
) {
  if (!type %in% c("sdid", "did")) {
    stop("The type should be either 'sdid' or 'did'.")
  }
  # Change the type name
  type_list <- list(sdid = "SSDiD", did = "DiD")
  type <- type_list[[type]]

   if (!(sim_name %in% c("CHC", "CPS"))) {
    stop("sim_name should be either 'CHC' or 'CPS'.")
  }

  lags <- c(lags)
  signal_perc <- signal_share * 100
  num_digits <- if (signal_perc > 0.5) floor(log10(round(signal_perc))) + 1 else 1
  if (sim_name == "CHC") {
    signal_char <- c(
      paste0(", ", as.character(signal_perc)[1:num_digits][1], " % signal"),
      paste0(", ", as.character(signal_perc)[1:num_digits][1], " percent signal")
    )
  } else {
    signal_char <- c("", "")
  }

  for (lag_show in lags) {
    hist(
      t_stat_array[lag_show + 1, ], breaks = breaks, probability = TRUE,
      xlab = "t-statistic",
      main = paste0("Histogram of ", type, " t-statistic", signal_char[1])
    )
    curve(
      dnorm(x, mean = 0, sd = 1),
      col = "darkblue", lwd = 2, add = TRUE, yaxt = "n"
    )
    if (save_pdf) {
      pdf_name <- paste0(
        "t-statistic ", type, " distibution ",
        lag_show," lag ", signal_char[2], ".pdf"
      )
      #if (length(path) > 0 & path[length(path)] != "/") {
      #  full_path = paste0(path, "/", pdf_name)
      #} else {
      #  full_path = paste0(path, pdf_name)
      #}
      pdf(pdf_name, width = width, height = height)
      hist(
        t_stat_array[lag_show + 1, ], breaks = breaks, probability = TRUE,
        xlab = "t-statistic",
        main = paste0("Histogram of ", type, " t-statistic", signal_char[1])
      )
      curve(
        dnorm(x, mean = 0, sd = 1),
        col = "darkblue", lwd = 2, add = TRUE, yaxt = "n"
      )
      dev.off()
    }
  }
}


#' Plot the empirical CDF of simuated adoption dates.
#'
#' @param cohort_array  Numeric array. The array N x M contains adoption dates in M simulations.
#' @param smoothed      Bool. If TRUE, plots smoothed ECDF.
#' @param T             Integer. The total number of time periods.
#' @param save_pdf      Bool. If TRUE, saves pdf beside generating a plot.
#' @export plot_cohorts
plot_cohorts <- function(
  cohort_array,
  smoothed,
  T = 40,
  save_pdf = FALSE
) {
  cohorts <- as.vector(cohort_array)
  emp_cdf <- stats::ecdf(cohorts)
  x  <- seq(0, T, length.out = T + 1)
  title_add <- ""
  if (smoothed) {
    y <- emp_cdf(x)
    emp_cdf <- approxfun(x, y, method = "linear")
    title_add <- "Smoothed "
  }
  plot(x, emp_cdf(x), type = "b", col = "blue", lwd = 2,
    main = paste0(title_add, "ECDF of adoption dates"),
    xlab = "Adoption dates", ylab = "Cumulative Probability",
	  ylim = c(0, 1), xlim = c(0, T))
  abline(h = 0, lty = 2, col = "grey")
  abline(h = 1, lty = 2, col = "grey")
  abline(h = emp_cdf(T), lty = 2, col = "grey")
  text(
    x = -1, y = emp_cdf(T) + 0.01,
    labels = substr(as.character(emp_cdf(T)), 1, 5),
    pos = 4, col = "black", cex = 1
  )
  if (save_pdf) {
    pdf(
      paste0(title_add, "ECDF of adoption dates.pdf"),
      width = 9, height = 9
    )
    plot(x, emp_cdf(x), type = "b", col = "blue", lwd = 2,
      main = paste0(title_add, "ECDF of adoption dates"),
      xlab = "Adoption dates", ylab = "Cumulative Probability",
      ylim = c(0, 1), xlim = c(0, T)
    )
    abline(h = 0, lty = 2, col = "grey")
    abline(h = 1, lty = 2, col = "grey")
    abline(h = emp_cdf(T), lty = 2, col = "grey")
    text(
      x = -1, y = emp_cdf(T) + 0.01,
      labels = substr(as.character(emp_cdf(T)), 1, 5),
      pos = 4, col = "black", cex = 1
    )
    dev.off()
  }

}
