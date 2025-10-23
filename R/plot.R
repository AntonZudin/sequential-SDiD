#' Plot estimated treatment effect by lag.
#'
#' @import ggplot2
#' @param object :          `sequential_estimate` class object.
#' @param se :              List or numeric vector.
#' @param lags_show         Integer or NULL. The number of lags to show.  If NULL, all lags are shown.
#' @param sig_level :       Numeric. Significance level for the confidence interval.
#' @param color :           String. The color to use for the plot.
#' @param title_font_size : Integer. The font size of the title.
#' @export

plot.sequential_estimate <- function(
  object,
  se,
  lags_show = NULL,
  sig_level = 0.05,
  color = c("#e87d72","#56bcc2"),
  ylim = NULL,
  xlab = "Lags",
  ylab = NULL,
  plot_title = "",
  error_bar_width = 0.1,
  point_size = 1.2,
  xgap = 1,
  legend = TRUE,
  ref_line = 0,
  theming = TRUE,
  title_font_size = 16,
  legend_position = "bottom",
  ...
) {

  if (is.null(dim(object))) {
    n_lags <- length(object[["tau_sdid"]])
  } else {
    # Only estimator is present
    n_lags <- dim(object)
    # Only one color is needed
    color <- color[1]
  }

  if (is.null(lags_show)) {
    lags_show <- n_lags
  } else if (n_lags < lags_show) {
    stop("The number of lags to show should not exceed the number of lags in the data.")
  }

  type <- attr(object, "type")
  if (type == "both") {
    types <- c("DiD", "SDiD")
  } else {
    types <- c(type)
  }

  results <- data.frame(
    treat_lag = c(),
    att = c(),
    lower_bound = c(),
    upper_bound = c(),
    type = c()
  )

  for (est_type in types) {
    tau_lag <- object[[paste0("tau_", tolower(est_type))]][1:lags_show]
    if ("tau_b" %in% names(se)) {
      tau_b <- se[["tau_b"]]
      if (type == "both") tau_b <- tau_b[[paste0("tau_b_", tolower(est_type))]]
      lower_bound <- apply(
        tau_b, 1,
        FUN = function(x) quantile(x, probs = 0.025)
      )[1:lags_show]
      upper_bound <- apply(
        tau_b, 1,
        FUN = function(x) quantile(x, probs = 0.975)
      )[1:lags_show]
    } else {
      if (type == "both") se_estim <- se[[paste0("se_", tolower(est_type))]]
      se_estim <- se_estim[1:lags_show]
      crit_value <- qnorm(1 - sig_level / 2)
      lower_bound <- (tau_lag - crit_value * se_estim)
      upper_bound <- (tau_lag + crit_value * se_estim)
    }

    results <- rbind(
      results,
      data.frame(
        treat_lag = 0:(lags_show - 1),
        att = tau_lag,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        type = rep(est_type, lags_show)
      )
    )
  }

  if (plot_title == "") {
    # get title right depending on which aggregation
    plot_title <- ("Treatment Effect by Lag")
  }

  plt <- ggplot(results,
                aes(x = treat_lag, y = att, colour = type,
                  ymin = lower_bound,
                  ymax = upper_bound)
                ) +
    geom_point(
      size = 1.5,
      position = position_dodge(width = 0.5)) +
    geom_errorbar(
      width = error_bar_width,
      position = position_dodge(width = 0.5)) +
    scale_x_continuous(breaks = seq(0, lags_show, xgap)) +
    labs(x = xlab, y = ylab, title = plot_title, color = NULL) +
    scale_color_manual(values = color)

  if (!is.null(ref_line)) {
    plt <- plt + geom_hline(aes(yintercept = ref_line), linetype = 'dashed')
  }

  if (theming) {
    plt <- plt + theme_classic() +
      theme(plot.title = element_text(color = "darkgrey", face="bold", size = title_font_size),
            axis.title = element_text(color = "black", size = 12),
            strip.background = element_rect(fill = 'white', color = 'white'),
            strip.text = element_text(color = 'darkgray', face = 'bold', size = 12, hjust = 0),
            legend.position = legend_position,
            legend.margin = margin(0, 0, 0, -10))
  }

  plt
}
