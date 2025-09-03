#' Plot average estimated effect.
#' @import ggplot2
#' @param object :    `sequential_estimate` class object.
#' @param se :        List or numeric vector.
#' @param lags_show   Integer or NULL. The number of lags to show.  If NULL, all lags are shown.
#' @param sig_level : Numeric. Significance level for the confidence interval.
#' @param color :     String. The color to use for the plot.
#' @export

plot.sequential_estimate <- function(
  object,
  se,
  lags_show = NULL,
  sig_level = 0.95,
  color = "#e87d72",
  ylim=NULL,
  xlab=NULL,
  ylab=NULL,
  title="",
  error_bar_width = 0.1,
  xgap=1,
  legend=TRUE,
  ref_line = 0,
  theming = TRUE,
  ...
) {

  if (is.null(lags_show)) {
    lags_show <- length(object)
  } else if (length(object) < lags_show) {
    stop("The number of lags to show should not exceed the number of lags in the data.")
  }

  type <- attr(object, "type")
  if (type == "both") {
    #tau_lag <-
  } else {
    tau_lag <- object[1:lags_show]
    se <- se[1:lags_show]
  }

  # TODO: for all cases

  results <- cbind.data.frame(
    treat_lag = 0:(lags_show - 1),
    att = tau_lag,
    att_se = se
  )
  crit_value <- qnorm(1 - (1 - sig_level)/2)


  if (title == "") {
    # get title right depending on which aggregation
    title <- ("Average Effect by Length of Exposure")
  }

  plt <- ggplot(results,
              ggplot2::aes(x = treat_lag, y = att,
                  ymin = (att - crit_value * att_se),
                  ymax = (att + crit_value * att_se))) +

    geom_errorbar(width = error_bar_width, color = color) +
    scale_x_continuous(breaks = seq(0, lags_show, xgap)) +
    labs(x = xlab, y = ylab, title = title, color = NULL)

  if (!is.null(ref_line)) {
    plt <- plt + geom_hline(aes(yintercept = ref_line), linetype = 'dashed')
  }

  if (theming) {
    plt <- plt + theme_classic() +
      theme(plot.title = element_text(color="darkgray", face="bold", size=12),
            axis.title = element_text(color="darkgray", face="bold", size=12),
            strip.background = element_rect(fill = 'white', color = 'white'),
            strip.text = element_text(color = 'darkgray', face = 'bold', size = 12, hjust = 0),
            legend.position = 'bottom')
  }

  plt
}
