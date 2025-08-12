#' Plot average estimated effect.
#' @import ggplot2
#' @param object :    `sequential_estimate` class object.
#' @param se :        List or numeric vector.
#' @param sig_level : Numeric. Significance level for the confidence interval.
#' @export

plot.sequential_estimate <- function(
  object,
  se,
  sig_level = 0.95,
  color = "#e87d72",
  ylim=NULL,
  xlab=NULL,
  ylab=NULL,
  title="",
  xgap=1,
  legend=TRUE,
  ref_line = 0,
  theming = TRUE,
  ...
) {
  type <- attr(object, "type")
  if (type == "both") {
    #tau_lag <-
  } else {
    tau_lag <- object[1:length(object)]
  }

  # TODO: for all cases

  results <- cbind.data.frame(
    treat_lag = 0:(length(tau_lag) - 1),
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

    geom_errorbar(width = 0.1, color = color) +
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
