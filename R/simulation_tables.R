#' Creates 2 tables: RMSE & Bias table and coverage table.
#' @param tau_array_sdid Numeric array. The max_lag x M array contains `sdid` estimates from M simulations.
#' @param tau_array_did  Numeric array. The max_lag x M array contains `did` estimates from M simulations.
#' @param t_stat_sdid    Numeric array. The max_lag x M array contains t-statistic of `sdid` estimate from M simulations.
#' @param t_stat_did     Numeric array. The max_lag x M array contains t-statistic of `did` estimate from M simulations.
#' @param sim_name       Character. The simulation name: "CHC" or "CPS".
#' @param signal_share   Numeric. The ratio of the interactive FE variance to the noise variace for 'CHC' simulation.
#' @param save_tables    Bool. If TRUE, the performance and coverage tables are saved.
#' @param path           Character.
#' @export
create_tables <- function(
    tau_array_sdid,
    tau_array_did,
    t_stat_sdid,
    t_stat_did,
    sim_name = "CHC",
    signal_share = 0.0,
    save_tables = FALSE,
    path = ""
) {
  if (!(sim_name %in% c("CHC", "CPS"))) {
    stop("sim_name should be either 'CHC' or 'CPS'.")
  }
  rmse_sdid <- apply(tau_array_sdid, 1, FUN = function(x) sqrt(mean(x^2)))
  rmse_did <- apply(tau_array_did, 1, FUN = function(x) sqrt(mean(x^2)))

  bias_sdid <- apply(tau_array_sdid, 1, FUN = function(x) mean(x))
  bias_did <- apply(tau_array_did, 1, FUN = function(x) mean(x))

  bias_to_sd_sdid <- bias_sdid / apply(tau_array_sdid, 1, FUN = function(x) sd(x))
  bias_to_sd_did <- bias_did / apply(tau_array_did, 1, FUN = function(x) sd(x))

  perf_table <- rbind(
    rmse_sdid, rmse_did, rmse_sdid / rmse_did,
    bias_sdid, bias_did, bias_sdid / bias_did,
    bias_to_sd_sdid, bias_to_sd_did, bias_to_sd_sdid / bias_to_sd_did
  )
  rownames(perf_table) <- c(
    "SSDiD RMSE", "DiD RMSE", "SSDiD to DiD RMSE ratio",
    "SSDiD Bias", "DiD Bias", "SSDiD to DiD Bias ratio",
    "SSDiD Bias/SD", "DiD Bias/SD", "SSDiD to DiD Bias/SD"
  )
  colnames(perf_table) <-  as.character(0:(length(rmse_did) - 1))

  signal_perc <- signal_share * 100

  num_digits <- if (signal_perc > 0.5) floor(log10(round(signal_perc))) + 1 else 1
  if (sim_name == "CHC") {
    signal_char <- c(
      paste0(", ", as.character(signal_perc)[1:num_digits], " % signal CHC"),
      paste0(", ", as.character(signal_perc)[1:num_digits], " percent signal")
    )
  } else {
    signal_char <- c(", CPS", "")
  }


  # RMSE and Bias table
  if (save_tables) {
    print(
      xtable::xtable(perf_table, caption = paste0("RMSE and Bias table", signal_char[1])),
      file = paste0(path, "RMSE and Bias table", signal_char[2], ".tex")
    )
  }
  print(
    xtable::xtable(perf_table, caption = paste0("RMSE and Bias table", signal_char[1]))
  )

  cat("\n\n")
  # Coverage table
  coverage_sdid <- apply(abs(t_stat_sdid) < qnorm(0.975), 1, mean)
  coverage_did <- apply(abs(t_stat_did) < qnorm(0.975), 1, mean)

  cov_table <- rbind(coverage_sdid, coverage_did)
  rownames(cov_table) <- c("SSDiD", "DiD")
  colnames(cov_table) <- as.character(0:(length(coverage_did) - 1))
  if (save_tables){
    print(
      xtable::xtable(cov_table, caption = paste0(
        "Coverage table of estimators", signal_char[1])
      ),
      file = paste0(path, "Coverage table", signal_char[2], ".tex")
    )
  }
  print(xtable::xtable(
    cov_table,
    caption = paste0("Coverage table of estimators", signal_char[1])
  )
  )
}
