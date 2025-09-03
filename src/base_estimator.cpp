#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

//' Penalty for SSDiD
 //' @description
 //' The regulazation term eta^2 for SSDiD estimator without s2.
 //' @param N    Double.
 //' @param deg  Double. The default value of 0.9 is proposed in the paper in Remark 3.2.
 //' @return     Numeric. The penalty term without s2.
 double penalty(double N, double pen = 0.9) {
   return std::pow(N, -pen);
 }

 // [[Rcpp::export]]
Rcpp::NumericVector base_estimator(
     const arma::mat& Y,
     const arma::vec& n_j,
     double penalty, // eta^2
     double s2 = -1.0, // Substitute for NULL
     std::string type = "sdid"
 ) {
   int j_c = Y.n_rows - 1;
   int t_c = Y.n_cols - 1;

   // Check dimensions
   if (static_cast<int>(Y.n_rows) != static_cast<int>(n_j.n_elem)) {
     Rcpp::stop("The number of cohorts does not coincide in Y matrix and n_j vector");
   }
   if (type != "did" && type != "sdid") {
     Rcpp::stop("The 'type' argument should be either 'sdid' or 'did'");
   }

   // s2 should be either -1.0 (NULL) or positive
   if (!(s2 == -1.0 || s2 > 0)) {
     Rcpp::stop("s2 should be either NULL (-1.0) or positive");
   }

   double N = arma::sum(n_j.subvec(0, j_c - 1));
   arma::vec pi = n_j.subvec(0, j_c - 1) / N;

   // Extract Y_c (j_c x t_c)
   arma::mat Y_c = Y.submat(0, 0, j_c - 1, t_c - 1);
   arma::vec Y_j0 = Y.submat(j_c, 0, j_c, t_c - 1).t();
   arma::vec Y_t0 = Y.submat(0, t_c, j_c - 1, t_c);
   double Y_j0_t0 = Y(j_c, t_c);

   arma::vec lambda_reg, gamma_reg;

   if (type == "did") {
     lambda_reg = arma::ones<vec>(t_c) / t_c;
     gamma_reg = pi;
   } else {
     arma::vec ones_t = arma::ones<vec>(t_c);
     arma::vec ones_j = arma::ones<vec>(j_c);

     // Gamma - unit weights
     arma::mat Sigma_tc = diagmat(penalty / pi); // It is not exactly a cov matix
     arma::vec grad_gamma1 = 2 * (Y_c * (-Y_j0));
     double grad_gamma2 = 2 * arma::sum(-Y_j0);

     arma::vec grad_reg = arma::zeros(j_c + 2);
     grad_reg.subvec(0, j_c - 1) = grad_gamma1;
     grad_reg(j_c) = grad_gamma2;
     grad_reg(j_c + 1) = -1;

     arma::mat block_1 = 2 * Y_c * Y_c.t() + 2 * Sigma_tc;
     arma::vec block_2 = 2 * (Y_c * ones_t);

     arma::mat hess_reg = arma::zeros(j_c + 2, j_c + 2);
     hess_reg.submat(0, 0, j_c - 1, j_c - 1) = block_1;
     hess_reg.submat(0, j_c, j_c - 1, j_c) = block_2;
     hess_reg.submat(0, j_c + 1, j_c - 1, j_c + 1) = ones_j;
     hess_reg.submat(j_c, 0, j_c, j_c - 1) = block_2.t();
     hess_reg(j_c, j_c) = 2 * t_c;
     hess_reg.submat(j_c + 1, 0, j_c + 1, j_c - 1) = ones_j.t();

     // This prevents the code from breaking,
     // forcing to use DiD when hessian of the objective function is singular.
     if (arma::rcond(hess_reg) < 2.5e-16) {
       Rcpp::warning("The hessian for unit weights is near singular. Doing DiD unit weights.");
       gamma_reg = pi;
     } else {
       arma::vec solution = arma::solve(hess_reg, -grad_reg, solve_opts::likely_sympd);
       gamma_reg = solution.subvec(0, j_c - 1);
     }

     // Lambda - time weights
     double diag_val_jc = penalty * (1.0 / j_c) * arma::sum(1.0 / pi);
     arma::mat Sigma_jc = diag_val_jc * arma::eye(t_c, t_c);

     arma::vec grad_lambda_1 = 2 * (Y_c.t() * (-Y_t0));
     double grad_lambda_2 = 2 * arma::sum(-Y_t0);

     arma::vec gradl_reg = arma::zeros(t_c + 2);
     gradl_reg.subvec(0, t_c - 1) = grad_lambda_1;
     gradl_reg(t_c) = grad_lambda_2;
     gradl_reg(t_c + 1) = -1;

     arma::mat blockl_1 = 2 * Y_c.t() * Y_c + 2 * Sigma_jc;
     arma::vec blockl_2 = 2 * (Y_c.t() * ones_j);

     arma::mat hessl_reg = arma::zeros(t_c + 2, t_c + 2);
     hessl_reg.submat(0, 0, t_c - 1, t_c - 1) = blockl_1;
     hessl_reg.submat(0, t_c, t_c - 1, t_c) = blockl_2;
     hessl_reg.submat(0, t_c + 1, t_c - 1, t_c + 1) = ones_t;
     hessl_reg.submat(t_c, 0, t_c, t_c - 1) = blockl_2.t();
     hessl_reg(t_c, t_c) = 2 * j_c;
     hessl_reg.submat(t_c + 1, 0, t_c + 1, t_c - 1) = ones_t.t();

     if (arma::rcond(hessl_reg) < 2.5e-16) {
       Rcpp::warning("The hessian for time weights is near singular. Doing DiD time weights");
       lambda_reg = arma::ones<vec>(t_c) / t_c;
     } else {
       arma::vec solution = arma::solve(hessl_reg, -gradl_reg, solve_opts::likely_sympd);
       lambda_reg = solution.subvec(0, t_c - 1);
     }
   }

   double term1 = Y_j0_t0 - arma::dot(Y_j0, lambda_reg);
   double term2 = arma::dot(Y_t0, gamma_reg) - arma::as_scalar(gamma_reg.t() * Y_c * lambda_reg);
   double tau = term1 - term2;

   Rcpp::NumericVector result = Rcpp::NumericVector(2);

   if (s2 > 0) {
     double est_variance = s2 *
       (1/ pi(j_c) + arma::sum(gamma_reg / pi.subvec(0, j_c - 1))) *
       (1 + arma::sum(lambda_reg));
     result[0] = tau;
     result[1] = est_variance;
   } else {
     result[0] = tau;
     result[1] = -1.0;
   }

   return result;
 }
