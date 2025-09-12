#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

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



//' Computes the base synthetic diff-in-diff or diff-in-diff estimate.
//' @description
//' The armadillo implementation of `base_estimator` algorithm.
//' The bottom right cell is the only cell being treated: W_it = 1.
//' @param Y :            Numeric arma::mat. A submatrix of outcomes with one treated obs in the bottom right corner.
//' @param coh_weights :  Numeric arma::vec. The vector of cohort weights. It might be the number of units in cohorts or
//'                                          the cohort weights, usually just cohort population.
//' @param penalty :      double.
//' @param s2 :           double or -1.0. The upper bound estimate of noise variance.
//' @param type :         std::string. Type of the estimator should be `sdid` or `did`.
//' @param fast:          bool. If true, fast matrix inversion is conducted.
//'
//' @return               NumericVector size of 2. The first element is the treatment effect,
//'                                         the second one is the asymptotic variance.
// [[Rcpp::export]]
Rcpp::NumericVector base_estimator(
     const arma::mat& Y,
     const arma::vec& coh_weights,
     double coh_weight_sum,
     double penalty, // eta^2
     double s2 = -1.0, // -1.0 is a substitute for NULL
     std::string type = "sdid",
     bool fast = false
 ) {
   int j_c = Y.n_rows - 1;
   int t_c = Y.n_cols - 1;

   // Check dimensions
   if (static_cast<int>(Y.n_rows) != static_cast<int>(coh_weights.n_elem)) {
     Rcpp::stop("The number of cohorts does not coincide in Y matrix and coh_weights vector");
   }
   if (type != "did" && type != "sdid") {
     Rcpp::stop("The 'type' argument should be either 'sdid' or 'did'");
   }

   // s2 should be either -1.0 (NULL) or positive
   if (!(s2 == -1.0 || s2 > 0)) {
     Rcpp::stop("s2 should be either NULL (-1.0) or positive");
   }

   double coh_weights_sum = arma::sum(coh_weights);
   arma::vec pi = coh_weights.subvec(0, j_c) / coh_weights_sum;
   double coh_weight_contr = arma::sum(coh_weights.subvec(0, j_c - 1));
   arma::vec adjusted_pi = coh_weights.subvec(0, j_c - 1) / coh_weight_contr; // This pi sums to 1

   // Extract Y_c (j_c x t_c)
   arma::mat Y_c = Y.submat(0, 0, j_c - 1, t_c - 1);
   arma::vec Y_j0 = Y.submat(j_c, 0, j_c, t_c - 1).t();
   arma::vec Y_t0 = Y.submat(0, t_c, j_c - 1, t_c);
   double Y_j0_t0 = Y(j_c, t_c);

   arma::vec lambda_weights, gamma_weights;

   if (type == "did") {
     lambda_weights = arma::ones<vec>(t_c) / t_c;
     gamma_weights = adjusted_pi.subvec(0, j_c - 1);
   } else {
     arma::vec ones_t = arma::ones<vec>(t_c);
     arma::vec ones_j = arma::ones<vec>(j_c);

     // Gamma - unit weights
     arma::mat Sigma_tc = diagmat(penalty / adjusted_pi); // It isn't exactly a cov matrix
     // Parts of gradient for gamma optimization problem
     arma::vec grad_gamma1 = 2 * (Y_c * (-Y_j0));
     double grad_gamma2 = 2 * arma::sum(-Y_j0);

     arma::vec gamma_grad = arma::zeros(j_c + 2);
     gamma_grad.subvec(0, j_c - 1) = grad_gamma1;
     gamma_grad(j_c) = grad_gamma2;
     gamma_grad(j_c + 1) = -1;

     arma::mat block_1 = 2 * Y_c * Y_c.t() + 2 * Sigma_tc;
     arma::vec block_2 = 2 * (Y_c * ones_t);

     arma::mat gamma_hess = arma::zeros(j_c + 2, j_c + 2);
     gamma_hess.submat(0, 0, j_c - 1, j_c - 1) = block_1;
     gamma_hess.submat(0, j_c, j_c - 1, j_c) = block_2;
     gamma_hess.submat(0, j_c + 1, j_c - 1, j_c + 1) = ones_j;
     gamma_hess.submat(j_c, 0, j_c, j_c - 1) = block_2.t();
     gamma_hess(j_c, j_c) = 2 * t_c;
     gamma_hess.submat(j_c + 1, 0, j_c + 1, j_c - 1) = ones_j.t();

     // This prevents the code from breaking,
     // forcing to use DiD when hessian of the objective function is singular.
     if (arma::rcond(gamma_hess) < 2.5e-16) {
       Rcpp::warning("The hessian for unit weights is near singular. Doing DiD unit weights.");
       gamma_weights = adjusted_pi;
     } else {
       arma::vec solution;

       if (fast) {
          solution = arma::solve(gamma_hess, -gamma_grad, solve_opts::fast);
       } else {
          solution = arma::solve(gamma_hess, -gamma_grad, solve_opts::likely_sympd);
       }
       gamma_weights = solution.subvec(0, j_c - 1);
     }

     // Lambda - time weights
     double diag_val_jc = penalty * (1.0 / j_c) * arma::sum(1.0 / adjusted_pi);
     arma::mat Sigma_jc = diag_val_jc * arma::eye(t_c, t_c);

     arma::vec grad_lambda_1 = 2 * (Y_c.t() * (-Y_t0));
     double grad_lambda_2 = 2 * arma::sum(-Y_t0);

     arma::vec lambda_grad = arma::zeros(t_c + 2);
     lambda_grad.subvec(0, t_c - 1) = grad_lambda_1;
     lambda_grad(t_c) = grad_lambda_2;
     lambda_grad(t_c + 1) = -1;

     arma::mat blockl_1 = 2 * Y_c.t() * Y_c + 2 * Sigma_jc;
     arma::vec blockl_2 = 2 * (Y_c.t() * ones_j);

     arma::mat lambda_hess = arma::zeros(t_c + 2, t_c + 2);
     lambda_hess.submat(0, 0, t_c - 1, t_c - 1) = blockl_1;
     lambda_hess.submat(0, t_c, t_c - 1, t_c) = blockl_2;
     lambda_hess.submat(0, t_c + 1, t_c - 1, t_c + 1) = ones_t;
     lambda_hess.submat(t_c, 0, t_c, t_c - 1) = blockl_2.t();
     lambda_hess(t_c, t_c) = 2 * j_c;
     lambda_hess.submat(t_c + 1, 0, t_c + 1, t_c - 1) = ones_t.t();

     if (arma::rcond(lambda_hess) < 2.5e-16) {
       Rcpp::warning("The hessian for time weights is near singular. Doing DiD time weights");
       lambda_weights = arma::ones<vec>(t_c) / t_c;
     } else {
       arma::vec solution;

       if (fast) {
          solution = arma::solve(lambda_hess, -lambda_grad, solve_opts::fast);
       } else {
          solution = arma::solve(lambda_hess, -lambda_grad, solve_opts::likely_sympd);
       }

       lambda_weights = solution.subvec(0, t_c - 1);
     }
   }

   double term1 = Y_j0_t0 - arma::dot(Y_j0, lambda_weights);
   double term2 = arma::dot(Y_t0, gamma_weights) - arma::as_scalar(gamma_weights.t() * Y_c * lambda_weights);
   double tau = term1 - term2;

   Rcpp::NumericVector result(2);

   if (s2 > 0.0) {
     // TODO: Make sure that I use the correct pi
     double est_variance = s2 *
       (1.0 /pi(j_c) + arma::sum(arma::pow(gamma_weights.subvec(0, j_c - 1), 2) / pi.subvec(0, j_c - 1))) *
       (1.0 + arma::sum(arma::pow(lambda_weights.subvec(0, t_c - 1), 2)));
     result[0] = tau;
     result[1] = est_variance;
   } else {
     result[0] = tau;
     result[1] = -1.0;
   }

   return result;
 }
