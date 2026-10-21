
#include <TMB.hpp>
// Hintzen et al. 2019, Ec. 1 — NB truncada en 0
// nb2 = 1: Var = mu + beta mu^2   (size = 1/beta)    [NB2]
// nb2 = 0: Var = mu (1 + beta)    (size = mu/beta)   [NB1, paper]
// log(beta_t) = mu_logbeta + eps_t ;  eps_t ~ N(0, sd)  (aleatorio por año)
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_IVECTOR(t);
  DATA_INTEGER(re);
  DATA_INTEGER(nb2);
  PARAMETER_VECTOR(log_mu);
  PARAMETER(mu_logbeta);
  PARAMETER(log_sd_logbeta);
  PARAMETER_VECTOR(eps);
  int nt = log_mu.size();
  Type nll = 0;
  vector<Type> log_beta(nt);
  for (int j = 0; j < nt; j++) log_beta(j) = mu_logbeta + eps(j);
  if (re == 1) for (int j = 0; j < nt; j++)
    nll -= dnorm(eps(j), Type(0), exp(log_sd_logbeta), true);
  for (int i = 0; i < y.size(); i++) {
    Type mu   = exp(log_mu(t(i)));
    Type beta = exp(log_beta(t(i)));
    Type var, logp0;
    if (nb2 == 1) {
      var   = mu + beta * mu * mu;
      logp0 = -(Type(1) / beta) * log(Type(1) + beta * mu);
    } else {
      var   = mu * (Type(1) + beta);
      logp0 = -(mu / beta) * log(Type(1) + beta);
    }
    Type ll = dnbinom2(y(i), mu, var, true);
    ll -= logspace_sub(Type(0), logp0);   // truncamiento: / (1 - P(Y=0))
    nll -= ll;
  }
  vector<Type> beta = exp(log_beta);
  vector<Type> mu = exp(log_mu);
  REPORT(beta); REPORT(mu);
  ADREPORT(log_beta);
  ADREPORT(log_mu);
  return nll;
}

