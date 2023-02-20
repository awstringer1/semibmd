/// @file ModelA.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
// (here it's ModelA)
template <class Type>
Type monotonesmoothing(objective_function<Type>* obj) {
  using namespace density;
  // Gaussian monotone GAM
  DATA_VECTOR(y); // Response
  DATA_MATRIX(Xmono); // Basis function matrix, dimension n x d
  DATA_VECTOR(c); // colSums(X), pass in for speed. Length d
  DATA_SPARSE_MATRIX(Dpmono); // Penalty matrix, dimension r x r
  DATA_MATRIX(Umono); // Eigenvectors of penalty, used to transform beta

  // PARAMETER(alpha); // intercept
  PARAMETER_VECTOR(betaRmono); // length r
  PARAMETER_VECTOR(betaFmono); // length d-r
  PARAMETER(alpha);


  PARAMETER(logprec); // log(1/sigma^2), log precision of y
  PARAMETER(logsmoothingmono); // log(lambda), log smoothing param

  // Constants
  int n = y.size();
  int r = betaRmono.size();
  int d = r + betaFmono.size();

  // Transformations
  vector<Type> betaC(d);
  betaC << betaRmono,betaFmono;
  vector<Type> beta = Umono*betaC;


  Type sd = exp(-0.5*logprec);
  Type lambda = exp(logsmoothingmono);
  // Dp --> lambda * Dp; have to loop over nonzero elements
  for (int k=0; k<Dpmono.outerSize(); k++)
    for (typename SparseMatrix<Type>::InnerIterator it(Dpmono,k); it; ++it)
      it.valueRef() = it.value() * lambda;

  vector<Type> gamma(d);
  gamma(0) = beta(0);
  for (int i=1;i<d;	i++)
    gamma(i) = gamma(i-1) - exp(beta(i));

  vector<Type> mu = Xmono*gamma;
  Type colsum = (c * gamma).sum();
  for (int i=0;i<n;i++) {
    mu(i) -= colsum / n;
    mu(i) += alpha;
  }

  // Objective function
  // Negative "log likelihood + log prior"
  // Initialize with prior
  Type logprior = GMRF(Dpmono)(betaRmono); // Negative! TODO: do this manually

  Type loglik = 0.;
  for (int i=0;i<n;i++)
    // loglik -= dnorm(y(i),mu(i),sd,true);
    loglik += 0.5*exp(logprec)*(y(i)-mu(i))*(y(i)-mu(i)) - 0.5*logprec; // This has been checked

  return logprior + loglik; // Actually minus loglik
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this


