/// @file ModelA.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
// (here it's ModelA)
template <class Type>
Type monotonesmoothingmulti(objective_function<Type>* obj) {
  using namespace density;
  // Gaussian monotone GAM
  DATA_VECTOR(y); // Response, n = y.size()
  DATA_MATRIX(Xmono); // Basis function matrix for monotone spline, dimension n x dm. eta = alpha + Xmono*gamma + ...
  DATA_SPARSE_MATRIX(Dpmono); // Penalty matrix for monotone smooths, dimension rm x rm
  DATA_MATRIX(Umono); // Eigenvectors of penalty, used to transform beta

  PARAMETER_VECTOR(betaRmono); // length sum(r)
  PARAMETER_VECTOR(betaFmono); // length sum(d-r)
  PARAMETER(alpha); // Global intercept
  PARAMETER(logprec); // log(1/sigma^2), log precision of y
  PARAMETER(logsmoothingmono); // log(lambda), log smoothing param

  // If there are non-monotone smooths, add them now
  DATA_INTEGER(nonmono);
  DATA_MATRIX(Xsmooth); // Basis function matrix for non-monotone smooths, dimension n x sum(ds). eta = ... + Xsmooth*betasmooth
  DATA_SPARSE_MATRIX(Dpsmooth); // Penalty matrix for monotone smooths, dimension rs x rs
  DATA_SPARSE_MATRIX(Usmooth); // Eigenvectors of penalty, used to transform beta
  PARAMETER_VECTOR(betaRsmooth); // length sum(rs)
  PARAMETER_VECTOR(betaFsmooth); // length sum(ds-rs)
  PARAMETER_VECTOR(logsmoothingsmooth); // log(lambda), log smoothing params for non-mono smooths, length = length(rs)

  // Constants
  int n = y.size();
  int r = betaRmono.size();
  int d = r + betaFmono.size();
  int ms = logsmoothingsmooth.size(); // Number of non-mono smooth terms
  // Dimensions have to be provided
  DATA_IVECTOR(rs); // Integer vector of dimensions of the penalized params for the non-mono smooths
  DATA_IVECTOR(ds); // Integer vector of dimensions of the full params for the non-mono smooths
  int dsm = ds.sum(), rsm=rs.sum();

  // Transformations
  vector<Type> betaCmono(d);
  betaCmono << betaRmono,betaFmono;
  vector<Type> betamono = Umono*betaCmono;
  vector<Type> gamma(d);
  gamma(0) = betamono(0);
  for (int i=1;i<d; i++)
    gamma(i) = gamma(i-1) - exp(betamono(i));
  // betaR and betaF contain all the smooth terms' coefficients sequentially
  // Segment them out into a full beta
  vector<Type> betaCsmooth(Xsmooth.cols());
  int tmpR=0,tmpF=0,tmp=0;
  vector<Type> betasmooth(betaCsmooth.size());
  betasmooth.setZero();
  // std::cout << "betaCsmooth.size() = " << betaCsmooth.size() << std::endl;
  if (nonmono==1) {
    for (int i=0;i<rs.size();i++) {
      // std::cout << "R iteration start " << tmp << ", length " << rs(i) << ", index " << tmpR << std::endl;
      // std::cout << "F iteration start " << tmp+rs(i) << ", length " << ds(i)-rs(i) << ", index " << tmpF << std::endl;
      betaCsmooth.segment(tmp,rs(i)) = betaRsmooth.segment(tmpR,rs(i));
      betaCsmooth.segment(tmp+rs(i),ds(i)-rs(i)) = betaFsmooth.segment(tmpF,ds(i)-rs(i));
      tmpR += rs(i);
      tmpF += ds(i)-rs(i);
      tmp += ds(i);
    }
    betasmooth = Usmooth * betaCsmooth;
  }
  Type sd = exp(-0.5*logprec);
  Type lambdamono = exp(logsmoothingmono);
  // Dp --> lambda * Dp; have to loop over nonzero elements
  for (int k=0; k<Dpmono.outerSize(); k++)
    for (typename SparseMatrix<Type>::InnerIterator it(Dpmono,k); it; ++it)
      it.valueRef() = it.value() * lambdamono;

  int tmpidx=0;
  int tmpcounter = rs(tmpidx);
  if (nonmono==1) {
    // Loop over the nonzero (diagonal) elements of Dpsmooth and multiply them by the appropriate element of lambda
    for (int k=0; k<Dpsmooth.outerSize(); k++)
      for (typename SparseMatrix<Type>::InnerIterator it(Dpsmooth,k); it; ++it) {
        // std::cout << "tmpcounter = " << tmpcounter << "| tmpidx = " << tmpidx << std::endl;
        it.valueRef() = it.value() * exp(logsmoothingsmooth(tmpidx));
        // std::cout << "Dp[j,j] = " << it.value() << std::endl;
        tmpcounter--;
        if (tmpcounter == 0) {
          tmpidx++;
          if (tmpidx < rs.size())
            tmpcounter = rs(tmpidx);
        }
      }
  }

  vector<Type> mu = Xmono*gamma;
  if (nonmono==1)
    mu += Xsmooth*betasmooth;
  Type mumean = mu.mean();
  for (int i=0;i<n;i++) {
    mu(i) -= mumean;
    mu(i) += alpha;
  }

  // Objective function
  // Negative "log likelihood + log prior"
  // Initialize with prior
  Type logprior = GMRF(Dpmono)(betaRmono); // Negative! Note: doing manually not faster
  if (nonmono==1)
    logprior += GMRF(Dpsmooth)(betaRsmooth);

  Type loglik = 0.;
  for (int i=0;i<n;i++)
    // loglik -= dnorm(y(i),mu(i),sd,true);
    loglik += 0.5*exp(logprec)*(y(i)-mu(i))*(y(i)-mu(i)) - 0.5*logprec; // This has been checked

  return logprior + loglik; // Actually minus loglik
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this


