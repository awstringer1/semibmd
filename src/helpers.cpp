#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// get the 0-based knot index
// [[Rcpp::export]]
int knotindex(double x,Eigen::VectorXd t) {
  int k=0;
  while(x>=t(k))
    k++;
  return k-1;
}

// B-splines
double weight(double x,Eigen::VectorXd t,int i,int k) {
  if (t(i+k-1) != t(i-1))
    return((x - t(i-1))/(t(i+k-1)-t(i-1)));
  return 0.;
}
// [[Rcpp::export]]
double Bspline(double x,int j,Eigen::VectorXd t,int p) {
  // Evaluate the jth B-spline
  // B_p(x) of order p (degree p-1) at x
  if (p==1)
    return(x>=t(j-1) && x<t(j+1-1));

  double w1 = weight(x,t,j,p-1);
  double w2 = weight(x,t,j+1,p-1);
  double b1 = Bspline(x,j,t,p-1);
  double b2 = Bspline(x,j+1,t,p-1);

  return w1*b1 + (1.-w2)*b2;
}
// [[Rcpp::export]]
Eigen::VectorXd Bsplinevec(double x,Eigen::VectorXd t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  int k = knotindex(x,t);
  for (int i=(k-(p-1));i<k+1;i++)
    b(i) = Bspline(x,i+1,t,p);
  return b;
}

// [[Rcpp::export]]
double BsplineD(double x,int j,Eigen::VectorXd t,int p) {
  if (p==1)
    return(0.);

  double out = (p-1) * (Bspline(x,j,t,p-1)/(t(j+(p-1)) - t(j)) - Bspline(x,j+1,t,p-1)/(t(j+p) - t(j+1)));
  return out;
}

// [[Rcpp::export]]
Eigen::VectorXd BsplinevecD(double x,Eigen::VectorXd t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  int k = knotindex(x,t);
  for (int i=(k-(p-1));i<k+1;i++)
    b(i) = BsplineD(x,i+1,t,p);
  return b;
}


// deBoor's algorithm for spline
// [[Rcpp::export]]
double deBoor(double x,int k,Eigen::VectorXd t,Eigen::VectorXd beta,int p) {
  // x: evaluation point
  // k: knot index (NOTE: keep it input, since will need it twice each iteration)
  // t: knots
  // beta: coefficients
  // p: spline ORDER (cubic = 4)
  //
  // outputs f(x) via deBoor's algorithm

  std::cout << "p=" << p << ", beta length=" << beta.size() << ", k=" << k << ", t length=" << t.size() <<  std::endl;
  Eigen::VectorXd d(p);
  for (int j=0;j<p;j++) {
    if (j+k-(p-1) < 0 || j+k-(p-1) > beta.size()-1) std::cout << "Accessing index " << j+k-(p-1) << "of beta" << std::endl; 
    d(j) = beta(j+k-(p-1));
  }

  double a =0.;
  for (int r=1;r<p;r++)
    for (int j=(p-1);j>r-1;j--) {
      if (j+k-(p-1) < 0 || j+k-(p-1) > t.size()-1 || j+1+k-r < 0 || j+1+k-r > t.size()-1) std::cout << "Accessing indices " << j+k-(p-1) << ", " << j+1+k-r << " of t" << std::endl;	    
      a = (x-t(j+k-(p-1))) / (t(j+1+k-r) - t(j+k-(p-1)));
      d(j) = (1. - a)*d(j-1) + a*d(j);
    }

    return d(p-1);
}
// deBoor's algorithm for spline derivatives.
// [[Rcpp::export]]
double deBoorDerivative(double x,int k,Eigen::VectorXd t,Eigen::VectorXd beta,int p) {
	// x: evaluation point
	// k: knot index 
	// t: knots
	// beta: ORIGINAL coefficients. will be transformed within the algorithm
	// p: spline ORDER (cubic = 4)
	//
	// outputs f(x) via deBoor's algorithm

	Eigen::VectorXd d(p-1);
	for (int j=0;j<p-1;j++)
		d(j) = (p-1) * (beta(j+k-(p-1)+1) - beta(j+k-(p-1))) / (t(j+k+1) - t(j+k-(p-1)+1));
	
	double a =0.;
	for (int r=1;r<p-1;r++)
		for (int j=(p-2);j>r-1;j--) {
			a = (x-t(j+k-(p-2))) / (t(j+1+k-r) - t(j+k-(p-2)));
			d(j) = (1. - a)*d(j-1) + a*d(j);
		}

	return d(p-2);
}
// Bounded Newton.
// [[Rcpp::export]]
double reflect(double xt,double lb,double ub) {
  double wt = std::fmod(abs(xt - lb),(2*(ub-lb)));
  double out = std::min(wt,2*(ub - lb) - wt) + lb;
  return out;
}
// [[Rcpp::export]]
double bounded_newton(Function g,Function gprime,Eigen::VectorXd bounds,double xt,double eps,int maxitr) {
  int itr=1;
  NumericVector gt = g(xt), gpt;
  while((itr < maxitr) && (abs(gt(0)) > eps)) {
    gt = g(xt);
    gpt = gprime(xt);
    xt -= gt(0) / gpt(0);
    itr++;
  }
  return xt;
}

// Function to get gamma from beta
// [[Rcpp::export]]
Eigen::VectorXd get_gamma(Eigen::VectorXd beta) {
  int d = beta.size();
  Eigen::VectorXd gamma(d);
  gamma(0) = beta(0);
  for (int i=1;i<d;i++)
    gamma(i) = gamma(i-1) - exp(beta(i));

  return gamma;
}

// Get BMD without function passing, full in C++
// [[Rcpp::export]]
double Ux_cpp(double x,Eigen::VectorXd beta,Eigen::VectorXd knots,int k,double fx0,double sigmaest,double A) {
  double fxb = deBoor(x,k,knots,beta,4);
  return (fx0 - fxb)/sigmaest - A;
}
// [[Rcpp::export]]
double Uxd_cpp(double x,Eigen::VectorXd beta,Eigen::VectorXd knots,int k,double sigmaest) {
  // Pass in DERIVATIVE knot sequence
  // UPDATE: no, pass in the original knot sequence, same signature as the spline function
  //double fxbp = deBoor(x,k,knots,beta,3); // OLD
  double fxbp = deBoor(x,k,knots,beta,4); // NEW: same signature as the function
  return -fxbp/sigmaest;
}

// Variance of U(x)
// [[Rcpp::export]]
double Vx_cpp(double x,Eigen::MatrixXd V,Eigen::VectorXd knots,Eigen::VectorXd bx0,double sigmaest) {
  // V: covariance matrix of gamma
  Eigen::MatrixXd L(V.llt().matrixL());
  L.transposeInPlace(); // Upper triangular
  int m = V.cols(), p = bx0.size();

  // spline
  Eigen::VectorXd bxdiff(m);
  bxdiff.setZero();
  bxdiff.segment(0,p) = bx0 - Bsplinevec(x,knots,4);

  // system solve
  Eigen::VectorXd outvec = L * bxdiff;

  return outvec.squaredNorm() / (sigmaest*sigmaest);
}

// Derivative of variance of U(x)
// [[Rcpp::export]]
double Vxd_cpp(double x,Eigen::MatrixXd V,Eigen::VectorXd knots,Eigen::VectorXd bx0,double sigmaest) {
  int m = V.cols(), p = bx0.size();
  Eigen::VectorXd bxdiff(m);
  bxdiff.setZero();
  bxdiff.segment(0,p) = bx0 - Bsplinevec(x,knots,4);
  Eigen::VectorXd bxD = BsplinevecD(x,knots,4);

  bxdiff = V * bxdiff;
  return -2.*bxdiff.dot(bxD) / (sigmaest*sigmaest);
}

// [[Rcpp::export]]
double Psix_cpp(double x,Eigen::VectorXd beta,Eigen::MatrixXd V,Eigen::VectorXd knots,int k,double fx0,Eigen::VectorXd bx0,double sigmaest,double A) {
  double Ux = Ux_cpp(x,beta,knots,k,fx0,sigmaest,A);
  double Vx = Vx_cpp(x,V,knots,bx0,sigmaest);
  double chisq = 3.841459;
  return (Ux*Ux) - Vx*chisq;
}

// [[Rcpp::export]]
double Psixd_cpp(double x,Eigen::VectorXd beta,Eigen::VectorXd betadiff,Eigen::MatrixXd V,Eigen::VectorXd knots,int k,double fx0,Eigen::VectorXd bx0,double sigmaest,double A) {
  double Ux = Ux_cpp(x,beta,knots,k,fx0,sigmaest,A);
  double Uxd = Uxd_cpp(x,betadiff,knots,k,sigmaest);
  double Vxd = Vxd_cpp(x,V,knots,bx0,sigmaest);
  double chisq = 3.841459;

  return 2*Ux*Uxd - Vxd*chisq;
}


// [[Rcpp::export]]
double get_bmd_cpp(Eigen::VectorXd beta,Eigen::VectorXd knots,Eigen::VectorXd bounds,double x0,double sigmaest,double A,double eps,int maxitr) {
  // Setup initial quantities
  Eigen::VectorXd gamma = get_gamma(beta);
  // std::cout << "gamma = " << gamma << std::endl << std::endl;
  double fx0 = deBoor(x0,knotindex(x0,knots),knots,gamma,4);
  // Differenced coefficients
  int d = gamma.size(), p=4;
  /** Eigen::VectorXd gammadiff(d-1);
  gammadiff(0) = 0;
  for (int i=1;i<d-1;i++)
    gammadiff(i) = (p-1)*(gamma(i) - gamma(i-1)) / (knots(i+p-1) - knots(i));
  **/
  // Newton
  int itr=1,k=0;
  double xt = (bounds(0)+bounds(1))/2.;
  double gt=1.+eps; // Make sure it's bigger than eps to start
  double gpt = 1;
  while((itr < maxitr) && (abs(gt) > eps)) {
    k = knotindex(xt,knots);
    gt = Ux_cpp(xt,gamma,knots,k,fx0,sigmaest,A);
    //gpt = Uxd_cpp(xt,gammadiff,knots,k,sigmaest);
    gpt = Uxd_cpp(xt,gamma,knots,k,sigmaest);
    // std::cout << "itr " << itr << " xt " << xt << " gt " << gt << " gpt " << gpt << " k " << k << std::endl;
    xt -= gt / gpt;
    xt = reflect(xt,bounds(0),bounds(1));
    itr++;
  }

  return xt;
}

// [[Rcpp::export]]
double get_score_cpp(Eigen::VectorXd beta,Eigen::MatrixXd V,Eigen::VectorXd knots,Eigen::VectorXd bounds,double x0,double sigmaest,double A,double eps,int maxitr) {
  // Setup initial quantities
  Eigen::VectorXd gamma = get_gamma(beta);
  // std::cout << "gamma = " << gamma << std::endl << std::endl;
  double fx0 = deBoor(x0,knotindex(x0,knots),knots,gamma,4);
  Eigen::VectorXd bx0 = Bsplinevec(x0,knots,4);
  // Differenced coefficients
  int d = gamma.size(), p=4;
  Eigen::VectorXd gammadiff(d-1);
  gammadiff(0) = 0;
  for (int i=1;i<d-1;i++)
    gammadiff(i) = (p-1)*(gamma(i) - gamma(i-1)) / (knots(i+p-1) - knots(i));

  // Newton
  int itr=1,k=0;
  double xt = (bounds(0)+bounds(1))/2.;
  double gt=1.+eps; // Make sure it's bigger than eps to start
  double gpt = 1;
  while((itr < maxitr) && (abs(gt) > eps)) {
    k = knotindex(xt,knots);
    gt = Psix_cpp(xt,gamma,V,knots,k,fx0,bx0,sigmaest,A);
    gpt = Psixd_cpp(xt,gamma,gammadiff,V,knots,k,fx0,bx0,sigmaest,A);
    // std::cout << "itr " << itr << " xt " << xt << " gt " << gt << " gpt " << gpt << " k " << k << std::endl;
    xt -= gt / gpt;
    xt = reflect(xt,bounds(0),bounds(1));
    itr++;
  }

  return xt;
}
