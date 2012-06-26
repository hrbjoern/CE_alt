/** \class Utilities::Statistics
 *
 * Statistics functions to properly calculate the significance
 * of an excess
 *
 * This is a class only containing static member functions. If you add
 * a member function to this class, make sure that it is static!!
 *
 * \author Jim Hinton, Mathieu de Naurois
 *
 * $Revision: 1.44 $
 * $Name: hap-12-03 $
 */
// $Source: /cvs/utilities/src/Statistics.C,v $
// $Author: deil $
// $Date: 2011/06/29 13:12:21 $

#include "Statistics.hh"

#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#ifdef USE_ROOT
#include <TMath.h>
#include <TGraph.h>
#endif

#include "StringTools.hh"

#define DEBUG 0 // change to 1 to enable debugging here. Then, use the
		// DEBUG_OUT macro defined in "debugging.hh" instead
		// of "std::cout" to write out debug info.
#include "debugging.hh"

/**
 * Returns Alpha (exposure_on/exposure_off) with checks on it blowing up
 */
double Utilities::Statistics::Alpha( double exposure_on, double exposure_off ){

  if (exposure_on < GetMinExposure()) {
    throw InvalidResult("ON Exposure is too low");
    return 0;
  }

  if (exposure_off < GetMinExposure()) {
    throw InvalidResult("OFF Exposure is too low for numerical precision");
    return 0;
  }

  return exposure_on/exposure_off;

}


/**
 * A safe function to calculate the significance (including checks on
 * input values). In general, this function calls the LiMa formula,
 * but it first checks for the validity of the results.
 *
 * \throws Utilities::Statistics::InvalidResult if the inputs are out of the
 * bounds where the Li&Ma formula is valid.
 */
double Utilities::Statistics::SafeSignificance(int non,int noff,double alpha) {

  
  // we don't want to return values that are crazy (e.g. when alpha is
  // very small or n1/n2 are small) so do some sanity checks:

  if (non < GetMinCounts() || noff < GetMinCounts()) {
    throw InvalidResult( "LiMa: N_on or N_off is out of bounds" );
  }
  
  if (alpha < GetMinAlpha() || alpha > GetMaxAlpha()) {
    throw InvalidResult( "LiMa: Alpha is out of bounds" );
  }

  return LiMa(non,noff,alpha); 

}

/** Significance of excess counts.
 *
 * According to Li&Ma statistics (ApJ \b 272, 317-324, 1983), 
 * the significance of an excess \f$ \mathrm{non} - \alpha\times \mathrm{noff} \f$
 * is given by the formula 
 *
 * \f[ S = \sqrt {-2\ln \lambda } \f]
 *
 * where \f$\lambda\f$ is the \b likelihood-ratio defined as
 *
\f[ \lambda = \frac{L(X|E_0,\hat T_c)}{L(X|\hat E,\hat T)}  =  \
\left[ \frac {\alpha}{1+\alpha} \left( \frac {\mathrm{non} + \mathrm{noff}}{\mathrm{non}}\right)\right]^\mathrm{non} \times \
\left[ \frac {1}{1+\alpha} \left( \frac {\mathrm{non} + \mathrm{noff}}{\mathrm{noff}}\right)\right]^\mathrm{noff} \f]
 */
double Utilities::Statistics::LiMa(int non, int noff, double alpha) {
    
  double excess = non - noff*alpha;
  double n1,n2,sign,a;

  if(alpha == 0) return sqrt(1.0*non);

  if (excess > 0) {
    n1 = non;
    n2 = noff;
    sign = 1.0;
    a = alpha;
  } else {
    n2 = non;
    n1 = noff;
    sign = -1.0;
    a = 1.0/alpha;
  }

  // A development rather gives this values for
  // non == 0 or noff == 0
  if(n2 == 0) 
    return sign * sqrt(2 * n1 * log((1+a)/a));
  if(n1 == 0)  
    return sign * sqrt(2 * n2 * log(1+a));

  // the standard Li & Ma formula:
  
  double nt = n1+n2;
  double pa = 1 + a;

  double t1 = n1*log((pa/a)*(n1/nt));
  double t2 = n2*log(pa*(n2/nt));
  double sig = sqrt(2)*sqrt(fabs(t1+t2));

  return sign*sig;  
}

/** the log of the probability of observing events non and noff for a given
 * signal and an unknown background.
 * The best estimate of the background is first derived.
 * Implemented by Mathieu
 *
 * \note These "Li & Ma errors" were implemented by Mathieu. They are not described
 * in the Li & Ma paper, which only deals with significance. But they do have
 * excellent coverage and work better than the simple error formula
 * sqrt(on + alpha * alpha * off) from Equation (4) in the Li & Ma paper.
 */ 
double Utilities::Statistics::LiMa_LogLikelihood(int non, int noff, double alpha,double signal)
{
  // This is a little trick to avoid problems for negative signal
  // (note that this method is often called with signal=excess).
  // In that case the Poisson signal + b0 can become negative
  // and log(signal + b0) isn't defined, i.e. the likelihood
  // cannot be computed.
  // To avoid this problem, the on and off regions are switched
  // (I'm being a little vague here since I don't 100% understand
  // what statistical model actually results from this switch.)
  double n1,n2,sign,a;
  if (signal > 0) {
    n1 = non;
    n2 = noff;
    sign = 1.0;
    a = alpha;
  } else {
    n2 = non;
    n1 = noff;
    sign = -1.0;
    a = 1.0/alpha;
    signal = -a * signal;
  }

  // Maximum likelihood background estimate b0
  double t = (1 + a) * signal - a * (n1 + n2);
  double delta = t * t + 4 * a * (1 + a) * signal * n2;
  double b0 = (sqrt(delta) - t)/ (2 + 2 * a);
  // b0 should always be >= 0.
  // Just to avoid possible problems with rounding errors:
  if(b0 <0) b0 = 0; 

  // Compute Poisson log likelihood r for two measurements: on / off
  double r = 0;
  if(n1)
    r += n1 * log(signal + b0);
  if(n2)
    r += + n2 * log(b0 / a);
  r += - (signal + b0)  -  b0/a;
  return r;
}

/** 
 * Statistical uncertainties on excess.
 *
 * Searches the average value of the signal
 * for which the log-likelihood ratio is nsig away.
 *
 * Implemented by Mathieu.
 *
 * \see Note in LiMa_LogLikelihood
 * \note The required step in log likelihood for "nsig sigma errors" is
 * offset = 0.5 * nsig * nsig
 * http://wwwasdoc.web.cern.ch/wwwasdoc/minuit/node30.html
 */
double Utilities::Statistics::LiMa_dExcess_Up(int non, int noff, double alpha,double nsig)
{
  // Optimal value of signal
  double excess0 = non - noff*alpha;
  double L,L0;
  L0 = LiMa_LogLikelihood(non,noff,alpha,excess0);

  // Search for signal such as -2 * log H > 1 * sig
  // i.e. find a signal > dExcess_Up so that we can
  // use an interval bisection algorithm to find dExcess_Up
  double offset = nsig * nsig *  0.5; 
  double signal = excess0 + 1;
  int count = 0;
  L = LiMa_LogLikelihood(non,noff,alpha,signal);
  while((L0-L) < offset && finite(L) && (++count)<100) {
    signal = (signal - excess0) * 2 + excess0;
    L = LiMa_LogLikelihood(non,noff,alpha,signal);
  }
  // Dichotomic search (i.e. interval bisection algorithm)
  double low = excess0;
  double up = signal;
  signal = (up + low)/2;
  L = LiMa_LogLikelihood(non,noff,alpha,signal);
  count = 0;
  while(fabs(L0-L-offset) > 1e-4 && finite(L) && (++count)<100) {
    if(L0-L > offset) up = signal;
    else low = signal;
    signal = (up + low)/2;
    L = LiMa_LogLikelihood(non,noff,alpha,signal);
  }
  return fabs(signal - excess0);
}

/** 
 * Statistical uncertainties on excess.
 *
 * Searches the average value of the signal
 * for which the log-likelihood ratio is nsig away.
 *
 * Implemented by Mathieu.
 *
 * \see Notes in LiMa_LogLikelihood and LiMa_dExcess_Up
 */
double Utilities::Statistics::LiMa_dExcess_Down(int non, int noff, double alpha,double nsig)
{
  // Optimal value of signal
  double excess0 = non - noff*alpha;
  double L,L0;
  L0 = LiMa_LogLikelihood(non,noff,alpha,excess0);
  
  // Search for signal such as -2 * log H > 1
  double offset = nsig * nsig * 0.5;
  double signal = excess0 - 1;
  int count = 0;
  L = LiMa_LogLikelihood(non,noff,alpha,signal);
  while((L0-L) < offset && finite(L) && (++count)< 100) {
    signal = (signal - excess0) * 2 + excess0;
     L = LiMa_LogLikelihood(non,noff,alpha,signal);
  }
  // Dichotomic search
  double up = excess0;
  double low = signal;
  signal = (up + low)/2;
  L = LiMa_LogLikelihood(non,noff,alpha,signal);
  count = 0;
  while(fabs(L0-L-offset) > 1e-4 && finite(L) && (++count)<100) {
    if(L0-L > offset) low = signal;
    else up = signal;
    signal = (up + low)/2;
    L = LiMa_LogLikelihood(non,noff,alpha,signal);   
  }
  //   if(!finite(L)) {
  //       std::cerr << "Infinite likelihood : non=" << non
  // 		<< ", noff=" << noff << ", signal=" << signal
  // 		<< std::endl;
  //   }
  return fabs(signal - excess0);
}

/** 
 * Arithmetic mean of the LiMa up and down errors on excess.
 *
 * LiMa_dExcess = (LiMa_dExcess_Up + LiMa_dExcess_Dn) / 2
 */
double Utilities::Statistics::LiMa_dExcess(int non, int noff, double alpha,double nsig)
{
  double err_up = LiMa_dExcess_Up(non, noff, alpha, nsig);
  double err_down = LiMa_dExcess_Down(non, noff, alpha, nsig);
  double err = (err_up + err_down) / 2;
  return err;
}

/**
 * Calculates the required excess to reach given significance, given the background alpha*noff
 *
 * Uses Newton method as implemented in GSL to find the root of the LiMa formula
 */
double Utilities::Statistics::LiMa_requiredExcess(double significance, int noff, double alpha) {

  // sanity checks
  if( significance <= 0.1 ) throw InvalidResult( "Requested significance too small!" );
  if( noff < GetMinCounts() ) throw InvalidResult( "N_off is too small!" );
  if( noff > 100000 ) throw InvalidResult( "N_off is too large!" );
  if( alpha < 1.e-5 ) throw InvalidResult( "alpha is too small!" );
  if( alpha > 100. ) throw InvalidResult( "alpha is too large!" );

  // prevent GSL from aborting the program in case of problems
  static bool error_handler_off = false;
  if( !error_handler_off ){
    gsl_set_error_handler_off();
    error_handler_off = true;
  }

  // GSL stuff
  int status;
  int iter = 0;
  int max_iter = 100;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double x0;
  gsl_function_fdf FDF;
  Utilities::lima_params params;
  params.noff = double(noff);
  params.alpha = alpha;
  params.target_sigma = significance;
     
  // start value for N_on, estimated from simple Significance=sqrt(excess)
  double x = significance*significance+alpha*noff;

  DEBUG_OUT << "parameters: " << significance << " " << noff << " " << alpha << std::endl;
  DEBUG_OUT << "start value: " << x << std::endl;


  FDF.f = &Utilities::lima;
  FDF.df = &Utilities::lima_deriv;
  FDF.fdf = &Utilities::lima_fdf;
  FDF.params = &params;
     
  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x);
     
  do {
    iter++;
    status = gsl_root_fdfsolver_iterate (s);
    if( status == GSL_EBADFUNC ){
      DEBUG_OUT << "bad function value in iter " << iter << " (x=" << x
		<< ") for parameters " << significance << " " << noff << " " << alpha << std::endl;
      throw InvalidResult("GSL_EBADFUNC");
    }
    if( status == GSL_EZERODIV ){
      WARN_OUT << "vanishing derivative in iter " << iter << " (x=" << x
	<< ") for parameters " << significance << " " << noff << " " << alpha << std::endl;
      throw InvalidResult("GSL_EZERODIV");
    }
    x0 = x;
    x = gsl_root_fdfsolver_root (s);
    status = gsl_root_test_delta (x, x0, 0, 1e-3);
  } while (status == GSL_CONTINUE && iter < max_iter);
     
  gsl_root_fdfsolver_free (s);
  if( status != GSL_SUCCESS ){
    WARN_OUT << "No convergence for parameters " << significance << " " << noff << " " << alpha << std::endl;
  }

  double non = x;
  double excess = non - alpha*noff;

  DEBUG_OUT << "excess: " << excess << " Non: " << non
	    << " => Significance: " << Utilities::Statistics::LiMa( int(non), noff, alpha ) << std::endl;

  return excess;
}


// 3 helper functions for LiMa_requiredExcess
double Utilities::lima (double x, void *params)
{
  struct Utilities::lima_params *p 
    = (struct Utilities::lima_params *) params;
     
  double noff = p->noff;
  double alpha = p->alpha;
  double ts = p->target_sigma;
     
  double S = sqrt(2*(x*log((1+alpha)/alpha*x/(x+noff))+noff*log((1+alpha)*noff/(x+noff))));

  return S - ts;
}
     
double Utilities::lima_deriv (double x, void *params)
{
  struct Utilities::lima_params *p 
    = (struct Utilities::lima_params *) params;
     
  double noff = p->noff;
  double alpha = p->alpha;
     
  double S = sqrt(2*(x*log((1+alpha)/alpha*x/(x+noff))+noff*log((1+alpha)*noff/(x+noff))));

  return 1.0/S * log((1.+alpha)/alpha*x/(x+noff));
}
     
void Utilities::lima_fdf (double x, void *params, 
			  double *y, double *dy)
{
  struct Utilities::lima_params *p 
    = (struct Utilities::lima_params *) params;
     
  double noff = p->noff;
  double alpha = p->alpha;
  double ts = p->target_sigma;
     
  double S = sqrt(2*(x*log((1+alpha)/alpha*x/(x+noff))+noff*log((1+alpha)*noff/(x+noff))));

  *y = S - ts;
  *dy = 1.0/S * log((1.+alpha)/alpha*x/(x+noff));
}


/** Calculates lower limits - 84% (1 sigma) for poisson event numbers.           
 * Derived from V. H. Regener Phys. Rev. 84, p161, 1951   
 * Good to better than 0.4% for any n.                    
 */
double Utilities::Statistics::ErrorDown(int n)
{
  double lbound[5] = { 0.0, 0.17, 0.71, 1.37, 2.09 };
  double d,l;
  
  if (n > 50) {
    l=(double)(n)-sqrt((double)(n));
  } else if (n < 5) {
    l=lbound[n];
  } else {
    d=0.125-0.065*log10((double)(n));
    l=(double)(n)-sqrt((double)(n))+d;
  }
  
  return l-(double)(n);
}

/** Calculates lower limits - 95% (1.645 sigma) for poisson event numbers.           
 * Derived from V. H. Regener Phys. Rev. 84, p161, 1951   
 * Good to better than 0.4% for any n.                    
 */
double  Utilities::Statistics::ErrorDown95Percent(int n)
{

  double lbound[5] = { 0.0, 0.05, 0.35, 0.82, 1.37 };
  double d,l;
  
  if (n > 50) {
    d=0.59;
    l=(double)(n)-1.645*sqrt((double)(n))+d;
  } else if (n < 5) {
    l=lbound[n];
  } else {
    d=0.7-0.064*log10((double)(n));
    l=(double)(n)-1.645*sqrt((double)(n))+d;
  }

  return l-(double)(n);
  
} 

/** Calculates upper limits - 84% (1 sigma) and  
 * 95% (1.645 sigma) for poisson event numbers.           
 * Derived from V. H. Regener Phys. Rev. 84, p161, 1951   
 * Good to better than 0.4% for any n.                    
 */
double Utilities::Statistics::ErrorUp(int n) 
{

  double ubound[5] = {1.84, 3.3, 4.64, 5.92, 7.16 };
  double u,d;

  if (n > 50) {
    d=1.0;
    u=(double)(n)+sqrt((double)(n))+d;
  } else if (n < 5) {
    u=ubound[n];
  } else {
    d=1.23-0.120*log10((double)(n));
    u=(double)(n)+sqrt((double)(n))+d;
  }
  
  return u-(double)(n);
} 

/** Calculates upper limits - 95% (1.645 sigma) for poisson event numbers.           
 * Derived from V. H. Regener Phys. Rev. 84, p161, 1951   
 * Good to better than 0.4% for any n.                    
 */
double Utilities::Statistics::ErrorUp95Percent(int n) 
{

  double ubound[5] = { 3.0, 4.74, 6.3, 7.75, 9.15 };
  double d,u;

  if (n > 50) {
    d=1.57;
    u=(double)(n)+1.645*sqrt((double)(n))+d;
  } else if (n < 5) {
    u=ubound[n];
  } else {
    d=1.98-0.205*log10((double)(n));
    u=(double)(n)+1.645*sqrt((double)(n))+d;
  }
  
  return u-(double)(n);
} 

/** Returns the value of UP that provides 1, 2 and 3 sigma coverage
 * for chi^2 = chi_min^2+UP for npar parameters in a fit.
 * See $HESSROOT/utilities/test/print_table.C for references. (Ulli Schwanke)
 */
double Utilities::Statistics::GetUP(int npar, int sigmas){
  static double errcont[10][3] = {
    { 1.00021,4.00009,9.0004 }, //npar=1 confidence levels = 68.27% 95.45% 99.73% 
    { 2.29568,6.18057,11.8294 }, //npar=2 confidence levels = 68.27% 95.45% 99.73% 
    { 3.5265,8.02517,14.1558 }, //npar=3 confidence levels = 68.27% 95.45% 99.73% 
    { 4.72031,9.71603,16.2514 }, //npar=4 confidence levels = 68.27% 95.45% 99.73% 
    { 5.88703,11.3134,18.2047 }, //npar=5 confidence levels = 68.27% 95.45% 99.73% 
    { 7.03869,12.8489,20.0621 }, //npar=6 confidence levels = 68.27% 95.45% 99.73% 
    { 8.17575,14.3374,21.8463 }, //npar=7 confidence levels = 68.27% 95.45% 99.73% 
    { 9.30405,15.789,23.5741 }, //npar=8 confidence levels = 68.27% 95.45% 99.73% 
    { 10.4241,17.2116,25.2574 }, //npar=9 confidence levels = 68.27% 95.45% 99.73% 
    { 11.5366,18.611,26.9003 } //npar=10 confidence levels = 68.27% 95.45% 99.73% 
  };
  if(npar<1){
    std::cout << "Utilities::Statistics::GetUP: number of parameters must be > 0 !" << std::endl;
    return -1;
  }
  if(npar>10){
    std::cout << "Utilities::Statistics::GetUP: run $HESSROOT/utilities/test/print_table.C to fix me!" << std::endl;
    return -1;
  }
  if(sigmas<1 || sigmas>3){
    std::cout << "Utilities::Statistics::GetUP: number of sigmas must be 1, 2, 3 not " << sigmas << std::endl;
    return -1;
  }
  return errcont[npar-1][sigmas-1];
}

#ifdef USE_ROOT

double Utilities::Statistics::GetHeleneConf(double excess, double err_excess, double conf_level1)
{
  // WRB

  double  conf_level2 = 1. - conf_level1/100.;

  if(err_excess <= 0.) return 0.;

  if(excess >= 0.){
    double zeta  = excess/err_excess;
    double value = zeta/sqrt(2.);
    double integral = (1.+erf(value))/2.;
    double integral2 = 1. -  conf_level2 * integral;
    
    double value_old = value;
    /*  The 1st Loop is for Speed  & 2nd For Precision*/
    double value_new = value_old + .01;

    if(integral > integral2) value_new = 0.;

    integral = (1.+erf(value_new))/2.;

    while(integral < integral2){
      value_old = value_new;
      value_new = value_new + .01;
      integral = (1.+erf(value_new))/2.;
    }  
    
    value_new = value_old + .0000001;
    
    integral = (1.+erf(value_new))/2.;
    
    while(integral < integral2){
      value_old = value_new;
      value_new = value_new + .0000001;
      integral = (1.+erf(value_new))/2.;
    }  
    
    value_new = value_new * sqrt(2.);  /* Due to Wierdness with erf */
    
    double conf_limit = (value_new + zeta) * err_excess;
    return conf_limit;
  }  /* if */
  else if(excess < 0.){
    
    double zeta  = -excess/err_excess;
    double value = zeta/sqrt(2.);
    double integral = 1 - (1.+erf(value))/2.;
    double integral2 = 1. - conf_level2 * integral;
    double value_old = value;

    /*  The 1st Loop is for Speed & 2nd for Precision*/

    double value_new = value_old + .01;
    
    integral = (1.+erf(value_new))/2.;

    while(integral < integral2){
      value_old = value_new;
      value_new = value_new + .01;
      integral = (1.+erf(value_new))/2.;
    }  
    
    value_new = value_old + .0000001;

    integral = (1.+erf(value_new))/2.;

    while(integral < integral2){
      value_old = value_new;
      value_new = value_new + .0000001;
      integral = (1.+erf(value_new))/2.;
    }  

    value_new = value_new* sqrt(2.);   /* Due to Wierdness with erf */
    
    double conf_limit = (value_new-zeta)*err_excess;
    return conf_limit;
  }  /* else if */
  return 0;
}

/** Gaussian approximation for upper limit
 * Requires Non >> 1
 * Basic stats but see - Helene NIM 212 (1983) 319
 */
double Utilities::Statistics::GetHeleneConf2(double excess, double error, double conf_level)
{
  // Gaussian approximation for upper limit
  // Requires Non >> 1
  // Basic stats but see - Helene NIM 212 (1983) 319

  if (error <= 0.) return -1.0;  
  
  double alpha = 1. - conf_level/100.; // alpha = I1 / I2
  
  // Guess ul - find I1 and iterate to give correct alpha  
  double ul = excess; // Only valid for alpha<0.5
  if (alpha>0.5) ul = excess-error*3.0;

  double sigma = excess/error;
  double sigma_ul = (ul-excess)/error; // Start at zero 
  double I1 = 1.0 - TMath::Freq(sigma_ul);
  double I2 = 1.0 - TMath::Freq(-sigma);
  double alpha_guess = I1/I2;

  while (alpha_guess > alpha) {
    double step = 0.03;
    if (alpha_guess-alpha < 0.05) step = 0.0003;
    sigma_ul += step;
    I1 = 1.0 - TMath::Freq(sigma_ul);
    alpha_guess = I1/I2;
  }
  
  return sigma_ul*error + excess;  
}

double Utilities::Statistics::GetHeleneConf(int rawon, int rawoff, double norm, double conf_level1)
{

  double excess = double(rawon) - norm * double(rawoff);
  double err_excess = sqrt(double(rawon) + norm * norm * double(rawoff));
  return Utilities::Statistics::GetHeleneConf(excess, err_excess, conf_level1);
}

double Utilities::Statistics::GetChisqProb(double chisq_fit, int dof_fit){

  if(dof_fit < 1) return -999.;
  if(dof_fit == 1 && chisq_fit > 30) return 0.;
  if(dof_fit == 2 && chisq_fit > 40) return 0.;
  if(dof_fit > 3 && (chisq_fit/dof_fit) > 10.) return 0.;
  if(dof_fit > 10 && (chisq_fit/dof_fit) > 5.) return 0.;
  if(dof_fit > 50 && (chisq_fit/dof_fit) > 3.) return 0.;

  /* This function evaluates the probability of getting a larger chisq */
  /* Uses the total Chisq and the Number of degrees of freedom */

  double temp_index = double(dof_fit) / 2.;
  double temp_prob = 1. / (exp(lgamma(temp_index)) * pow(2.,temp_index));

  double temp_step = 0.00001;
  if(dof_fit == 3) temp_step = 0.00005;
  if(dof_fit == 2) temp_step = 0.000001;
  if(dof_fit == 1) temp_step = 0.0000001;
  double temp_val_int = chisq_fit - temp_step/2.;

  double temp_func = 0.;
  while(temp_val_int > 0){
    temp_func += pow(temp_val_int,(temp_index-1.)) * exp((-temp_val_int/2.)) * temp_step;
    temp_val_int -= temp_step;
  }
  temp_prob = temp_prob * temp_func;

  double value = 1.-temp_prob;

  return value;
}

double Utilities::Statistics::GetSigmaProb(double sigma){

  /* This function takes a sigma and returns the corresponding probability */

  double value = 0.;
  double temp_step = 1.e-7;
  double temp_val = 20. - temp_step/2.;

  if(sigma < 0.){
    value = 0.5;
    temp_step = 1.e-7;
    temp_val = - temp_step / 2.;
  }

  if(sigma >= 20.){
    std::cout << "\nThis could take a while" << std::endl;
    temp_val = 40. - temp_step/2.;
    if(sigma >= 40.) temp_val = 100. - temp_step/2.;
    if(sigma >= 100) return 0;
  }

  if(sigma == 0.) return 0.5;
  else if(sigma > 0.){
    while(temp_val > sigma){
      value += (1. / sqrt(2*TMath::Pi())) * exp(-temp_val*temp_val/2.) * temp_step;
      temp_val -= temp_step; 
    }
  }
  else if(sigma < 0.){
    while(temp_val < fabs(sigma)){
      value -= (1. / sqrt(2*TMath::Pi())) * exp(-temp_val*temp_val/2.) * temp_step;
      temp_val += temp_step;
    }
  }

  if(sigma > 0.) return value;
  else return 1.-value;
}
    
void Utilities::Statistics::GetBinomErrorBars(int n, int m, int conf_level, int num_steps, 
					      double &f, double &f_low, double &f_up)
{
  // for given n successes out of m trials, each with a probability p
  // find f which maximises the probability distribution  
  // P(p, n, m) = p^n * (1-p)^(m-n) * (m! / (n!*(m-n)!))
  // and calculate the upper and lower error bar for a certain confidence level
  // (choice of confidence levels: 68%, 95% and 99%)

  int max_step = 0;
  double p = 0.;
  double prob = 0.;
  double ln_prob_level = 0.;
  double max_prob = 0.;
  double low_prob = 0.;
  double up_prob = 0.;

  //confidence level
  if(conf_level != 68 && conf_level != 95 && conf_level != 99) {
    std::cout << "Choose confidence level 68, 95 or 99!" << std::endl;
    return;
  }
  if (conf_level == 68) ln_prob_level = 1./2.;
  if (conf_level == 95) ln_prob_level = 4./2.;//95.45%
  if (conf_level == 99) ln_prob_level = 6.63/2.;

  std::vector<double> xp;
  std::vector<double> yp;

  //maximum probability
  for (int i=0; i<num_steps; i++) {
    p=(1.+i)/num_steps;
    prob = TMath::Power(p,n)*TMath::Power(1-p,m-n)*TMath::Binomial(m,n);
    if (prob>1e-10) {
      xp.push_back(p);
      yp.push_back(-log(prob));
    }
    if(prob > max_prob) {
      max_prob = prob;
      max_step = i;
      f = p;
    }
  }
  std::cout << "Calculated max prob to: " << max_prob << 
    "   -ln(max prob): " << -log(max_prob) << "   at f=" << f  << std::endl;

  //level for lower/upper error bars
  ln_prob_level += -log(max_prob);
  std::cout << "-ln(level prob) =  " << ln_prob_level << 
               "   for confidence level: " << conf_level << std::endl; 

  TGraph *graph = new TGraph(xp.size(),&xp[0],&yp[0] );
  graph->Draw("AL*");

  f_low = 0;
  //lower error bar
  for (int i=0; i<max_step; i++) {
    p=(1.+i)/num_steps;
    prob = TMath::Power(p,n)*TMath::Power(1-p,m-n)*TMath::Binomial(m,n); 
    if(-log(prob) > ln_prob_level) {
      low_prob = prob;
      f_low = p;
    }
  }
  std::cout << "Calculated low prob to: " << low_prob <<            
    "   -ln(low prob): " << -log(low_prob) << "   at f_low=" << f_low << std::endl;

  //upper error bar
  for (int i=max_step; i<num_steps; i++) {
    p=(1.+i)/num_steps;
    prob = TMath::Power(p,n)*TMath::Power(1-p,m-n)*TMath::Binomial(m,n); 
    if(-log(prob) < ln_prob_level) {
      up_prob = prob;
      f_up = p;
    }
  }
  std::cout << "Calculated up prob to: " << up_prob << 
    "   -ln(up prob): " << -log(up_prob) << "   at f_up=" << f_up <<  std::endl;

  return;
}

// //////////////////////////////////////////////////////////////////////////////
void Utilities::Swap (double& aA, double& aB) 
{ double temp = aA; aA = aB; aB = temp; }
// //////////////////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////////////////////////////
double Utilities::Sign (const double aA, const double aB) 
{return ( aB >= 0.0 ? fabs (aA) : - fabs(aA)); }
// //////////////////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////////////////////////////
double Utilities::select (const int aK, std::vector<double>& aV)
// //////////////////////////////////////////////////////////////////////////////
// NR
// Returns the kth smallest value in the array arr[1..n]. The input array will be rearranged
// to have this value in location arr[k], with all smaller elements moved to arr[1..k-1] (in
// arbitrary order) and all larger elements in arr[k+1..n] (also in arbitrary order).
{
  int i,ir,j,l,mid;
  double a;
  l = 0;
  ir = aV.size () - 1;
  
  if (aK < 0 || aK > ir)
    {
      std::cout << "Select ERROR requested bin out of range: " << aK << " , ir " << ir << std::endl;
      return 1.0;
    }
  
  for (;;) 
    {

      //      for (int p = 0; p <= ir; p++) std::cout << aV[p] << std::endl;

      // Active partition contains 1 or 2 elements.
      if (ir <= l+1) 
	{
	  //  Case of 2 elements.
	  if (ir == l+1 && aV[ir] < aV[l]) 
	    Swap (aV[l], aV[ir]);
	    if (! std::isfinite (aV[aK]) )
	      std::cout << "Select pb: " << aK << " " << aV[aK] << std::endl;
	  //std::cout << "select: " << aK << " " << aV[aK] << endl;
	  return aV[aK];
	} 
      else 
	{
	  // Choose median of left, center, and right elements as partitioning element a. Also
	  // rearrange so that arr[l] = arr[l+1], arr[ir] = arr[l+1].
	  mid = (l+ir) >> 1;
	  Swap (aV[mid], aV[l+1]);
	  if (aV[l] > aV[ir]) 
	    Swap (aV[l], aV[ir]);
	  if (aV[l+1] > aV[ir]) 
	    Swap (aV[l+1], aV[ir]);
	  if (aV[l] > aV[l+1]) 
	    Swap (aV[l], aV[l+1]);
	  
	  //Initialize pointers for partitioning.
	  i = l+1; 
	  j = ir;
	  // Partitioning element.
	  a = aV[l+1];
	  //Beginning of innermost loop.
	  for (;;) 
	    { 
	      // Scan up to find element > a.
	      do i++; 
	      while (aV[i] < a && i <= ir);
	      //  Scan down to find element < a.
	      do j--;
	      while (aV[j] > a && j >= 0);
	      // Pointers crossed. Partitioning complete.
	      if (j < i) 
		break;
	      Swap (aV[i], aV[j]);
	    } 
	  //End of innermost loop.
	  // Insert partitioning element.
	  aV[l+1] = aV[j]; 
	  aV[j] = a;
	  // Keep active the partition that contains the kth element.
	  if (j >= aK) 
	    ir = j-1; 
	  if (j <= aK) 
	    l = i;
	}
    }
  
  if (! std::isfinite (aV[aK]) )
   std::cout << "Select pb: " << aK << " " << ir << aK << " " << aV[aK] << std::endl;

  return aV[aK];
}


// //////////////////////////////////////////////////////////////////////////////
double Utilities::rofunc (const std::vector<double>& aX, 
			  const std::vector<double>& aY, 
			  double& aA, const double aB, 
			  double& aAbDev, const int aNZ)
// //////////////////////////////////////////////////////////////////////////////
// NR 
// Evaluates sum _i=1 ^N x_i sign(y_i - a - b*x_i) for a given value of b. 
{
  double EPS = 1.0e-7;
  double d;
  double sum = 0.0;
  std::vector<double> vec (aNZ);
  int ndata = aX.size ();
  int j, k;
  int usedpoints = 0;
  
  k = 0;
  for (j = 0; j < aNZ; j++) vec[j] = 0.;
  for (j = 0; j < ndata; j++) 
    if (aX[j] != 0.0)
      {
	vec[k] = aY[j] - aB*aX[j];
	k++;
      }
  
  // If ndata odd
  if (aNZ & 1) 
    {
      j = (aNZ+1)>>1;
      j--;
      aA = select (j, vec);
    }
  else 
    // ndata even
    {
      j = aNZ >> 1;
      j--;
      aA = 0.5* (select (j, vec) + select (j+1, vec) );
    }
  
  aAbDev = 0.0;
  for (j=0; j < ndata; j++) 
    // Adaptation for selected range
    if (aX[j] != 0.0)
      {
	d = aY[j] - (aB * aX[j] + aA);
	aAbDev += d*d;
	
	usedpoints++;
	//aAbDev += fabs (d);
	if (aY[j] != 0.0) 
	  d /= fabs (aY[j]);
	if (fabs (d) > EPS) 
	  sum += (d >= 0.0 ? aX[j] : -aX[j]);
      }

  if (usedpoints)
    aAbDev /= (double) usedpoints;
  
  if (! std::isfinite (sum) )
    std::cout << "Rofunc pb: " << sum << std::endl;

  return sum;
}

// //////////////////////////////////////////////////////////////////////////////
double Utilities::medfit (const std::vector<double>& aX, const std::vector<double>& aY, 
			  double& aA, double& aB)
// //////////////////////////////////////////////////////////////////////////////
// NR
// Fits y = a + bx by the criterion of least absolute deviations. The arrays aX[1..ndata] and
// aY[1..ndata] are the input experimental points. The fitted parameters aA and aB are output,
// along with return value abdev, which is the mean absolute deviation (in y) of the experimental points 
// from the fitted line.
{
  int ndata = aX.size ();
  int zerNb = 0;
  int ndatanz = ndata;
  
  double aa, bb, b1, b2, del, f, f1, f2, sigb, temp;
  double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0, chisq = 0.0, abDev = 0.0;
  
  //  As a first guess for a and b, we will find the least-squares fitting line.
  for (int j=0; j < ndata; j++) 
    if (aX[j] == 0.0)
      zerNb++;
    else
      {
	sx += aX[j];
	sy += aY[j];
	sxy += aX[j] * aY[j];
	sxx += aX[j] * aX[j];
      }
  ndatanz = ndata - zerNb;
  
  aA = 0.;
  aB = 0.;
  if (ndatanz == 0)
    return 0.; 

  del = ndatanz * sxx - sx * sx;
  //  Least-squares solutions.
  aa = (sxx*sy-sx*sxy) / del;
  bb = (ndatanz*sxy-sx*sy) / del;
  for (int j=0; j < ndata; j++)
    if (aX[j] != 0.0)
      {
	temp = aY[j] - (aa + bb * aX[j]);
	chisq += temp*temp;
      }
  //  std::cout << "Medfit: chi2 " << chisq << " / " << del << std::endl;
  //  The standard deviation will give some idea of how big an iteration step to take.
  sigb = sqrt (chisq / del);
  b1 = bb;
  f1 = rofunc (aX, aY, aa, b1, abDev, ndatanz);
  
  if (sigb > 0.0) 
    {
      // Guess bracket as 3-s away, in the downhill direction known from f1.
      b2 = bb + Sign (3.0*sigb, f1);
      f2 = rofunc (aX, aY, aa, b2, abDev, ndatanz);
      if (b2 == b1) 
	{
	  aA = aa;
	  aB = bb;
	  return sqrt (abDev); // / (double) ndatanz;
	}
      // Bracketing.
      //      std::cout << "Medfit: Bracketing" << std::endl;
      while (f1*f2 > 0.0) 
	{
	  bb = b2 + 1.6 * (b2-b1);
	  b1 = b2;
	  f1 = f2;
	  b2 = bb;
	  f2 = rofunc (aX, aY, aa, b2, abDev, ndatanz);
	}
      //  Refine until error a negligible number of standard deviations.
      //      std::cout << "Medfit: Refine" << std::endl;
      sigb = 0.01*sigb;
      while ( fabs (b2-b1) > sigb) 
	{
	  // Bisection.
	  bb = b1 + 0.5 * (b2 - b1); 
	  if (bb == b1 || bb == b2) 
	    break;
	  f =  rofunc (aX, aY, aa, bb, abDev, ndatanz);
	  if (f*f1 >= 0.0) 
	    {
	      f1 = f;
	      b1 = bb;
	    } 
	  else 
	    {
	      f2 = f;
	      b2 = bb;
	    }
	}
    }
  
  aA = aa;
  aB = bb;
  
  //  std::cout << "Medfit: " << sqrt(abDev) << " / " << ndatanz << std::endl;

  return sqrt(abDev); // / (double) ndatanz;
}


// //////////////////////////////////////////////////////////////////////////////
double Utilities::medfit (const std::vector<double>& aX, const std::vector<double>& aY, 
			  double& aA)
// //////////////////////////////////////////////////////////////////////////////
// NR
// Fits y = a by the criterion of least absolute deviations. The arrays aX[1..ndata] and
// aY[1..ndata] are the input experimental points. The fitted parameter aA is output,
// along with return value abdev, which is the mean absolute deviation (in y) of the experimental points 
// from the fitted line.
{
  int ndata = aX.size ();
  int zerNb = 0;
  int ndatanz = ndata;
  
  double aa, bb, /*b1, b2,*/ del, /*f,*/ f1,/* f2,*/ sigb, temp;
  double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0, chisq = 0.0, abDev = 0.0;
  
  //  As a first guess for a and b, we will find the least-squares fitting line.
  for (int j=0; j < ndata; j++) 
    if (aX[j] == 0.0)
      zerNb++;
    else
      {
	sx += aX[j];
	sy += aY[j];
	sxy += aX[j] * aY[j];
	sxx += aX[j] * aX[j];
      }
  ndatanz = ndata - zerNb;
  
  aA = 0.;
  if (ndatanz == 0)
    return 0.; 

  del = ndatanz * sxx - sx * sx;
  //  Least-squares solutions.
  aa = (sxx*sy-sx*sxy) / del;
  bb = (ndatanz*sxy-sx*sy) / del;
  for (int j=0; j < ndata; j++)
    if (aX[j] != 0.0)
      {
	temp = aY[j] - aa;
	chisq += temp*temp;
      }
  //  std::cout << "Medfit: chi2 " << chisq << " / " << del << std::endl;
  //  The standard deviation will give some idea of how big an iteration step to take.
  sigb = sqrt (chisq / del);
  //  b1 = bb;
  f1 = rofunc (aX, aY, aa, 0., abDev, ndatanz);
  
//   if (sigb > 0.0) 
//     {
//       // Guess bracket as 3-s away, in the downhill direction known from f1.
//       b2 = bb + Sign (3.0*sigb, f1);
//       f2 = rofunc (aX, aY, aa, 0., abDev, ndatanz);
//       if (b2 == b1) 
// 	{
// 	  aA = aa;
// 	  return sqrt (abDev); // / (double) ndatanz;
// 	}
//       // Bracketing.
//       //      std::cout << "Medfit: Bracketing" << std::endl;
//       while (f1*f2 > 0.0) 
// 	{
// 	  bb = b2 + 1.6 * (b2-b1);
// 	  b1 = b2;
// 	  f1 = f2;
// 	  b2 = bb;
// 	  f2 = rofunc (aX, aY, aa, 0., abDev, ndatanz);
// 	}
//       //  Refine until error a negligible number of standard deviations.
//       //      std::cout << "Medfit: Refine" << std::endl;
//       sigb = 0.01*sigb;
//       while ( fabs (b2-b1) > sigb) 
// 	{
// 	  // Bisection.
// 	  bb = b1 + 0.5 * (b2 - b1); 
// 	  if (bb == b1 || bb == b2) 
// 	    break;
// 	  f =  rofunc (aX, aY, aa, 0., abDev, ndatanz);
// 	  if (f*f1 >= 0.0) 
// 	    {
// 	      f1 = f;
// 	      b1 = bb;
// 	    } 
// 	  else 
// 	    {
// 	      f2 = f;
// 	      b2 = bb;
// 	    }
// 	}
//     }
  
  aA = aa;
  
  //  std::cout << "Medfit: " << sqrt(abDev) << " / " << ndatanz << std::endl;

  return sqrt(abDev); // / (double) ndatanz;
}

#endif










