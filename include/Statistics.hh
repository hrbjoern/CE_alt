#ifndef ROOT_Utils_Statistics
#define ROOT_Utils_Statistics

#ifdef USE_ROOT
#include <TROOT.h>
#endif

#include <string>
#include <stdexcept>
#include <vector>

namespace Utilities {

  class Statistics {
  public:
    virtual ~Statistics() {}
    static double LiMa(int non, int noff, double alpha = 1.0);
    static double Alpha( double exposure_on, double exposure_off );
    static double Significance(int non, int noff, double alpha = 1.0){return LiMa(non,noff,alpha);}
    static double SafeSignificance(int non, int noff, double alpha = 1.0);
    static double LiMa_LogLikelihood(int non, int noff, double alpha,double signal);
    static double LiMa_dExcess_Up(int non, int noff, double alpha = 1.0,double nsig = 1.0);
    static double LiMa_dExcess_Down(int non, int noff, double alpha = 1.0,double nsig = 1.0);
    static double LiMa_dExcess(int non, int noff, double alpha = 1.0,double nsig = 1.0);
    static double LiMa_requiredExcess(double significance, int noff, double alpha = 1.0);

    // Real poisson error bars for small statistics
    static double ErrorDown(int n);
    static double ErrorDown95Percent(int n);
    static double ErrorUp(int n);
    static double ErrorUp95Percent(int n);

    static double GetChisqProb(double chisq_fit, int dof_fit);
    static double GetSigmaProb(double sigma);
    static double GetUP(int npar,int sigmas);

#ifdef USE_ROOT    
    static double GetHeleneConf(double excess, double err_excess, double conf_level1 = 99.0);
    static double GetHeleneConf2(double excess, double error, double conf_level = 99.0);
    static double GetHeleneConf(int rawon, int rawoff, double norm, double conf_level1 = 99.0);
    static void GetBinomErrorBars(int n, int m, int conf_level, int num_steps, double &max_prob, 
			   double &lower_prob,  double &upper_prob );
#endif

    static double GetMinExposure(){ return 1.0e-8; }
    static double GetMinAlpha(){ return 1.0e-6; }
    static double GetMaxAlpha(){ return 1.0e6; }
    static int    GetMinCounts(){ return 5; }


    /**
     * This exception is thrown by statistical methods if the result
     * is not valid (due to inputs out of bounds, for example). For
     * instance, if the LiMa formula is called with too low statistics
     * where it is not a valid measure of significance, this should
     * be thrown and caught by calling functions.
     */
    class InvalidResult : public std::runtime_error {
    public:
      InvalidResult(std::string why) : runtime_error(why.c_str()){ ; }
    };
    
#ifdef USE_ROOT
    ClassDef(Utilities::Statistics,0);
#endif
  };
  
  double medfit (const std::vector<double>& aX, const std::vector<double>& aY, 
			  double& aA, double& aB);

  double medfit (const std::vector<double>& aX, const std::vector<double>& aY, 
		 double& aA);

  void Swap (double& aA, double& aB);
    
  double Sign (const double aA, const double aB);
    
  double select (const int aK, std::vector<double>& aV);
    
  double rofunc (const std::vector<double>& aX, 
		 const std::vector<double>& aY, 
		 double& aA, const double aB, 
		 double& aAbDev, const int aNZ);


  // needed for LiMa_requiredExcess
  struct lima_params
  {
    double noff, alpha, target_sigma;
  };
     
  double lima (double x, void *params);
  double lima_deriv (double x, void *params);
  void lima_fdf (double x, void *params, 
		 double *y, double *dy);


}

#endif
