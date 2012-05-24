#ifndef _ROOT_TBOUNDED_GAUSSIAN_
#define _ROOT_TBOUNDED_GAUSSIAN_

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TRolke.h>

namespace Utilities {

  enum { STAT_VERBOSE=1, STAT_DEBUG=2,STAT_QUICK=4, STAT_HILL=0,STAT_LIMA=8,STAT_ROLKE=16 };

  class TBoundedGaussian : public TObject {
  public:
    TBoundedGaussian(Double_t CL=0.6827,Double_t syst=0.0,int flags=STAT_QUICK | STAT_ROLKE);
    ~TBoundedGaussian();
    //Get confidence interval for event counts: Non-norm*Noff
    //Confidence interval boundaries are returned in min & max 
    int GetConfidenceInterval(int Non,int Noff,Double_t norm,Double_t& min,Double_t& max);
    //Get confidence interval for quantity x with error sigmax
    //Confidence interval boundaries are returned in min & max 
    int GetConfidenceInterval(Double_t x,Double_t sigmax,Double_t& min,Double_t& max);
    Double_t GetCL() const { return fCL; }    
    Double_t GetSystErr() const { return fSyst; }
    Double_t GetXLimit() const { return fXlimit; }
    Double_t GetErrLimit() const { return fErrLimit; }
  private:
    void WriteToFile(const char* name,bool update=true);  // development only
    void PlotPDF(const char* fn);                         // development only
  private:
    Double_t fCL;               //confidence level 0.0 < CL < 1.0
    Double_t fSyst;             //systematic error 0.0 <= err < 0.5 
    Double_t fMu;               //true value. Used temporarily.
    Double_t fSigma;            //width of Gaussian PDF. Used temporarily.
    Int_t fFlag;                //flags 
    TH1D* fLower;               //left,
    TH1D* fUpper;               //right limit of confidence belts
    Double_t fXlimit;           //transition between UL and measurement   
    Double_t fErrLimit;         //error for x>>1
    TRolke* fRolke;             //pointer to TRolke object for  calculation of confidence intervals
                                //when Non and Noff are small (if STAT_TROLKE is set (default))
    
    int GetConfidenceInterval(Double_t x,Double_t& min,Double_t& max);    
    void Calc(Double_t& x1,Double_t& x2);
    Double_t Interpolate(TH1D* h,int a,Double_t xx,bool IsX=true);
    Double_t Interpolate(Double_t x[],Double_t y[],Double_t xx);
    Double_t GetOrdering(Double_t p,Double_t x);
    Double_t GetLogOrdering(Double_t logp,Double_t x);
    Double_t GetProb(Double_t x);
    Double_t GetLogProb(Double_t x);
    Double_t FindZero(Double_t (*ff)(Double_t,Double_t), Double_t CL,Double_t min,Double_t max);
    Double_t TranslateX(TH1D*,Double_t x,Double_t sign=1.0);
    Double_t DoubleGauss(Double_t x,Double_t mu,Double_t eps);
    Double_t IntegrateDeps(Double_t x);
    Double_t GetEffectiveSigma(Double_t mu);
    void Fit(TH1D* h);
    void UpdateOnFile(TFile& f,TObject* obj,bool update=true);
    TString MakeName(const char *s);
    double IntegrateGaus(Double_t x1,Double_t x2,Double_t mu=0.0,Double_t sigma=1.0);
    double IntegrateProb(Double_t x1,Double_t x2);
    void SetMu(Double_t tru);
    ClassDef(Utilities::TBoundedGaussian,1)
  };
}

#endif
 
 
 
