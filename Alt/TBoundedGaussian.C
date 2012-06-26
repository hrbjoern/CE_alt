/** \class Utilities::TBoundedGaussian
 *
 * \deprecated After discussions by Ulli Schanke, Christoph Deil and
 * Vincent Marandon we decided to use TRolke for any number of on or
 * off counts intead of this class.
 * Results are practically identical and TRolke has the advantage that
 * it is a published, well-studied, fast method one can reference in papers.
 * This class was left unchanged and is simply no longer used in the HESS
 * analysis software.
 *
 * \note You can no longer state in your paper "Limits were computed using
 * the Feldman & Cousins method." This was never really true anyways because
 * you should have been more explicitly stating that the model from Section
 * IV.b (Gaussian with a boundary at the origin) was used, plus if you had
 * less than 10 on or off counts the Rolke method was used before.
 * Now you can state: "Limits were computed using the profile likelihood
 * method" and refer to TRolke or Rolke's paper.
 * If you want to be very specific: Model 4 in TRolke or Rolke's paper was
 * used: Poisson background, known efficiency.
 *
 * \brief Confidence Intervals according to Feldman & Cousins (the one-sided
 * Gaussian model described in section IV.b, not the more widely known Poisson 
 * model from Section IV.a) with systematic errors according to Hill. 
 * For low-statistics situations (number of on or off counts smaller than 10), 
 * the profile likelihood method (as implemented in Root's TRolke) is used. 
 *
 * \author Ulli Schwanke
 *
 * See the following internal H.E.S.S. note for detailed explanations
 * http://www-eep.physik.hu-berlin.de/hess/protected/limits.pdf
 * and consult http://lanl.arxiv.org/abs/physics/0403059 
 * and http://root.cern.ch/root/html/TRolke.html 
 * for more information on TRolke and the profile likelihood.
 *
 * Here is an example how to use the class.
 \code
 #include <utilities/TBoundedGaussian.hh>
 
 Utilities::TBoundedGaussian fc(0.99); //99% confidence level
 double upper(0),lower(0);
 double flux = ...;
 double sigma_flux = ...;
 std::cout << "Measured flux is: " << flux << std::endl;
 std::cout << "Error of flux is: " << sigma_flux << std::endl;
 fc.GetConfidenceInterval(flux,sigma_flux,lower,upper);
 if( lower<=0 )
   std::cout << fc.GetCL()*100.0 << "% CL upper limit on flux: " << upper << std::endl;
 else {
   std::cout << "Confidence interval for true flux: " << lower << " .. " << upper << std::endl;
   std::cout << "Flux is (at " << fc.GetCL()*100.0 << " confidence level): " << std::endl;
   std::cout << flux << " -" << flux-lower << " +" << upper-flux << std::endl;
 }
 \endcode
 *
 */
#include "TBoundedGaussian.hh"
#include "Statistics.hh"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <TMath.h>
#include <TF1.h>
#include <TROOT.h>



/** Calculate for given CL where the transition from UL to measurement occurs.
 * (Why is there no ErfcInverse in our current ROOT version??)
 * Formula: \int_{-\infty}^x exp(-0.5*t^2) dt = CL
 */
static Double_t Func1(Double_t x,Double_t CL){
  return TMath::Freq(x)-CL;
}

/** Calculate for given CL the +/- n sigma interval.
 * (Why is there no ErfInverse in our current ROOT version??)
 * Formula: \int_{-s}^{+s} exp(-0.5*t^2) dt = CL
 */
static Double_t Func2(Double_t x,Double_t CL){
  const Double_t sqrt2 = sqrt(2.0);
  return TMath::Erf(x/sqrt2)-CL;
}

/*! Constructor of a Utilities::TBoundedGaussian object.
 *
 * The first argument is the desired confidence level (default: 0.6827 i.e. 68.27%),
 * the second arguments is the fractional systematic error (default: 0),
 * the third argument can be used to enable debugging.
 */
Utilities::TBoundedGaussian::TBoundedGaussian(Double_t cl,Double_t syst,int flag) : 
  fCL(cl),fSyst(syst),fMu(0),fSigma(0),fFlag(flag),fLower(0),fUpper(0),fXlimit(0),fErrLimit(-1),fRolke(0) {
  if( cl<=0.0 || cl>=1.0 ){
    throw 
      std::range_error("TBoundedGaussian: confidence level must be in (0,1)");
  }
  if( syst>0.5 ){
    throw 
      std::range_error("TBoundedGaussian: syst. error must be <= 50%");
  }
  Double_t alpha = 1.0-fCL;
  TString upnam = MakeName("upper");
  TString lonam = MakeName("lower");   
  fLower = dynamic_cast<TH1D*>(gDirectory->Get(lonam));
  fUpper = dynamic_cast<TH1D*>(gDirectory->Get(upnam));
  if( fFlag & STAT_VERBOSE ){
    if(fLower) std::cout << "Found " << lonam << std::endl;
    if(fUpper) std::cout << "Found " << upnam << std::endl;
  }
  bool run(fLower==0 || fUpper==0);
  const Double_t muMax = 10.0; //calculate confidence intervals for mu = 0.0..10.0 in
  const int N = 1000;          //steps of 0.01
  if( fUpper==0) fUpper = new TH1D(upnam,upnam,N,0,muMax);
  if( fLower==0) fLower = new TH1D(lonam,lonam,N,0,muMax);  
  if(run){
    int last = -1;
    for(int i=1;i<=N;i++){
      Double_t tru = fUpper->GetBinCenter(i);
      SetMu(tru);
      Double_t x1(0),x2(0);    
      Calc(x1,x2);
      fLower->SetBinContent(i,x1);
      fUpper->SetBinContent(i,x2);
      fLower->SetBinError(i,1e-4);
      fUpper->SetBinError(i,1e-4);
      if( i>1 && (fFlag & STAT_DEBUG) ){
	if( x1<fLower->GetBinContent(i-1) ){
	  std::cerr << "lo: Unexpected behaviour: " << fLower->GetBinCenter(i-1) << " : " 
		    << fLower->GetBinContent(i-1) << " " << GetEffectiveSigma(fLower->GetBinCenter(i-1)) << std::endl;
	  std::cerr << "lo:                     : " << fLower->GetBinCenter(i-0) << " : " 
		    << fLower->GetBinContent(i-0) << " " << GetEffectiveSigma(fLower->GetBinCenter(i-0)) << std::endl;
	}
	if( x2<fUpper->GetBinContent(i-1) ){
	  std::cerr << "up: Unexpected behaviour: " << fUpper->GetBinCenter(i-1) << " : " << fUpper->GetBinContent(i-1) << std::endl;
	  std::cerr << "up:                       " << fUpper->GetBinCenter(i) << " : " << fUpper->GetBinContent(i) << std::endl;
	}
      }	
      if( fFlag & STAT_VERBOSE ){
	int percent = int(100.0*Double_t(i)/Double_t(N)+0.5);
	if( percent != last ){
	  last = percent;
	  std::cout << last << " %" << std::endl;
	}
      }
      if(tru<1 && (fFlag & STAT_DEBUG)) 
	printf("%5.3f: %5.3f..%5.3f: %.8f\n",tru,x1,x2,100.0*IntegrateProb(x1,x2));
    }
    //For mu>6.0, use a linear interpolation of the tables. 
    //Fit between mu=5.0 and mu=6.0
    Fit(fLower);
    Fit(fUpper);
  }
  //Calculation transition from UL to measurement 
  if( fSyst>0.0 ){
    fXlimit = Interpolate(fUpper,2,0.0);
  } else {
    fXlimit = FindZero(Func1,cl,-100,100);
    fErrLimit = FindZero(Func2,cl,0,100);
  }
  if( fFlag & STAT_VERBOSE ){
    const int n = fLower->GetNbinsX();    
    std::cout << fUpper->GetName() << ' ' << fLower->GetName() << std::endl;
    std::cout << " CL,alpha = " << fCL << ", " << alpha << std::endl;
    std::cout << " xlimit   = " << fXlimit << std::endl;
    std::cout << " errlimit = " << fErrLimit << std::endl;
    std::cout << " syst     = " << fSyst << std::endl;
    std::cout << " Flag     = " << fFlag << std::endl;
    std::cout << " xmax     = " << fLower->GetBinContent(n) << ", " << fUpper->GetBinContent(n) << std::endl;
  }
}

/** Integrate normalized Gaussian between x1 and x2
 */
double Utilities::TBoundedGaussian::IntegrateGaus(Double_t x1,Double_t x2,Double_t mu,Double_t sigma){
  const Double_t sqrt2 = sqrt(2.0);
  return 1.0-TMath::Freq((x1-mu)/sigma)-0.5*TMath::Erfc((x2-mu)/sigma/sqrt2);
}

/** Integrate PDF between x1 and x2
 */
double Utilities::TBoundedGaussian::IntegrateProb(Double_t x1,Double_t x2){
  return IntegrateGaus(x1,x2,fMu,fSigma);  
}

void Utilities::TBoundedGaussian::SetMu(Double_t tru){
  Double_t newfSigma = GetEffectiveSigma(tru); 
  if( fSigma>0.1 && (fFlag & STAT_DEBUG) && (newfSigma<fSigma) ){
    printf("Warning in SetMu(mu,sigma): (%.7f,%.7f) -> (%.7f,%.7f)\n",fMu,fSigma,tru,newfSigma);
  }
  fMu = tru;
  fSigma = newfSigma;
  //printf("%.7f %.7f\n",fMu,fSigma);
}

/** Fit histogram with straight line between endpoint and endpoint-1.0. 
 * The resulting fit will be stored with the histogram.
 */
void Utilities::TBoundedGaussian::Fit(TH1D* h){
  const int n2 = h->GetNbinsX();
  Double_t mu2 = h->GetBinCenter(n2);
  Double_t mumax = mu2+0.5*h->GetBinWidth(n2);
  Double_t x2 = h->GetBinContent(n2);  
  const int n1 = h->FindBin(mu2-0.5);
  Double_t mu1 = h->GetBinCenter(n1);
  Double_t mumin = mu1-0.5*h->GetBinWidth(n1);
  Double_t x1 = h->GetBinContent(n1);  
  Double_t a = (x2*mu1-x1*mu2)/(mu1-mu2);
  Double_t b = (x2-x1)/(mu2-mu1);
  TString nam(h->GetName());
  nam += "_fit";
  TF1* f = new TF1(nam,"pol1",mumin,mumax);
  f->SetParameters(a,b);
  h->Fit(f,"Q0","N",mumin,mumax);
  //if( fFlag & STAT_VERBOSE ){    
  //std::cout << "Fit : " << h->GetName() << ": " << mumin << ".." << mumax << "   " << n1 << ' ' << n2 << std::endl;
  //std::cout << "init: " << a << ' ' << b << std::endl;
  //std::cout << "resu: " << f->GetParameter(0) << ' ' << f->GetParameter(1) << std::endl;
  //std::cout << "erro: " << f->GetParError(0) << ' ' << f->GetParError(1) << std::endl;
  //}
}

Utilities::TBoundedGaussian::~TBoundedGaussian() {
  if(fRolke) delete fRolke;
  if(fLower) delete fLower;
  if(fUpper) delete fUpper;
}

/** Find the null of a function.
 * Find the null of the passed function and return it. It is
 * assumed that there is indeed a null between min and max. 
 */
Double_t Utilities::TBoundedGaussian::FindZero(Double_t (*ff)(Double_t,Double_t), Double_t CL,Double_t min,Double_t max){
  Double_t fmin = ff(min,CL);
  Double_t fmax = ff(max,CL);
  Double_t lastx = 0.0;
  bool run(true);
  for(int steps=0;run;steps++){
    Double_t x = 0.5*(min+max);
    Double_t f = ff(x,CL);
    //std::cout << "x: " << min << ' ' << x << ' ' << max << std::endl;
    //std::cout << "f: " << fmin << ' ' << f << ' ' << fmax << std::endl;
    if( f*fmin<0 ){
      fmax = f;
      max = x;
    } else if( f*fmax<0 ){
      fmin = f;
      min = x;
    }
    if( steps>0 && fabs(lastx-x)<1e-5 ){
      run = false;
    } else if( steps>100){
      std::cout << "convergence error " << std::endl;
      run = false;
    }
    lastx = x;
  }
  return lastx;
}

/** Return the likelihood ordering.
 * Returns R(x,mu)/R(x,mubest)
 */
Double_t Utilities::TBoundedGaussian::GetOrdering(Double_t p,Double_t x){
  Double_t pbest = 1.0/sqrt(TMath::TwoPi());
  if( x<0 ){
    pbest *= exp(-0.5*x*x);
  }
  return p/pbest;
}

/** Return log() of the likelihood ordering.
 * Same as function above but log() for better numerical behaviour.
 */
Double_t Utilities::TBoundedGaussian::GetLogOrdering(Double_t logp,Double_t x){
  Double_t pbest = -0.5*log(TMath::TwoPi()); 
  if( x<0.0 ){
    pbest -= 0.5*x*x;
  }
  return logp-pbest;
}

/** Look up confidence interval.
 * The first two arguments are the measured value and its estimated stat. error.
 * The boundaries of the confidence interval are returned in the 3rd and 4th argument.
 * A zero lower bound means that the result is an upper limit.
 */
int Utilities::TBoundedGaussian::GetConfidenceInterval(Double_t x,Double_t sigmax,Double_t& min,Double_t& max){
  if( sigmax<=0.0 ){
    std::cerr << "TBoundedGaussian: sigmax must be positive, you passed " << sigmax << " !" << std::endl;
    return 0;
  }
  int ok = GetConfidenceInterval(x/sigmax,min,max);
  if(ok){
    min *= sigmax;
    max *= sigmax;
  }
  return ok;  
}

/** Convolution of two Gaussians.
 * This is needed for taking systematic errors into account.
 */
Double_t Utilities::TBoundedGaussian::DoubleGauss(Double_t x,Double_t mu,Double_t eps){
  const Double_t c = TMath::TwoPi();
  return exp(-0.5*pow(x-eps*mu,2))*exp(-0.5*pow((1.0-eps)/fSyst,2))/c/fSyst;
}

/** Convolution of two Gaussians and numeric integration. 
 * Beware of the very primitive integration algorithm.
 */
Double_t Utilities::TBoundedGaussian::IntegrateDeps(Double_t x){
  const Double_t LIMIT = 1e-10;
  if( fabs(x-fMu)<LIMIT ){
    //fMu==x case. Analytical solution is easy.
    return 1.0/sqrt(TMath::TwoPi()*(1+x*x*fSyst*fSyst));
  }
  if( fabs(fMu)<LIMIT ){
    //fMu==0 case. Analytical solution is easy.
    return exp(-0.5*x*x)/sqrt(TMath::TwoPi());
  }
  //calculate maximum point of PDF
  const Double_t s2 = fSyst*fSyst;
  const Double_t emax = (x*fMu*s2+1.0)/(fMu*fMu*s2+1.0);
  const Double_t fmax = DoubleGauss(x,fMu,emax);
  const Double_t estep = 0.1*fSyst;  
  const Double_t limit = 1e-8*fmax;
  Double_t integ = 0.0;  
  int n = 0;
  for(Double_t e=emax+0.5*estep;;e+=estep,n++){
    Double_t fe = DoubleGauss(x,fMu,e);
    integ += (fe*estep);
    if( fe <= limit ){
      //std::cout << "Stopping + after " << n << " steps: " << emax << " .. " << e << std::endl;
      break;
    }
  }
  n = 0;
  for(Double_t e=emax-0.5*estep;;e-=estep,n++){
    Double_t fe = DoubleGauss(x,fMu,e);
    integ += (fe*estep);
    if( fe <= limit ){
      //std::cout << "Stopping - after " << n << " steps: " << e << " .. " << emax << std::endl;
      break;
    }
  }  
  return integ;
}

/** Return the confidence interval for a measurement.
 * The first (second) argument is the number of on (off) counts. The third
 * argument is the on-off normalization. The lower (upper) limit of the confidence 
 * interval is returned in the 4th (5th) argument.
 * A zero lower bound means that the result is an upper limit. TRolke
 * is used when the number of on or off counts is smaller than 10.
 */
int Utilities::TBoundedGaussian::GetConfidenceInterval(int Non,int Noff,Double_t norm,Double_t& min,Double_t& max){
  if(norm<=0.0){
    std::cerr << "TBoundedGaussian: norm must be greater than 0, you passed " 
	      << norm << " !" << std::endl;
    return 0;
  }  
  if((fFlag & STAT_ROLKE) && (Non<10 || Noff<10)){
    if(fRolke==0){
      fRolke = new TRolke;
      fRolke->SetCL(fCL);
    }
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
    fRolke->SetPoissonBkgKnownEff(Non,Noff,1.0/norm,1);
    fRolke->GetLimits(min,max);
#else
    max = fRolke->CalculateInterval(Non,Noff,0,0.0,0.0,1.0,4,0.0,0.0,1.0/norm,0.0,0);
    min = fRolke->GetLowerLimit();
#endif
    return 1;
  } 
  if(Non<=0 || Noff<=0){
    std::cerr << "TBoundedGaussian: on and off counts must be greater than 0, you passed " 
              << Non << " and " << Noff << " !" << std::endl;
    return 0;
  }
  Double_t e = Double_t(Non)-norm*Double_t(Noff);
  Double_t s = sqrt(Double_t(Non)+norm*norm*Double_t(Noff));  
  return GetConfidenceInterval(e,s,min,max);
}

/** Convert an x-value into a mu-value.
 * The function uses a linear fit (for mu>6.0)
 * or interpolation with 2nd order polynomial (for mu<6.0).
 */
Double_t Utilities::TBoundedGaussian::TranslateX(TH1D* h,Double_t x,Double_t sign){
  if(h==0) {
    std::cout << "NULL histogram in Utilities::TBoundedGaussian::TranslateX" << std::endl;
    return -1;
  }
  const int n = h->GetNbinsX();
  if( x>h->GetBinContent(n) ){
    if( fErrLimit>0 ){
      return x+sign*fErrLimit;      
    }
    TString nam(h->GetName());
    nam += "_fit";
    TF1* f = h->GetFunction(nam);
    if( f==0 ){
      std::cerr << "Fit failed." << std::endl;
      return 0;
    }
    Double_t a = f->GetParameter(0);
    Double_t b = f->GetParameter(1);
    return (x-a)/b;
  }
  int besti = 1;
  Double_t val = fabs(h->GetBinContent(besti)-x);
  for(int i=2;i<=n;i++){
    Double_t newval = fabs(h->GetBinContent(i)-x);
    if( newval<val ){
      val = newval;
      besti = i;
    }
  }  
  //Double_t old = h->GetBinCenter(besti);
  int off = besti;
  if(off==1) off++;
  if(off==n) off--;
  Double_t neff = Interpolate(h,off,x,false);
  //std::cout << "old = " << old << " new = " << neff << " off=" << off << std::endl;
  return neff;
}

/** Interpolation with 2nd order polynomial.
 */
Double_t Utilities::TBoundedGaussian::Interpolate(TH1D* h,int off,Double_t xx,bool IsX){
  Double_t x[3] = {h->GetBinCenter(off-1),h->GetBinCenter(off),h->GetBinCenter(off+1)};
  Double_t y[3] = {h->GetBinContent(off-1),h->GetBinContent(off),h->GetBinContent(off+1)};
  if(IsX) return Interpolate(x,y,xx);  
  return Interpolate(y,x,xx);
}

/** Interpolation with 2nd order polynomial.
 */
Double_t Utilities::TBoundedGaussian::Interpolate(Double_t x[],Double_t y[],Double_t xx){      
  Double_t a0 = ((xx-x[1])*(xx-x[2])) / ((x[0]-x[1])*(x[0]-x[2]));
  Double_t a1 = ((xx-x[0])*(xx-x[2])) / ((x[1]-x[0])*(x[1]-x[2]));
  Double_t a2 = ((xx-x[0])*(xx-x[1])) / ((x[2]-x[0])*(x[2]-x[1]));
  return a0*y[0]+a1*y[1]+a2*y[2];
}

/* Return confidence interval.
 * Used internally. Returns confidence interval for unit resolution (i.e. sigma=1)
 */
int Utilities::TBoundedGaussian::GetConfidenceInterval(Double_t x,Double_t& min,Double_t& max){  
  if( x<=fXlimit ){
    //This is going to be an upper limit. 
    max = TranslateX(fLower,x,+1);
    min = 0.0; 
  } else {
    min = TranslateX(fUpper,x,-1);
    max = TranslateX(fLower,x,+1);    
  }
  return 1;
}

/** Evaluate PDF at x for fixed mu=tru.
 */
Double_t Utilities::TBoundedGaussian::GetProb(Double_t x){
  if( fSyst<=0.0 || (fFlag & STAT_QUICK) ){
    return TMath::Gaus(x,fMu,fSigma,kTRUE);
  }
  return IntegrateDeps(x);
}

/** Evaluate log(PDF) at x for fixed mu=tru.
 * Same as function above, but log() for better numerical behaviour.
 */
Double_t Utilities::TBoundedGaussian::GetLogProb(Double_t x){
  if( fSyst<=0.0 ){   
    return -0.5*(log(TMath::TwoPi())+pow(x-fMu,2));
  } else if( fFlag & STAT_QUICK ){
    return -0.5*(log(TMath::TwoPi()*fSigma*fSigma)+pow((x-fMu)/fSigma,2));
  }
  return log(IntegrateDeps(x));
}

/** Calculate confidence interval for given mu=tru.
 * Return boundaries of confidence interval in x1 and x2
 */
void Utilities::TBoundedGaussian::Calc(Double_t& x1,Double_t& x2){  
  const Double_t dx = 0.002;
  Double_t dxlo = dx;
  Double_t dxhi = dx;
  const Double_t fak = dx/GetLogProb(fMu);
  Double_t xhi= fMu+0.5*dxhi;
  Double_t xlo= fMu-0.5*dxlo;
  Double_t plo = GetLogProb(xlo);
  Double_t phi = GetLogProb(xhi);
  Double_t rlo = GetLogOrdering(plo,xlo);
  Double_t rhi = GetLogOrdering(phi,xhi);
  Double_t rmax = (rlo>rhi) ? rlo+1 : rhi+1;
  Double_t sum = 0.0;
  //std::cout << "start: " << fMu << ' ' << dxlo << ' ' << dxhi << std::endl;
  do {
    //std::cout << "lo: " << xlo << ' ' << rlo << ' ' << plo << ' ' << dxlo << std::endl;
    //std::cout << "hi: " << xhi << ' ' << rhi << ' ' << phi << ' ' << dxhi << std::endl;   
    if( rlo>rhi ){
      if( rlo>rmax && (fFlag & STAT_DEBUG)) 
	std::cout << "-Wrong assumptions or numerical problems! " << rlo << ' ' << rmax << ' ' << xlo << std::endl;      
      rmax = rlo;
      Double_t ds = exp(plo)*dxlo;
      sum += ds;     
      x1 = xlo-0.5*dxlo;
      dxlo = fak*plo;
      double poss = (fCL-sum)/2/exp(plo);
      if(dxlo>poss){
	//std::cout << "Truncate- " << dxlo << " to " << dx << " since poss is=" << poss << std::endl;
	dxlo=dx;
      }
      xlo = x1-0.5*dxlo;
      plo = GetLogProb(xlo);
      rlo = GetLogOrdering(plo,xlo);
      //std::cout << " - Add " << ds << " total=" << sum << std::endl;
      //std::cout << ds << ' ' << dxlo << std::endl;
    } else {
      if(rhi>rmax && (fFlag & STAT_DEBUG)) 
	std::cout << "+Wrong assumptions or numerical problems! " << rhi << ' ' << rmax << ' ' << rhi << std::endl;      
      rmax = rhi;
      Double_t ds = exp(phi)*dxhi;
      sum += ds;
      x2 = xhi+0.5*dxhi;
      dxhi = fak*phi;
      double poss = (fCL-sum)/2/exp(plo);
      if(dxhi>poss){
	//std::cout << "Truncate+ " << dxhi << " to " << dx << " since poss is=" << poss << std::endl;
	dxhi=dx;
      }
      xhi = x2+0.5*dxhi;
      phi = GetLogProb(xhi);
      rhi = GetLogOrdering(phi,xhi);
      //std::cout << " + Add " << ds << " total=" << sum << std::endl;
      //std::cout << ds << ' ' << dxhi << std::endl;
    }
    if( fSyst<=0.0 || (fFlag & STAT_QUICK )){
      sum = IntegrateProb(x1,x2);
    }
  } while( sum<fCL );
  //std::cout << "end  : " << tru << ' ' << dxlo << ' ' << dxhi << std::endl;
  //printf("%5.3f: %5.3f .. %5.3f\n",float(tru),float(x1),float(x2));
}

/** A utility function.
 * Just for development. 
 */
void Utilities::TBoundedGaussian::UpdateOnFile(TFile& f,TObject* obj,bool update){
  if(update){
    std::ostringstream b;
    b << obj->GetName() << ";1";
    f.Delete(b.str().c_str());
  }
  obj->Write();
  if( fFlag & STAT_VERBOSE ) std::cout << "Updated " << obj->GetName() << " on file " << f.GetName() << std::endl;
}

/** A utility function.
 * Just for development. 
 */
void Utilities::TBoundedGaussian::WriteToFile(const char* fname,bool update){
  const char* mode = (update) ? "UPDATE" : "RECREATE";
  TFile* f = new TFile(fname,mode);
  if(f){
    if( f->IsOpen() ){
      UpdateOnFile(*f,fLower,update);
      UpdateOnFile(*f,fUpper,update);
      f->Close();
    }
    delete f;
  }
}

/** Make object name.
 * Make object name depending on CL, systematic error and flags.
 */
TString Utilities::TBoundedGaussian::MakeName(const char *s){
  char tmp[128];
  int f = 0;
  if( fSyst>0 ) f = (fFlag & STAT_QUICK) ? 1 : 2;
  sprintf(tmp,"%s%08i_%04i_%i",s,int(1e7*fCL+0.5),int(fSyst*1e3+0.5),f);
  return TString(tmp);
}

/** When including syst. errors, the resulting PDF is again
 * a Gaussian but with sigma>1.0. The effective sigma is a function
 * of both the true value mu and the syst. error. We use a parameterization
 * sigma = 1.0 + mu*c1 + mu^2*c2 + mu^3*c3 + mu^4*c4 for a systematic
 * error of 1..50%
 */
Double_t Utilities::TBoundedGaussian::GetEffectiveSigma(Double_t mu){
  static Double_t coeff[50][5] = {
    { 1.0, 9.43526e-11, 5e-05, -3.32853e-11, -1.24322e-09 }, // 1 % systematic. error
    { 1.0, -5.3332e-09, 0.000200008, -3.57057e-09, -1.94606e-08 }, // 2 % systematic. error
    { 1.0, -6.08035e-08, 0.000450089, -3.99677e-08, -9.52039e-08 }, // 3 % systematic. error
    { 1.0, -3.28747e-07, 0.000800484, -2.17533e-07, -2.869e-07 }, // 4 % systematic. error
    { 1.0, -1.19156e-06, 0.00125176, -7.96834e-07, -6.59065e-07 }, // 5 % systematic. error
    { 1.0, -3.3458e-06, 0.00180498, -2.2669e-06, -1.26915e-06 }, // 6 % systematic. error
    { 1.0, -7.86223e-06, 0.00246178, -5.40793e-06, -2.15518e-06 }, // 7 % systematic. error
    { 1.0, -1.61876e-05, 0.00322445, -1.13268e-05, -3.32615e-06 }, // 8 % systematic. error
    { 1.0, -3.00962e-05, 0.00409586, -2.14625e-05, -4.7559e-06 }, // 9 % systematic. error
    { 1.0, -5.1589e-05, 0.00507935, -3.75595e-05, -6.38098e-06 }, // 10 % systematic. error
    { 1.0, -8.26986e-05, 0.00617853, -6.15908e-05, -8.10455e-06 }, // 11 % systematic. error
    { 1.0, -0.000125667, 0.0073974, -9.58147e-05, -9.7883e-06 }, // 12 % systematic. error
    { 1.0, -0.00018223, 0.00873955, -0.000142473, -1.12806e-05 }, // 13 % systematic. error
    { 1.0, -0.000253999, 0.0102085, -0.000203915, -1.24035e-05 }, // 14 % systematic. error
    { 1.0, -0.000342198, 0.0118073, -0.000282445, -1.29685e-05 }, // 15 % systematic. error
    { 1.0, -0.000447631, 0.0135387, -0.000380287, -1.27803e-05 }, // 16 % systematic. error
    { 1.0, -0.000570669, 0.0154046, -0.000499542, -1.16422e-05 }, // 17 % systematic. error
    { 1.0, -0.000711163, 0.0174065, -0.000642118, -9.36432e-06 }, // 18 % systematic. error
    { 1.0, -0.000868642, 0.0195452, -0.000809761, -5.76301e-06 }, // 19 % systematic. error
    { 1.0, -0.0010421, 0.0218209, -0.00100398, -6.67556e-07 }, // 20 % systematic. error
    { 1.0, -0.0012303, 0.0242331, -0.00122611, 6.0812e-06 }, // 21 % systematic. error
    { 1.0, -0.00143181, 0.0267812, -0.00147735, 1.46364e-05 }, // 22 % systematic. error
    { 1.0, -0.00164462, 0.0294635, -0.00175852, 2.51196e-05 }, // 23 % systematic. error
    { 1.0, -0.00186675, 0.0322783, -0.00207039, 3.7648e-05 }, // 24 % systematic. error
    { 1.0, -0.00209605, 0.0352233, -0.00241353, 5.23236e-05 }, // 25 % systematic. error
    { 1.0, -0.0023302, 0.0382959, -0.0027883, 6.92308e-05 }, // 26 % systematic. error
    { 1.0, -0.00256691, 0.0414935, -0.00319501, 8.84462e-05 }, // 27 % systematic. error
    { 1.0, -0.00280385, 0.0448131, -0.00363378, 0.000110031 }, // 28 % systematic. error
    { 1.0, -0.00303858, 0.0482514, -0.00410458, 0.000134033 }, // 29 % systematic. error
    { 1.0, -0.00326889, 0.0518054, -0.00460741, 0.000160497 }, // 30 % systematic. error
    { 1.0, -0.00349254, 0.0554718, -0.00514209, 0.000189455 }, // 31 % systematic. error
    { 1.0, -0.00370723, 0.0592469, -0.00570834, 0.000220924 }, // 32 % systematic. error
    { 1.0, -0.00391082, 0.0631276, -0.00630585, 0.000254919 }, // 33 % systematic. error
    { 1.0, -0.00410131, 0.0671102, -0.00693427, 0.000291448 }, // 34 % systematic. error
    { 1.0, -0.00427679, 0.0711916, -0.00759319, 0.000330512 }, // 35 % systematic. error
    { 1.0, -0.00443552, 0.0753684, -0.00828224, 0.000372115 }, // 36 % systematic. error
    { 1.0, -0.00457579, 0.0796375, -0.00900096, 0.000416249 }, // 37 % systematic. error
    { 1.0, -0.00469605, 0.0839957, -0.00974887, 0.000462907 }, // 38 % systematic. error
    { 1.0, -0.00479501, 0.0884401, -0.0105256, 0.00051208 }, // 39 % systematic. error
    { 1.0, -0.00487134, 0.0929679, -0.0113306, 0.000563763 }, // 40 % systematic. error
    { 1.0, -0.00492384, 0.0975761, -0.0121635, 0.000617937 }, // 41 % systematic. error
    { 1.0, -0.0049516, 0.102262, -0.0130238, 0.000674594 }, // 42 % systematic. error
    { 1.0, -0.00495359, 0.107024, -0.0139112, 0.000733723 }, // 43 % systematic. error
    { 1.0, -0.00492916, 0.111859, -0.0148252, 0.000795312 }, // 44 % systematic. error
    { 1.0, -0.00487755, 0.116764, -0.0157654, 0.000859348 }, // 45 % systematic. error
    { 1.0, -0.00479843, 0.121739, -0.0167316, 0.000925839 }, // 46 % systematic. error
    { 1.0, -0.00469114, 0.12678, -0.0177233, 0.000994763 }, // 47 % systematic. error
    { 1.0, -0.00455532, 0.131886, -0.0187402, 0.00106611 }, // 48 % systematic. error
    { 1.0, -0.00439065, 0.137056, -0.0197821, 0.00113988 }, // 49 % systematic. error
    { 1.0, -0.00419702, 0.142287, -0.0208486, 0.00121607 } // 50 % systematic. error
  };
  if( fSyst<=0.0 ) return 1.0;
  Double_t xx = 100.0*fSyst;
  Double_t x[2]={floor(xx),ceil(xx)};
  int k[2] = {int(x[0]),int(x[1])};
  if( k[0]==k[1] ){
    if( k[0]<1 ) return 1.0;
    if( k[0]>50 ) throw std::runtime_error("Talk to Ulli Schwanke about (513)");
    Double_t* p = &coeff[k[0]-1][0];
    return p[0]+mu*(p[1]+mu*(p[2]+mu*(p[3]+mu*p[4])));
  }
  Double_t y[2]={0,0};
  for(int i=0;i<2;i++){
    if(k[i]<1){
      y[i] = 1.0;
    } else if(k[i]>50){
      if( k[0]>50 ) throw std::runtime_error("Talk to Ulli Schwanke about (519)");
    } else {
      Double_t* p = &coeff[k[i]-1][0];
      y[i] = p[0]+mu*(p[1]+mu*(p[2]+mu*(p[3]+mu*p[4])));
    } 
  }
  //interpolate between different syst. errors
  return y[0]+(xx-x[0])/(x[1]-x[0])*(y[1]-y[0]);
}

/** A utility function.
 * Calculate PDF with systematic error
 * and fit it with a Gaussian. Store width.
 */
void Utilities::TBoundedGaussian::PlotPDF(const char* fn){
  const int N = 100;
  const int NN = 1000;
  char tmp[128];
  std::vector<TH1D*> vec;
  TH1D* htemp = 0;  
  for(Double_t sigma=0.010;sigma<=0.505;sigma+=0.010){
    fSyst = sigma;    
    sprintf(tmp,"h%03i",int(sigma*1000.0+0.5));
    std::cout << tmp << std::endl;
    TH1D* h = new TH1D(tmp,"",N,0,6);
    vec.push_back(h);
    sprintf(tmp,"j%03i",int(sigma*1000.0+0.5));
    std::cout << tmp << std::endl;
    TH1D* hh = new TH1D(tmp,"",N,0,6);
    vec.push_back(hh);
    for(int i=1;i<=N;i++){
      Double_t mu = h->GetBinCenter(i);
      SetMu(mu);
      const Double_t minx = mu-6;
      const Double_t maxx = mu+6;
      if(htemp) delete htemp;      
      htemp = new TH1D("htemp","",NN,minx,maxx);
      for(int k=1;k<=NN;k++){
	Double_t x = htemp->GetBinCenter(k);
	Double_t e = IntegrateDeps(x);
	htemp->SetBinContent(k,e);
	htemp->SetBinError(k,0.1);
      }
      htemp->Fit("gaus","Q");
      const int p = 2;
      Double_t e = htemp->GetFunction("gaus")->GetParameter(p);
      Double_t ee = htemp->GetFunction("gaus")->GetParError(p);
      h->SetBinContent(i,e);
      h->SetBinError(i,ee);
      hh->SetBinContent(i,GetEffectiveSigma(mu));
      hh->SetBinError(i,0);
    }
  }
  if(htemp) delete htemp;
  TFile f(fn,"RECREATE");
  for(std::vector<TH1D*>::iterator q=vec.begin();q!=vec.end();q++){
    (*q)->Write();    
  }    
  f.Close();
  std::cout << "Closed " << fn << std::endl;
}

ClassImp(Utilities::TBoundedGaussian)
