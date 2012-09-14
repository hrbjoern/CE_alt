#from PhotonSpectra import *
import PhotonSpectra
#from math import pi, exp, log, log10, sqrt
from numpy import pi, exp, log, log10, sqrt
from scipy.integrate import quad, quadrature
from scipy import interpolate
from scipy.stats import poisson, lognorm
import numpy as np
import fnmatch
from Limits import UpperLimit

class Object(object):
    """This class is supposed to hold the different \"objects\", i.e. 
    specific (publications of) DM searches towards a certain
    astrophysical object, and the formulae to calculate upper limits
    from their results."""
    
    def __init__(self, Name, Aeff, Tobs, Jbar, Non, Noff, alpha, spectrum):
        """ Initializes \"object\", including some specific object modifications."""
        #        print "\nHallo", Name
        self.Name = Name
        self.Eth = 30.
        if Name=="Segue1M":
            self.Eth=100.
        elif fnmatch.fnmatch(self.Name,'Scu*'):
            self.Eth=200. #Sculptor fudge factor. Carina?!

        self.Aeffdata = np.genfromtxt(Aeff)
        # Consider all energies in GeV, all Aeffs in cm**2:
        if Aeff=="Aeffs/VERITAS-Aeff_20deg.dat":
            pass
        elif Aeff=="Aeffs/MAGIC_Gaug_Aeff.dat":
            self.Aeffdata[:,1]*=100. # Factor 100 in MAGIC Aeff table
            self.Aeffdata[:,1]*=1e4  # m^2 to cm^2
            self.Eth = 100.
        else:
            self.Aeffdata[:,1]*=1e4 # m^2 to cm^2
            if Aeff!="Aeffs/VERITAS_Segue1_Aeff_TrueE.dat":
                self.Aeffdata[:,0]*=1000. # TeV to GeV
            #elif fnmatch.fnmatch(self.Name,'Scu*'):
                #self.Aeffdata[:,1]*=0.8 # Fudge efficiency factor

        self.Tobs = Tobs*3600.
        # Introducing Jbar PDF: 
        self.JbarError = 0.01 # Testing ...
        #self.JbarError = 0.5 # half a dec error for the moment
        self.Jbar = Jbar
        self.Non = Non
        self.Noff = Noff
        self.alpha = alpha
        #TRolke calculation:
        self.Nul = float(UpperLimit(self.Non, self.Noff, self.alpha))
        # Get spectrum function:
        # self.Spectrum = eval("PhotonSpectra."+spectrum) # works, but 
        self.Spectrum = getattr(PhotonSpectra, spectrum)  # seems better
        # To be used for sensitivity integration/Exclusion.py:
        self.SensiIntegralArray = np.zeros(1000000)
        # To be used for likelihood:
        self.logPois2 = poisson.logpmf(self.Noff, self.Noff)

                
    def printObject(self):
        """ Prints object, just for testing. """
        print
        print self.Name
        print self.Tobs, ", Tobs in sec"
        print self.Jbar, "GeV^2 cm^-5"
        print self.Nul, ", Nul at 95% CL"
        print self.Aeffdata[0,0], ", Emin(Aeff) in GeV"
        print

    def Aeff(self,Energy):
        """ Aeff function from the Aeffdata stored in each object."""
        if (self.Aeffdata[0,0] < Energy < self.Aeffdata[-1,0]):
            return float(interpolate.interp1d(self.Aeffdata[:,0],self.Aeffdata[:,1])(Energy))
        else:
            return 0.
        
    def Sensi(self, E, mchi): 
       return self.Spectrum(E, mchi)*self.Aeff(E)

    def SensiIntegral_scalar(self, mchi):
        return quad(self.Sensi, self.Eth, 1.01*mchi, args=(mchi),
                        limit=50,full_output=1)[0]
    
    def SensiIntegral(self, mchi):
        vectorresult = np.vectorize(self.SensiIntegral_scalar)
        return vectorresult(mchi)

    def logJbarPDF(self, logJbartestvalue):
        """ The (quasi-lognormal) distribution of log10(Jbar).
        Parameters:
        self.Jbar -> central/expectation value
        self.JbarError -> Error (in dex?!)
        logJbartestvalue -> log(variable of PDF)
        """
        return log((1./(sqrt(2*pi)*log(10)*(10.**logJbartestvalue)*self.JbarError) *
                exp(-0.5*(logJbartestvalue - log10(self.Jbar))**2/self.JbarError**2) ))
        #return 1.

    # MOST important:
    def logLhood(self, sigmav, mchi, jbar=None):
    #def logLhood(self, ParamArray)
        """Most important! This is each object's log likelihood function.

        Note: The data parameters Non, Noff, alpha are class members - hence they 
        need not be called as function parameters.
        Poissonian PMF: poisson.pmf(k,mu) = exp(-mu) * mu**k / k!

        ## ParamArray:
        ## [0] = sigmav
        ## [1] = mchi
        ## [2] = jbarArray
        """
        # Oldish:
        #Ns = ((sigmav / (8.*pi*mchi**2)) * self.Tobs * self.Jbar *
        #      self.SensiIntegralArray[mchi] )
        # Now: make logLhood a function of jbar!
        if jbar is None: # because self.Jbar as default argument doesn't work
            jbar=self.Jbar
        if np.isnan(mchi):
            print 'mchi is nan:', mchi
            sys.exit()
        ## if self.SensiIntegralArray[int(mchi)] < 1.:
        ##     print 'mchi = ', mchi
        ##     print 'sensi: ', self.SensiIntegralArray[int(mchi)]
            
        Ns = ((sigmav/(8.*pi*mchi**2)) * self.Tobs * jbar *
              self.SensiIntegralArray[int(mchi)] )

        # Use log PMF for faster (?) calculation:
        logPois1 = poisson.logpmf(self.Non, (Ns+self.alpha*self.Noff))
        #logPois2 = poisson.logpmf(self.Noff, self.Noff)
        # (moved to __init__ for faster running
        return logPois1+self.logPois2
        
    
    def ULsigmav(self,mchi):
        """ Calculates UL on sigmav with the object's member spectrum."""
        prefactor = 8.*pi*(mchi**2)*self.Nul/(self.Tobs*self.Jbar)
        result = prefactor / quad(lambda E: (self.Spectrum(E,mchi)*self.Aeff(E)),
                                  self.Eth, 1.01*mchi,limit=50,full_output=1)[0]
        return result 

    
class Pub(object):
    """This class holds the actual limits from publications."""

    def __init__(self, name, tablefile, legend):
        self.Name = name
        self.table = np.genfromtxt(tablefile)
        self.legend = legend
        print self.Name+" imported\n"
        self.mchis=self.table[:,0]
        self.ul=self.table[:,1]

    
