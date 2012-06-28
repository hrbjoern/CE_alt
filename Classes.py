from PhotonSpectra import *
from math import pi, exp
from scipy.integrate import quad
from scipy import interpolate
import numpy as np
import fnmatch
from Limits import UpperLimit

class Object(object):
    """This class is supposed to hold the different \"objects\", i.e. 
    specific (publications of) DM searches towards a certain
    astrophysical object, and the formulae to calculate upper limits
    from their results."""
    
    def __init__(self, Name, Aeff, Tobs, Jbar, Non, Noff, alpha):
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
        self.Jbar = Jbar
        self.Non = Non
        self.Noff = Noff
        self.alpha = alpha
        self.Nul = float(UpperLimit(self.Non, self.Noff, self.alpha))
                
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

    def ULsigmav_tautau(self,mchi):
        """ Calculates UL on sigmav with a Cembranos tautau spectrum."""
        prefactor = 8.*pi*(mchi**2)*self.Nul/(self.Tobs*self.Jbar)
#        result = prefactor / quad(lambda E: (Bergstrom1998(E,mchi)*self.Aeff(E)),
#                                  self.Eth, 1.01*mchi,limit=50,full_output=1)[0]
        result = prefactor / quad(lambda E: (tautau(E,mchi)*self.Aeff(E)),
                                  self.Eth, 1.01*mchi,limit=50,full_output=1)[0]
#        return np.minimum(result,1.)
        return result
    
    def ULsigmav_bbbar(self,mchi):
        """ Calculates UL on sigmav with a Cembranos bbbar spectrum."""
        prefactor = 8.*pi*(mchi**2)*self.Nul/(self.Tobs*self.Jbar)
        result = prefactor / quad(lambda E: (bbbar(E,mchi)*self.Aeff(E)),
                                  self.Eth, 1.01*mchi,limit=50,full_output=1)[0]
        return result     
    
    def ULsigmav_WW(self,mchi):
        """ Calculates UL on sigmav with a standard Bergstrom spectrum."""
        prefactor = 8.*pi*(mchi**2)*self.Nul/(self.Tobs*self.Jbar)
        result = prefactor / quad(lambda E: (Bergstrom1998(E,mchi)*self.Aeff(E)),
                                  self.Eth, 1.01*mchi,limit=50,full_output=1)[0]
        return result

    
    
class Pub(object):
    """This class holds the actual limits from publications."""

    def __init__(self, name, tablefile, legend):
        self.Name = name
        self.table = np.genfromtxt(tablefile)
        self.legend = legend
        print self.Name+" imported\n"
        self.energies=self.table[:,0]
        self.ul=self.table[:,1]

    
