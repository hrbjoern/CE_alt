from PhotonSpectra import *
from math import pi, exp
from scipy.integrate import quad
from scipy import interpolate
import numpy as np


class Object(object):
    """This class is supposed to hold the different \"objects\", i.e. 
    specific (publications of) DM searches towards a certain
    astrophysical object, and the formulae to calculate upper limits
    from their results."""
    
    def __init__(self, name, aeff, tobs, jbar,nul):
        """ Initializes \"object\", including some specific object modifications."""
        print "\nHallo", name
        self.Name = name
        self.Eth = 30.
        self.Aeffdata = np.genfromtxt(aeff)
        # Consider all energies in GeV, all Aeffs in cm**2:
        self.Aeffdata[:,1]*=1e4
        if self.Name == "Segue1M":
            self.Aeffdata[:,1]*=100.
            self.Eth = 100.
        elif self.Name == "Segue1V" or "Sculptor" or "Sgr": # Works! :)
            self.Aeffdata[:,0]*=1000.
        self.Tobs = tobs*3600.
        self.Jbar = jbar
        self.Nul = nul
                
    def printObject(self):
        """ Prints object, just for testing. """
        print self.Name
        print self.Tobs, ", Tobs in sec"
        print self.Jbar, "GeV^2 cm^-5"
        print self.Nul, ", Nul at some CL"
        print self.Aeffdata[0,0], ", Eth in GeV"


    def Aeff(self,Energy):
        """ Aeff function from the Aeffdata stored in each object."""
        if (self.Aeffdata[0,0] < Energy < self.Aeffdata[-1,0]):
            return float(interpolate.interp1d(self.Aeffdata[:,0],self.Aeffdata[:,1])(Energy))
        else:
            return 0.

    def ULsigmav(self,mchi):
        """ Calculates UL on sigmav with a standard Bergstrom spectrum."""
        prefactor = 8*pi*mchi**2*self.Nul/(self.Tobs*self.Jbar)
        result = prefactor / quad(lambda E: (Bergstrom1998(E,mchi)*self.Aeff(E)), self.Eth, 1.1*mchi,
                                  limit=50,full_output=1)[0]
        #        return np.minimum(result,1.)
        return result
        
    
        
class Pub(object):
    """This class holds the actual limits from publications."""

