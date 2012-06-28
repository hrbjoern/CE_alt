#!/usr/bin/python

# Hello world 
print 'Hekki hekki hekki patang.\n'
#print "IMPORTANT: Please enter all numbers as floats!\n"


############################################################
# Imports:
############################################################

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cPickle as pickle
from math import pi
from scipy.integrate import quad

from Classes import Object, Pub
import PhotonSpectra
import Plotting
import Limits

print "Importing complete\n"

############################################################
# Objects:
#    def __init__(self, Name, Aeff, Tobs(h), Jbar, Non, Noff, alpha):
#    Numbers: See Log VII, p.26
############################################################

#Segue1M = Object("Segue1M", "Aeffs/MAGIC_Gaug_Aeff.dat",29.42,1.14e19,365.5)
# Including Eth = 100 GeV:
Segue1M = Object("Segue1M", "Aeffs/MAGIC_Gaug_Aeff.dat",
                 Tobs=29.42,
                 Jbar=1.14e19,
                 Non=52978,
                 Noff=53301,
                 alpha=1.0)
Segue1M.printObject()

Segue1V = Object("Segue1V","Aeffs/VERITAS_Segue1_Aeff_TrueE.dat", #"Aeffs/VERITAS-Aeff_20deg.dat",
                 Tobs=47.8,
                 Jbar=7.7e18,
                 Non=1082,
                 Noff=12479,
                 alpha=0.084)
Segue1V.printObject()

Segue1V_tautau = Object("Segue1V_tautau","Aeffs/VERITAS_Segue1_Aeff_TrueE.dat", #"Aeffs/VERITAS-Aeff_20deg.dat",
                        Tobs=47.8,
                        Jbar=7.7e18,
                        Non=1082,
                        Noff=12479,
                        alpha=0.084)
Segue1V.printObject()

SculptorIso = Object("SculptorIso","Aeffs/HESS_Aeffs_Crab_20deg.dat",
                     Tobs=11.8,
                     Jbar=2.9e17,
                     Non=117,
                     Noff=2283,
                     alpha=0.04)
SculptorIso.printObject()

SculptorNFW = Object("SculptorNFW","Aeffs/HESS_Aeffs_Crab_20deg.dat",
                     Tobs=11.8,
                     Jbar=2.9e17,
                     Non=117,
                     Noff=2283,
                     alpha=0.04)                    
SculptorNFW.printObject()

Sgr = Object("Sgr","Aeffs/HESS_Aeffs_Crab_20deg.dat",
             Tobs=11.,
             Jbar=1e18,
             Non=437,
             Noff=4270,
             alpha=(1./10.1))
Sgr.printObject()

Willman1V = Object("Willman1V", "Aeffs/VERITAS-Aeff_20deg.dat",
                   Tobs=13.68,
                   Jbar=84.3e17,
                   Non=326,
                   Noff=3602,
                   alpha=(1./11.))
Willman1V.printObject()

ObjList = (Segue1M, Segue1V, Segue1V_tautau, SculptorIso, Willman1V, Sgr)
#ObjList = (Segue1M,Segue1V)
#ObjList = (Segue1M,Segue1M)

# Consider all energies in GeV!!!
print "\nAll energies are in GeV\n"

# IMPORTANT. Basic energy value spacing.
try:
    esteps = int(sys.argv[1])
except:
    esteps = 100
print '\nenergy steps =', esteps
print
energies = np.logspace(2,5,esteps)

# EVEN MORE important: Combined upper limit(s).
UL_bbbar = []

SumOfTobsTimesJbar = 0.
SumOfNul = 0.
SumOfStuff = 0.
SumOfNon = 0.
SumOfNoff = 0.
SumOfAlphaNoff = 0.

CombinationList = (Segue1M, Segue1V, SculptorIso, Willman1V,Sgr)
#CombinationList = (Segue1M, Segue1V, SculptorIso, Willman1V)
CombinationList = (Segue1V, Segue1V, Segue1V, Segue1V)
#CombinationList = (Segue1M, Segue1M)

for o in CombinationList:
     SumOfNul += o.Nul
     SumOfNon += o.Non
     SumOfNoff += o.Noff
     SumOfAlphaNoff += o.alpha*o.Noff

# Calculate "mean" alpha:
alphabar = SumOfAlphaNoff / SumOfNoff

# Calculate "mean" sum of Nul:
SumOfSignal = SumOfNon - alphabar*SumOfNoff

print 'SumOfNon = ', SumOfNon
print 'SumOfNoff = ', SumOfNoff
print 'alphabar = ', alphabar
print 'SumOfNul = ', SumOfNul
print 'SumOfSignal = ', SumOfSignal

# Calculate "mean" total upper limit:
SumOfNulbar = Limits.UpperLimit(int(SumOfNon), int(SumOfNoff), alphabar)
print 'SumOfNulbar = ', SumOfNulbar

#sys.exit()
    
def IntegrandSum(E):
    soa = 0.
    for o in CombinationList:
        soa += o.Tobs*o.Jbar*o.Aeff(E)
    return soa

#SumOfStuff = SumOfNul/SumOfTobsTimesJbar


for mchi in energies:
    UL_bbbar.append(8.*pi*mchi**2*SumOfNulbar/
                    quad(lambda E: (PhotonSpectra.bbbar(E,mchi)*IntegrandSum(E)),
                         30., 1.01*mchi, limit=50,full_output=1)[0])

## print
## print 'UL_bbbar = ', UL_bbbar
## print

CombList = (UL_bbbar)

#sys.exit()

# Calculate individual upper limits for different spectra:
for o in ObjList:
    o.ul_tautau=[]
    o.ul_bbbar=[]
    o.ul_WW=[]
    for mchi in energies:
        o.ul_tautau.append(o.ULsigmav_tautau(mchi))
        o.ul_bbbar.append(o.ULsigmav_bbbar(mchi))
        o.ul_WW.append(o.ULsigmav_WW(mchi))




########################################################
# Data from publications:
########################################################
Segue1M_bb = Pub("Segue1M_bb","Publications/MAGIC_Segue1_ULsigmav_bb.dat",
                 "Segue1/MAGIC, bbar")

Segue1V_bb = Pub("Segue1V_bb","Publications/VERITAS_Segue1_ULsigmav_bb.dat",
                 "Segue1/VERITAS, bbar")

Segue1V_tautau = Pub("Segue1V_tautau","Publications/VERITAS_Segue1_ULsigmav_tautau.dat",
                 "Segue1/VERITAS, tautau")

Sculptor_IsoWW = Pub("Sculptor_IsoWW", "Publications/HESS_Scu-Limits_IsoBeta05.dat",
                   "Sculptor/HESS, Iso profile")

Willman1V_bb = Pub("Willman1V_bb", "Publications/VERITAS_Dwarfs_Willman1.dat",
                   "Willman 1/VERITAS")
                   

PubList = (Segue1M_bb, Segue1V_bb, Segue1V_tautau, Sculptor_IsoWW, Willman1V_bb)


########################################################
# Pickle stuff:
########################################################

print 'Start pickling:'
print

pickle.dump(energies, open('saveE.p','wb'))
## for i in ObjList:
##     i.printObject()
##     print
pickle.dump(ObjList, open('saveObj.p', 'wb'))
pickle.dump(PubList, open('savePub.p', 'wb'))
pickle.dump(CombList, open('saveComb.p', 'wb'))

########################################################
# Plot stuff:
########################################################

#Plotting.plotAeffs()
Plotting.plotObjects()
    

print "\nund tschuess \n"
sys.exit()

##########################################################


############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################


