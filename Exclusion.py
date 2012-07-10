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
#import timeit
import time

from math import pi
from scipy.integrate import quad
from scipy.optimize import fmin
from pprint import pprint

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
                 alpha=1.0,
                 spectrum="bbbar")
#Segue1M.printObject()

Segue1V = Object("Segue1V","Aeffs/VERITAS_Segue1_Aeff_TrueE.dat", #"Aeffs/VERITAS-Aeff_20deg.dat",
                 Tobs=47.8,
                 Jbar=7.7e18,
                 Non=1082,
                 Noff=12479,
                 alpha=0.084,
                 spectrum="bbbar")
#Segue1V.printObject()

Segue1V_tautau = Object("Segue1V-tautau","Aeffs/VERITAS_Segue1_Aeff_TrueE.dat", #"Aeffs/VERITAS-Aeff_20deg.dat",
                        Tobs=47.8,
                        Jbar=7.7e18,
                        Non=1082,
                        Noff=12479,
                        alpha=0.084,
                        spectrum="tautau")
#Segue1V.printObject()

SculptorIso = Object("SculptorIso","Aeffs/HESS_Aeffs_Crab_20deg.dat",
                     Tobs=11.8,
                     Jbar=2.9e17,
                     Non=117,
                     Noff=2283,
                     alpha=0.04,
                     spectrum="bbbar")
#SculptorIso.printObject()

SculptorNFW = Object("SculptorNFW","Aeffs/HESS_Aeffs_Crab_20deg.dat",
                     Tobs=11.8,
                     Jbar=2.9e17,
                     Non=117,
                     Noff=2283,
                     alpha=0.04,
                     spectrum="bbbar")                    
#SculptorNFW.printObject()

Sgr = Object("Sgr","Aeffs/HESS_Aeffs_Crab_20deg.dat",
             Tobs=11.,
             Jbar=1e18,
             Non=437,
             Noff=4270,
             alpha=(1./10.1),
             spectrum="bbbar")
#Sgr.printObject()

Willman1V = Object("Willman1V", "Aeffs/VERITAS-Aeff_20deg.dat",
                   Tobs=13.68,
                   Jbar=84.3e17,
                   Non=326,
                   Noff=3602,
                   alpha=(1./11.),
                   spectrum="bbbar")
#Willman1V.printObject()

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
print '\n # energy steps =', esteps
print
mchis = np.logspace(2,5,esteps)

## # TESTING:
## print 'SensiIntegral(mchis): '
## for o in ObjList:
## ##     #o.Sensi = np.vectorize(o.Sensi_scalar)
## ##     #pprint (o.Sensi(mchis))
##     pprint (o.SensiIntegral(mchis))
## sys.exit()


# Calculate individual upper limits:
for o in ObjList:
    o.ul = []
    for mchi in mchis:
        o.ul.append(o.ULsigmav(mchi))
        ## o.ul_tautau.append(o.ULsigmav_tautau(mchi)) ... deprecated ;)

    #o.ul2 = (np.vectorize(o.ULsigmav))(mchis) #... funktioniert leider nicht
    #o.ul2 = o.ULsigmav(mchis)
    #pprint(o.ul2)
    # pprint (vars(o)) # Tolltolltoll!!! :)

# Add sensitivity curve S(E) to all objects
# - but only for selected masses!

#for o in ObjList:
#    def o.S(E)
    
    


########################################################
# "Simple" Combination of the limits:
########################################################

UL_bbbarComb = []

SumOfTobsTimesJbar = 0.
SumOfNul = 0.
SumOfStuff = 0.
SumOfNon = 0.
SumOfNoff = 0.
SumOfAlphaNoff = 0.

CombinationList = (Segue1M, Segue1V, SculptorIso, Willman1V,Sgr)
#CombinationList = np.array([Segue1M, Segue1V, SculptorIso, Willman1V,Sgr], dtype=np.object)
# ... doesnt make a difference.
#CombinationList = (Segue1M, Segue1V, SculptorIso, Willman1V)
#CombinationList = (Segue1V, Segue1V, Segue1V, Segue1V)
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

# Calculate "mean" total upper limit of photons:
SumOfNulbar = Limits.UpperLimit(int(SumOfNon), int(SumOfNoff), alphabar)
print 'SumOfNulbar = ', SumOfNulbar


# Calculate "mean" total upper limit on sigmav:
def IntegrandSum(E):
    soa = 0.
    for o in CombinationList:
        soa += o.Tobs*o.Jbar*o.Aeff(E)
    return soa

## UL_bbbarComb = [] # ... unnecessary! See array-ified version below.
## for mchi in mchis:
##     UL_bbbarComb.append(8.*pi*mchi**2*SumOfNulbar/
##                     quad(lambda E: (PhotonSpectra.bbbar(E,mchi)*IntegrandSum(E)),
##                          30., 1.01*mchi, limit=50,full_output=1)[0])
##     # Testing:
##     #for o in ObjList:
##     #    print 'logLikelihood(sigmav=0, mchi) = ', o.logLhood(0., mchi)

# Vectorize the integration:
def integral(mchi):
    return quad(lambda E: (PhotonSpectra.bbbar(E,mchi)*IntegrandSum(E)),
                30., 1.01*mchi, limit=50,full_output=1)[0]

vec_integral = np.vectorize(integral)

# Array-ify the calculation: works! :) Replaces "for mchi in mchis". 
UL_bbbarComb = (8.*pi*mchis**2*SumOfNulbar/ vec_integral(mchis) )



#print 'UL_bbbarComb = ', UL_bbbarComb


########################################################
# "Likelihood" Combination of the limits:
########################################################
# (still needs "CombinationList" from above)
print '\nStart Combination: \n'

def ComblogLhood(sigmav, mchi):
    """ Combined log likelihood. THE master function."""
    return -sum(o.logLhood(sigmav, mchi) for o in CombinationList)

# preliminary maximization:
#CLmax = fmin(ComblogLhood, 10., args=[1000.])
#print 'CLmax = ', CLmax

# Timing test:
class Timer():
   def __enter__(self): self.start = time.time()
   def __exit__(self, *args): print time.time() - self.start

print 'CL calling time: '
with Timer():
    ComblogLhood(1e-22, 1000.)
print 
    
#t = timeit.Timer('ComblogLhood(1e-22, 1000.)', "from __main__ import ComblogLhood")
#print '\n time it:', t.timeit()
#sys.exit()

# LEFTOVER: move elsewhere!
# List of combined limits: Add more later!
CombList = (UL_bbbarComb)


# Vectorize the comb. lhood function? Yes!
CL_vec = np.vectorize(ComblogLhood)
#CL_vec = ComblogLhood

# np.array for sigmav values: (over which to interpolate?)
sigmav_steps = esteps # for the moment
sigmavs = np.logspace(-26, -20, sigmav_steps)
print 'sigmavs = ', sigmavs

# 2D-valued array of sigmavs and mchis: Works! :)
sv, mv = np.meshgrid(sigmavs, mchis)

# Resulting array of CL values
#print '\n\n CL-Array = CL_vec(sigmavs, mchis):'
#CLArray = CL_vec(sv, mv)
#print CLArray


print '\nCL Array calling time: '
with Timer():
    CLArray = CL_vec(sv, mv)
print CLArray


## print '\nCL Array calling time, v2: '
## with Timer():
##     ComblogLhood(sv, mv)

sys.exit()


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

pickle.dump(mchis, open('saveE.p','wb'))
## for i in ObjList:
##     i.printObject()
##     print
pickle.dump(ObjList, open('saveObj.p', 'wb'))
pickle.dump(PubList, open('savePub.p', 'wb'))
pickle.dump(CombList, open('saveComb.p', 'wb'))

########################################################
# Plot stuff:
########################################################

print 'Start plotting: '
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


