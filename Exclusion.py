#!/usr/bin/python

# Hello world 
print 'Hekki hekki hekki patang.'
print
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
from scipy.optimize import fmin, minimize, fmin_tnc, brute # minimization (scipy >=0.11!)
from scipy.optimize import fsolve, root, brentq # root finding
from scipy.interpolate import interp2d, interp1d
from scipy import interp
from pprint import pprint

#import pyximport; pyximport.install(pyimport = True)
#import pyximport; pyximport.install()
from Classes import Object, Pub
import PhotonSpectra
import Plotting
import Limits

print "Importing complete:", time.ctime()
print

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


# Calculate individual upper limits:
for o in ObjList:
    o.UL = []
    for mchi in mchis:
        o.UL.append(o.ULsigmav(mchi))
        ## o.ul_tautau.append(o.ULsigmav_tautau(mchi)) ... deprecated ;)



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

CombinationList = (Segue1M, Segue1V, SculptorIso, Willman1V, Sgr)
## print '\nCombinationList =', CombinationList
## for o in CombinationList:
##     print CombinationList.index(o)

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

# Vectorize the integration:
def integral(mchi):
    return quad(lambda E: (PhotonSpectra.bbbar(E,mchi)*IntegrandSum(E)),
                30., 1.01*mchi, limit=50,full_output=1)[0]

vec_integral = np.vectorize(integral)

# Array-ify the calculation: works! :) Replaces "for mchi in mchis". 
UL_bbbarComb = (8.*pi*mchis**2*SumOfNulbar/ vec_integral(mchis) )

print '\nUL_bbbarComb = ', UL_bbbarComb


########################################################
# "Likelihood" Combination of the limits:
########################################################
# (still needs "CombinationList" from above!)
print '\nStart Combination: \n'

# Timing tests:
class Timer():
   def __enter__(self): self.start = time.time()
   def __exit__(self, *args): print time.time() - self.start


print 'Sensi integration time:'
with Timer():
    # Define array of integrated sensitivities:
    for o in CombinationList:
        for mchi in mchis:
            o.SensiIntegralArray[int(mchi)] = o.SensiIntegral(mchi)
            #print 'o.SensiIntegralArray[%.2f] = %.2f' %(mchi, o.SensiIntegralArray[mchi])

    
# Jbar testing list:
JbarList = range(len(CombinationList))
for i in JbarList:
    #JbarList[i] = 1e21
    JbarList[i] = 19. # (using log values instead!)
    

# IMPORTANT:
def ComblogLhood(ParArray, mchi):
    """ Combined (minus) log likelihood. THE master function.
    ## ParArray:
    ## [0] = log10(sigmav)
    ## [1...] = log10jbarArray -> Note: Index+1 in sum!
    """
    # Using logarithmic parameters:
    return -sum(o.logLhood(10.**ParArray[0], mchi, 10.**ParArray[CombinationList.index(o)+1])
                +o.logJbarPDF(ParArray[CombinationList.index(o)+1]) for o in CombinationList)

def ComblogLhood_sv(JbarArray, log10sigmav, mchi):
    """ Combined (minus) log likelihood. The master function, as a function of sigmav.
    ## JbarArray:
    ## [0...] = log10jbarArray -> Note: NO index+1 in sum!
    """
    # Using logarithmic parameters:
    return -sum(o.logLhood(10.**log10sigmav, mchi, 10.**JbarArray[CombinationList.index(o)])
                +o.logJbarPDF(JbarArray[CombinationList.index(o)]) for o in CombinationList)

TestArray = np.append(-21., JbarList)
print 'TestArray = ', TestArray
    
print 'CL calling time: '
with Timer():
    clout1 = ComblogLhood(TestArray, mchis[-2])
print 'clout1 = ', clout1 
print
#sys.exit()


########################################################
# Start minimization:
########################################################

sigmavTestRange = np.arange(-26, -18, (26.-18.)/esteps)
JBA = np.array(JbarList) # 19's
JBbds = []
for num in JbarList:
    JBbds.append((16., 22.))



print '\n Testing (and plotting) the profile likelihood:', time.ctime()
#print '(not right now)\n'
# Tolerance for minimization attempts:
mintol = 0.2

# 2-dim. CL array:
CLmins = np.empty((len(mchis), len(sigmavTestRange)))

for mchi in mchis:
    #print 'mchi = ', mchi
    mindex = (mchis.tolist()).index(mchi)
    for sv in sigmavTestRange:
        #print 'sigmav = ', sv
        svindex = (sigmavTestRange.tolist()).index(sv)
        minresult4 = minimize(ComblogLhood_sv, JBA, args=(sv, mchi,), bounds=JBbds, 
                              method='TNC').fun
        minresult5 = minimize(ComblogLhood_sv, 0.9*JBA, args=(sv, mchi,), bounds=JBbds, 
                              method='TNC').fun
        if abs(minresult4-minresult5) > mintol:
            minresult6 = minimize(ComblogLhood_sv, 1.1*JBA, args=(sv, mchi,), bounds=JBbds, 
                                  method='TNC').fun
            minimum = min((minresult4, minresult5, minresult6))
        else:
            minimum = min((minresult4, minresult5))
                          
        CLmins[mindex][svindex] = minimum


print '\nCLmins array:'
print '(not printed anymore)\n'
#pprint(CLmins)

pickle.dump(CLmins, open('saveCLmins.p','wb'))


#sys.exit()

print
print


# Find sigmav values where  Delta(-2 ln PL) = 2.706 :
print 'Find sigmav values where  Delta(-2 ln PL) = 2.706 : PRUEFEN!!!\n'
CLmin = minimize(ComblogLhood_sv, JBA, args=(-26, 2000.,), bounds=JBbds, 
                  method='TNC').fun
print 'For logsv =-26, mchi=2000:'
print 'CLmin = ', CLmin


def findRoot(sigmav, mchi):
    CLaktuell = minimize(ComblogLhood_sv, JBA, args=(sigmav, mchi,), bounds=JBbds, 
                         method='TNC').fun
    # Test the minimization result:
    CLaktuell2 = minimize(ComblogLhood_sv, JBA*0.9, args=(sigmav, mchi,), bounds=JBbds, 
                         method='TNC').fun
    if abs(CLaktuell-CLaktuell2) > mintol:
        CLaktuell3 = minimize(ComblogLhood_sv, JBA*1.1, args=(sigmav, mchi,), bounds=JBbds, 
                              method='TNC').fun
        # Take the real minimum:
        CLaktuell = min((CLaktuell, CLaktuell2, CLaktuell3))
    else:
        # Take the real minimum:
        CLaktuell = min((CLaktuell, CLaktuell2))
    CLdiff = CLaktuell-CLmin-(2.706/2.)
    ## if abs(CLdiff) < 0.1:
    ##     print 'CLaktuell = ', CLaktuell
    return CLdiff

# Array of sigmav limit values:
UL_CombLhood = np.zeros_like(mchis)

print '\nStarting to find the upper limits:', time.ctime()

for mchi in mchis:
    #if np.log10(mchi)%1.==0:
        #print 'mchi = ', mchi
    sigmav_UL = brentq(findRoot, -26, -17, # bracketing values
                       args=(mchi,),
                       #full_output=True,
                       #maxfev = 1000, 
                       # factor=0.1
                       )
    #print 'sigmav_UL = ', sigmav_UL

    # Fill array of sigmav limits:
    UL_CombLhood[(mchis.tolist()).index(mchi)] = 10.**(sigmav_UL)

print '\nResult: Limits from full likelihood', time.ctime()
print 'UL_CombLhood = ', UL_CombLhood


# List of combined limits: entries are arrays 
CombList = (UL_bbbarComb, UL_CombLhood)





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



## ## 2D:
## # Interpol function:
## CLinterpol = interp2d(sigmavs, mchis, CLA_cut, kind='linear')
## # Test array:
## #print 'CLinterpol(sigmavs, mchis):'
## CL_Test = CLinterpol(sigmavs, mchis)
## #print CL_Test
## #print

## TestDiff = np.where((np.abs(CL_Test-CLA_cut)/CLA_cut)>0.1,
##                     np.abs(CL_Test-CLA_cut), 0.)
## print '\nLarge differences (>10%) in interpolation:'
## #print np.nonzero(TestDiff)
## print
## Fehler = len(np.nonzero(TestDiff)[0])
## von = len(CLA_cut.flat)
## Prozent = 100.*float(Fehler)/float(von)
## print 'This occured %i (out of %i) times, i.e. %.1f%% of the time.\n' \
##       % (Fehler, von, Prozent)


############################################################
# Alter Minimierungsgedoens:
############################################################

## print '\nCL Array calling time: '
## with Timer():
##     CLArray = CL_vec(mv, sv, JbarList)
## ## print '\n CLArray(mv, sv):'
## ## print CLArray
## ## print

## # Interpolation:

## # Cut off infinities in CLArray:
## CLA_cut = np.where(CLArray < 100., CLArray, 100.)
## ## print '\n CLArray, neu: CLA_cut'
## ## print CLA_cut
## ## print 'CLA_cut.shape: ', CLA_cut.shape
## ## print

## # Array of sigmav limit values:
## UL_CombLhood = np.zeros_like(mchis)


## # Bounds:
## bds = (np.atleast_1d(np.array([0.])), np.atleast_1d(np.array([1e-15])))
## print '\nMinimization bds = ', bds

## ## 1D along each mchi:
## for m in np.arange(mchis.size): # loop over "mchis" indices, not entries
##     CLmin_m = CLA_cut[:, m].min()
##     #CLmax_m = CLA_cut[:, m].max()

##     def CLinterpol(sigmav):
##         return interp(sigmav, sigmavs, CLA_cut[:, m]) # this is the right index!
##     ## for s in sigmavs:
##     ##     print 'sigmav, CLinterpol(sigmav):', s, CLinterpol(s)
##     #sigmav_min = fmin(CLinterpol, 1e-23) # returns x, not f(x)!
##     #sigmav_min = minimize(CLinterpol, 1e-23, method="L-BFGS-B",

##     sigmav_min = minimize(CLinterpol, np.atleast_1d(np.array(1e-23)), method="TNC")#,
##                           #bounds=bds) # returns x, not f(x)!

##     def CLsolve(s):
##         DeltalogL = 2.71 #..???????
##         return CLinterpol(s)-DeltalogL-CLmin_m # fsolve for f(x) = 0
##     sigmav_minplus2 = fsolve(CLsolve, 1e-23)

##     # Fill array of sigmav limits:
##     UL_CombLhood[m] = sigmav_minplus2



## # List of combined limits: entries are arrays 
## CombList = (UL_bbbarComb, UL_CombLhood)


############################################################
# Nicht ganz so altes Minimierungsgedoens:
############################################################

# Minimization test:
#print '\nMinimization test:'
#print '(not any longer)\n'
## print 'Starting values:'
## print TestArray

## # Bounds: List of tuples!
## bds = [(-26., -20.)] # sigmav bounds
## for num in JbarList:
##     bds.append((16., 22.))
## print 'Bounds = ', bds

## minresult = minimize(ComblogLhood, TestArray, args=(mchis[-2],), bounds=bds, 
##                      #options={'maxiter':int(1e2), 'disp': True},
##                      method='TNC')
## print '\nminresult.x = ', minresult.x
## print 'minresult.message = ', minresult.message
## print 'minresult.fun = ', minresult.fun

## print '\nMinimization test 2:'
## print 'Starting values:'
## TestArray2 = np.append(-22., JbarList)
## print TestArray2

## minresult2 = minimize(ComblogLhood, TestArray2, args=(mchis[-2],), bounds=bds, 
##                      #options={'maxiter':int(1e2), 'disp': True},
##                      method='TNC')
## print '\nminresult2.x = ', minresult2.x
## print 'minresult2.fun = ', minresult2.fun
## ## minresult2 = fmin_tnc(ComblogLhood, TestArray2, args=(mchis[-3],), approx_grad=1,
## ##                       epsilon=1., bounds=bds, disp=5, accuracy=0.01, ftol = 0.1)
## ## print 'minresult2 = ', minresult2

## print '\nMinimization test 3:'
## print
## print 'For sigmav in (-26, -19), mchi=mchis[-2] ...:\n'

## for sv in sigmavTestRange:
##     print 'sigmav = ', sv
##     minresult3 = minimize(ComblogLhood_sv, JBA, args=(sv, mchis[-2],), bounds=JBbds, 
##                   method='TNC')
##     print 'minresult3.x = ', minresult3.x
##     print 'minresult3.fun = ', minresult3.fun
##     print


## print '\nMinimization test 4:' ## Resultat: in 2012-08-28_Exclusion_Minimize-out.txt
## print
## print 'For sigmav in (-26, -19), mchi in ...:\n'


############################################################
    
## # Vectorize the comb. lhood function? Yes!
## CL_vec = np.vectorize(ComblogLhood)

## # np.array for sigmav values: (over which to interpolate)
## sigmav_steps = esteps # make same-dim array
## sigmavs = np.logspace(-26, -20, sigmav_steps)
## #print 'sigmavs = ', sigmavs

## # 2D-valued array of sigmavs and mchis: Works! :)
## mv, sv = np.meshgrid(mchis, sigmavs) # Note: indexing rules!


## ## print 'meshgrid(s): mv, sv'
## ## print mv
## ## print sv
## ## print



#sys.exit()
