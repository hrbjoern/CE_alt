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

from Classes import Object, Pub

import Plotting

print "Importing complete\n"

############################################################
# Objects:
#    def __init__(self, name, aeff, tobs(h), jbar,nul):
#    Numbers: See Log VII, p.26
############################################################

#Segue1M = Object("Segue1M", "Aeffs/MAGIC_Gaug_Aeff.dat",29.42,1.14e19,365.5)
# Including Eth = 100 GeV:
Segue1M = Object("Segue1M", "Aeffs/MAGIC_Gaug_Aeff.dat",29.42,1.14e19,453)
Segue1M.printObject()

Segue1V = Object("Segue1V","Aeffs/VERITAS-Aeff_20deg.dat",47.8,7.7e18,135.9)
Segue1V.printObject()

SculptorIso = Object("SculptorIso","Aeffs/HESS_Aeffs_Crab_20deg.dat",11.8,2.9e17,32.4)
SculptorIso.printObject()

SculptorNFW = Object("SculptorNFW","Aeffs/HESS_Aeffs_Crab_20deg.dat",11.8,32.9e17,32.4)
SculptorNFW.printObject()

Sgr = Object("Sgr","Aeffs/HESS_Aeffs_Crab_20deg.dat",11.,1e18,56.)
Sgr.printObject()

Willman1V = Object("Willman1V", "Aeffs/VERITAS-Aeff_20deg.dat",13.68,84.3e17,18.8)
Willman1V.printObject()

ObjList = (Segue1M, Segue1V, SculptorIso, Willman1V, Sgr)
#ObjList = (Segue1M,Segue1V)

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
                 "Segue1/MAGIC, bbar, published")

Segue1V_bb = Pub("Segue1V_bb","Publications/VERITAS_Segue1_ULsigmav_bb.dat",
                 "Segue1/VERITAS, bbar, published")

Sculptor_IsoWW = Pub("Sculptor_IsoWW", "Publications/HESS_Scu-Limits_IsoBeta05.dat",
                   "Sculptor/HESS, Iso profile")

Willman1V_bb = Pub("Willman1V_bb", "Publications/VERITAS_Dwarfs_Willman1.dat",
                   "Willman 1/VERITAS")
                   

PubList = (Segue1M_bb,Segue1V_bb, Sculptor_IsoWW, Willman1V_bb)


########################################################
# Pickle stuff:
########################################################

pickle.dump(energies, open('saveE.p','wb'))
pickle.dump(ObjList, open('saveObj.p', 'wb'))
pickle.dump(PubList, open('savePub.p', 'wb'))


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


