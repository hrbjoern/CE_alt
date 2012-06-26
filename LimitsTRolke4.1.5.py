#!/usr/bin/env python
#
# Python example for Limit Calculations and hypothesis tests
# with the Profile Likelihood method, using TRolke2
# Johan Lundberg, johan.lundberg@cern.ch -- June 2009, CERN
#
#
# Adapted by HB: 2012-05-02
#  - using efficiency = 1.0 at the moment
#

import sys

print "TRolke limits from python"
print "Poissonian background, known efficiency \'eff\' "
print "-------------------------------------"

try:
    import ROOT
    print "pyRoot loaded" 
except:
    print "Could not import pyRoot !"
    print "   Make sure you added root to your python path "
    print '   export PYTHONPATH="$ROOTSYS"/lib:"$PYTHONPATH"      '
    print "   Also, make sure your root installation contains pyRoot"
    print "   e.g. you should have the file $ROOTSYS/lib/libPyROOT.so"  
    sys.exit()

try:
    ROOT.gSystem.Load("libTRolke.1.so")
    print "TRolke library loaded" 
except:
    print "TRolke library not found"
    print "  ... this may not work, going on."

######################################################
# declare variables for return types
# -------
# this trick is needed since in python integers and 
# floats are not mutable (i.e. they can not change in place)
#
lower=ROOT.Double() 
upper=ROOT.Double()  
ncrit=ROOT.Long() 
LDS=ROOT.Double()  
#
######################################################


######################################################
#
# make trolke object named tr

tr = ROOT.TROLKE2.TRolke()
#
# if you would like to use the TRolke version delivered
# with Root: change this to 
#   tr = ROOT.TRolke()
#
######################################################

print " "
print "////////////////////////////////////////////////////////"
print "// Limit calculation "
print "////////////////////////////////////////////////////////"

#-------------------------------
## Prepare the model parameters.
#-------------------------------

# Defaults:
x = 160       # events in the signal region
y = 122       # events in background region
tau = 1.      # ratio between size of signal/background region
alpha = 0.95; # Confidence Level
e = 1.0       # efficiency
# (Notation is (c) J. Lundberg..!)


# Command line arguments:
if len(sys.argv)!=6:
    print "Please provide the parameters: "
    print "python LimitsRolkePoisson5.py Non Noff ConfLevel alpha eff"
    print
    print "Using default values:"
else:
    x = int(sys.argv[1])
    y = int(sys.argv[2])
    alpha = float(sys.argv[3])
    tau = float(sys.argv[4])
    e = float(sys.argcv[5])
print
print 'Non = ', x
print 'Noff = ', y
print 'C.L. = ', alpha
print 'alpha = ', tau
print 'eff = ', e
print 


# Select model no. 5:
# ------------
tr.SetPoissonBkgKnownEff(x,y,tau,e)
tr.SetCL(alpha);

# calculate limits
# ----------------
status = 0
status=tr.GetLimits(lower,upper);

# print limits
# ------------
if(status):
    print "The Profile Likelihood interval is :" 
    print "[%.4f,%.4f]" % (lower,upper)

else:
    print "Calculation failed. (This should not be)"



print
print
