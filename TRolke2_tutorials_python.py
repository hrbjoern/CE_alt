#!/usr/bin/env python
#
# Python example for Limit Calculations and hypothesis tests
# with the Profile Likelihood method, using TRolke2
#
# Johan Lundberg, johan.lundberg@cern.ch -- June 2009, CERN

import sys

print "Tutorial for TRolke usage from python"
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
print "// Example of limit calculation, Gauss, Gauss           "
print "////////////////////////////////////////////////////////"

#
# prepare the model parameters. We are going to use Gauss,Gauss

bm = 5;      # expected background
x = 10;      # events in the signal region
sdb = 0.5;   # standard deviation in background estimate
em = 0.9;    # measured efficiency
sde = 0.05;  # standard deviation of efficiency
alpha =0.99; # Confidence Level

# select model
# ------------
tr.SetGaussBkgGaussEff(x,bm,em,sde,sdb); 
tr.SetCL(alpha);

# calculate limits
# ----------------
status = 0
status=tr.GetLimits(lower,upper);

# print limits
# ------------
if(status):
    print "For model 3 : Gaussian / Gaussian" 
    print "the Profile Likelihood interval is :" 
    print "[%.4f,%.4f]" % (lower,upper)

else:
    print "Calculation failed. (This should not be)"



print " "
print "////////////////////////////////////////////////////////"
print "// Example of hypothesis tests                          "
print "////////////////////////////////////////////////////////"

bm=3.7;
sdb=0.51;
em=0.5;
sde=0.08;
alpha=0.99;
beta=0.5;
anything=0

tr.SetGaussBkgGaussEff(anything, bm, em, sde, sdb);
tr.SetCL(alpha);
tr.GetCriticalNumber(ncrit);
tr.GetLeastDetectableSignal(LDS, beta);


if(status):
    print " ---------------------------- " 
    print "  critical number:         " , ncrit
    print "  least detectable signal: " , LDS
    print " "
else:
    print "Calculation failed. (This should not be)"

print "  [Done!]"
print " "



