from math import pi, exp, log


############################################################
# Cembranos Spectra:
##############################################

# General:
def dNdx(x,a1,b1,b2,c1,c2,d1,d2,q,p,n1,n2):
    return (x**-1.5) * (a1*exp(-b1*x**n1 -b2*x**n2
                               -c1/x**d1 +c2/x**d2)
                        +(q*(x**1.5)*log(p*(1.-x))
                          *((x**2-2*x+2)/x)))
                                                
# tautau:
def tautau(E, mchi):
    if (E < 10. or E>=mchi):
        return 0.
    else:
        x = E/mchi
        a1=14.7
        b1=5.4
        b2=5.31
        n2=1.4
        c1=2.54
        d1=0.295
        c2=0.373
        d2=0.47
        q=0.0026
        if (25 < mchi < 1e4):
            n1 = 10.6*mchi**-0.0148
            p = 0.773*mchi**1.75
        elif mchi >= 1e4:
            n1 = -7*mchi**-1.99 + 179.*mchi**-0.763 + 9.09
            p = 3.07*mchi**1.55
        return mchi**-1*dNdx(x,a1,b1,b2,c1,c2,d1,d2,q,p,n1,n2)

# bbbar:
def bbbar(E, mchi):
    if (E < 10. or E>=mchi):
        return 0.
    else:
        x = E/mchi
        a1=10.
        b1 = 152.*mchi**-0.462
        b2=11.0
        n1 = 18.7*mchi**-0.248
        c1 = 0.328*mchi**0.0447
        d1=0.295
        c2=0.0151
        d2=0.55
        q=2.6e-4
        p = 11.8*mchi**0.641
        if (25 < mchi < 1000):
            n2 = 0.805*mchi**-0.0319
        elif mchi >= 1e3:
            n2 = 0.707*mchi**-0.0129
        if (25 < mchi < 600):
            d1 = 0.474*mchi**-0.0639 + 37.1*mchi**-1.87
        elif mchi>=600:
            d1 = 0.449*mchi**-0.0552
        return mchi**-1*dNdx(x,a1,b1,b2,c1,c2,d1,d2,q,p,n1,n2)




############################################################
# "Old" Spectra:
##############################################

def Bergstrom1998(E,mchi):
    if E>mchi:
        return 0.
    else:
        dNdE = (1./mchi) * 0.73 * (E/mchi)**-1.5 * exp(-7.8*E/mchi) 
        #    dNdE = 0.73 * (E/mchi)**-1.5 * exp(-7.8*E/mchi) 
    return dNdE


def Bringmann2008(E,mchi):
    if E>mchi:
        return 0.
    else:
        x = E/mchi
        eps = 80./mchi
        dndE = ( (1./mchi) * (alpha/pi) * 4.*(1-x+x**2)**2 /
                 ((1-x+eps/2.)*x)*(log(2.*(1-x+eps/2.)/eps) -0.5 +x -x**3 ) )
        if dndE<0:
            dndE=0.

        dndE = dndE+( (1./mchi)*0.73*(E/mchi)**-1.5 *exp(-7.8*E/mchi) )
    return dndE

def DGLAP2009(E,mchi):
    if E>mchi:
        return 0.
    else:
        x = E/mchi
## Cirelli / Bergstrom: me vs. mmu?
##         dndE = (1./mchi) * ( (alpha/(2.*pi)) *
##                              (-1 + np.log((4.*mchi**2/me**2)*(1.-x)))
##                              *(1.+(1.-x)**2)/x )
##         dndE = (1./mchi) * ( (alpha/(2.*pi)) *
##                              (-1 + np.log((4.*mchi**2/mmu**2)*(1.-x)))
##                              *(1.+(1.-x)**2)/x )

## Bovy:
        dndE = ( (alpha/pi) * ((1.+(1.-x)**2)/E) *
                 (-1. + np.log(4.*(1.-x)) - 2.*np.log(mmu/mchi) ) )
                 

        if dndE<0:
            dndE=0.
    return dndE



def Delta(E, mchi):
    if ( fabs(E-mchi)<1.1):
        return 1.
    else:
        return 0.001


    

##############################################
# Integrated annihilation spectra:
##############################################


def BergInt(mchi):
    return quad(lambda E: Bergstrom1998(E,mchi), Eth, 1.1*mchi)[0]
#    result = quad(lambda E: BergstromSmoothed(E,mchi), Eth, 2.*mchi)[0]
#    print ('\ntype(lambda E: BergstromSmoothed(E,mchi)) = ',
#           type(lambda E: BergstromSmoothed(E,mchi)), '\n')
#    result = quadrature(lambda E: BergstromSmoothed(E,mchi), Eth, 2.*mchi,
#                        vec_func=False)[0]

def BringInt(mchi):
    timeBring_i = time()
    result = quad(lambda E: Bringmann2008(E,mchi), Eth, 1.1*mchi,
                  limit=150,full_output=1)[0]
#        result = quad(lambda E: BringmannSmoothed(E,mchi), Eth, 2.*mchi,
    timeBring_f = time()
#    print 'BringmannIntegralTime= %.1f sec' % (timeBring_f-timeBring_i)
    return result

## def BringIntTest(mchi):
##     timeBring_i = time()
##     result = fixed_quad(lambda E: BringmannSmoothedTest(E,mchi), Eth, 2.*mchi)
##     #,                        vec_func=False)[0]
##     timeBring_f = time()
## #    print 'BringmannIntegralTime= %.1f sec' % (timeBring_f-timeBring_i)
##     return result

# CHB: Testing ...
#print '\n CHB: Testing ... Integration\n'
#print 'BringInt(2000.) = ', BringInt(2000.)
#print 'BergInt(2000.) = ', BergInt(2000.)
#BringIntTest(1000.)
#sys.exit()


# 'old' limit calculation: Now replaced by integration including Aeffs!
def BergSigmavUL(mchi):
    result = 8.*pi*FluxUL*mchi**2 / (Jdeltaomega*BergInt(mchi))
    return result

def BringSigmavUL(mchi):
    result = 8.*pi*FluxUL*mchi**2 / (Jdeltaomega*BringInt(mchi))
    return result


