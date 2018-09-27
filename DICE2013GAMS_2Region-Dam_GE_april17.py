"""

"""
from scipy.stats import lognorm
import numpy as np
import random as rn
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import csv
from numpy import genfromtxt
#import xlsxwriter as xlwt
# =========================================== Global Parameters =========================================== #
#%reset
global al, sigma, etree, l, cost1, rr, forcoth
# =========================================== Parameters =========================================== #
tstep = 5                                               #Years per Period
hzn = 1                      
rng_ADP = int(20/tstep + hzn)                           #number of steps ahead for ADP                                     
T = int(300/tstep)                                      #Modeling time horizon
T0 = T - 10                                             #Graph time horizon
niter = 50                                               #Number of iterations
stepsize = 0.1
zx = 4
## Preferences  ##
elasmu = [1.45, 1.45]                                   #Elasticity of marginal utility of consumption
prstp = [0.015, 0.015]                                  #Initial rate of social time preference per year

## == Population and technology == ##
gama = [0.3, 0.3]                                       #Capital elasticity in production function 
pop0 = [6838.00/2, 6838.00/2]                               #Initial world population (millions)
popadj = [0.134*tstep/5, 0.134*tstep/5]                 #Growth rate to calibrate to 2050 pop projection
popasym = [10500.00/2, 10500.00/2]                          #Asymptotic population (millions)
dk = [0.10, 0.1]                                        #Depreciation rate on capital (per year)
q0 = [63.69/2, 63.69/2]                                     #Initial world gross output (trill 2005 USD)
k0 = [135.00/2, 135.00/2]                                   #Initial capital value (trill 2005 USD)
a0 = [3.80, 3.80]                                       #Initial level of total factor productivity      
ga0 = [0.0790, 0.0790]                                  #Initial growth rate for TFP per 5 years 
dela = [0.0060, 0.0060]                                 #Decline rate of TFP per 5 years

## == Emissions parameters == ##
gsigma1 = [-0.010, -0.010]                              #Initial growth of sigma (per year)
dsig = [-0.0010, -0.0010]                               #Decline rate of decarbonization (per period)
eland0 = [3.3/2, 3.3/2]                                     #Carbon emissions from land 2010 (GtCO2 per year)
deland = [0.2, 0.2]                                     #Decline rate of land emissions (per period)
e0 = [33.61, 33.61]                                     #Industrial emissions 2010 (GtCO2 per year)
miu0 = [0.039, 0.039]                                   #Initial emissions control rate for base case 2010 

## == Carbon cycle == ##
# == Initial Conditions == #
mat0 = 830.4                                   #Initial Concentration in atmosphere 2010 (GtC)
mu0 = 1527.00                                #Initial Concentration in upper strata 2010 (GtC)
ml0 = 10010.00                              #Initial Concentration in lower strata 2010 (GtC)
mateq = 588.00                                #Equilibrium concentration atmosphere  (GtC)
mueq = 1350.00                              #Equilibrium concentration in upper strata (GtC)
mleq = 10000.00                             #Equilibrium concentration in lower strata (GtC)

# == Flow paramaters == #
b12 = 0.088*tstep/5                                     #atmosphere to biosphere/shallow oceans (b12)
b11 = 1 - b12                                           #atmosphere to atmosphere (b11)
b21 = b12 * mateq/ mueq                                 #biosphere/shallow oceans to atmosphere (b21)
b23 = 0.0025*tstep/5                                    #biosphere/shallow oceans to deep oceans (b23)
b22 = 1 - b21 - b23                                     #biosphere/shallow oceans to biosphere/shallow oceans (b22)
b32 = b23 * mueq/ mleq                                  #deep oceans to biosphere/shallow oceans (b32)
b33 = 1 - b32                                           #deep oceans to deep oceans (b33)
sig0 = [0, 0]                               
for i in range(2):
    sig0[i] = e0[i]/(q0[i]*(1 - miu0[i]))               #Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)

## == Climate model parameters == ##
t2xco2 = 2.9                                #Equilibrium temp impact (oC per doubling CO2)
fex0 = 0.25                                 #2010 forcings of non-CO2 GHG (Wm-2) 
fex1 = 0.70                                 #2100 forcings of non-CO2 GHG (Wm-2) 
tocean0 = 0.0068                            #Initial lower stratum temp change (C from 1900)
tatm0 = 0.80                                #Initial atmospheric temp change (C from 1900)

c10 = 0.098                                 #Initial climate equation coefficient for upper level
c1beta = 0.01243                            #Regression slope coefficient(SoA~Equil TSC) 

c1 = 0.098                                  #Climate equation coefficient for upper level
c3 = 0.088                                  #Transfer coefficient upper to lower stratum
c4 = 0.025                                  #CTransfer coefficient for lower level
fco22x = 3.8                                #Forcings of equilibrium CO2 doubling (Wm-2)
c1 =  c10 + c1beta*(t2xco2 - 2.9)

## == Climate damage parameters == ##
a10 = 0                                     #Initial damage intercept
a20 = 0.00267                               #Initial damage quadratic term
a1 = 0                                      #Damage intercept
a2 = 0.00267                                #Damage quadratic term
a3 = 2.00                                   #Damage exponent

## == Abatement cost == ##
expcost2 = [2.8, 2.8]                               #Exponent of control cost function
pback = [344, 344]                                  #Cost of backstop 2005$ per tCO2 2010
gback = [0.025, 0.025]                              #Initial cost decline backstop cost per period
limmiu = [1.2, 1.2]                                 #Upper limit on control rate after 2150
tnopol = [45, 45]                                   #Period before which no emissions controls base
cprice0 = [1.0, 1.0]                                #Initial base carbon price (2005$ per tCO2)
gcprice = [0.02, 0.02]                              #Growth rate of base carbon price per year 

# == Other parameters == # 
lam = fco22x/ t2xco2                                #Climate model parameter
optlrsav = [0, 0]                         #Optimal long-run savings rate used for transversality
for j in range(2):
    optlrsav[j] = (dk[j] + .004)/(dk[j] + .004*elasmu[j] + prstp[j])*gama[j]

# == Geoengineering parameters == # 
Gcost = [0.0027/20, 0.0027/20]
Gexpcost2 = [2, 2]

origdamage = a1 * tatm0 + a2 * tatm0**a3

#   In new specification, 80% of this should be directly from temperature
a2temp = [0.80 * origdamage/tatm0**2, 0.80 * origdamage/tatm0**2]
#   10% should be directly from upper ocean concentration
a2ocean = [0.10 * origdamage/(mu0 - mueq)**2, 0.10 * origdamage/(mu0 - mueq)**2]
#   10% should be directly from atmospheric concentration
a2atmos = [0.10 * origdamage/(mat0 - mateq)**2, 0.10 * origdamage/(mat0 - mateq)**2]

nuG = [0.03/20, 0.03/20]
# ============================================ State Variables ============================================ #

#Capital ($trill, 2005$)
Kt_CM = np.zeros((zx, 2, T))                            #Competition
Kt_CM[:, :, 0] = k0
     
#Atmospheric concentration of carbon (GTC)
MatCM = mat0 + np.zeros((zx, T))                        #Competition

#Concentration in biosphere and upper oceans (GTC)       
MloCM = ml0 + np.zeros((zx, T))                       #Competition

#Concentration in deep oceans (GTC) Cooperation
MupCM = mu0 + np.zeros((zx, T))                       #Competition

#Atmospheric temperature (degrees Celsius above preindustrial)
TatCM = tatm0 + np.zeros((zx, T))                     #Competition

#Lower ocean temperature (degrees Celsius above preindustrial)
TloCM = tocean0 + np.zeros((zx, T))                   #Competition

#Total increase in radiative forcing since preindustrial (Watts per square meter)
RFCM = np.zeros((zx, T))                                #Competition

# ============================================ Exogenous Variables ============================================ #

# == Population == #
#Level of population and labor
l = np.zeros((2, T))                                   

# == Technology == #
#Level of total factor productivity
al = np.zeros((2, T))                                  

#Growth rate of productivity
ga = np.zeros((2, T))

 # == Forcing == #
#Exogenous forcing for other greenhouse gases
forcoth = [1]*(T + 1)

# == Abetment Cost == #
#Adjusted cost for backstop
theta1 = np.zeros((2, T))

#Backstop price
pbacktime = np.zeros((2, T))
 
# == Emissions == #
#Abatement cost function coefficient
sigma = np.zeros((2, T))

#Change in sigma (cumulative improvement of energy efficiency)
gsig = np.zeros((2, T))

#Emissions from deforestation
etree = np.zeros((2, T))

# == Abatement == #
#Adjusted cost for backstop
cost1 = np.zeros((2, T))

#Average utility social discount rate
rr = 1 + np.zeros((2, T))

 # == Results == #
XoptCM = np.zeros((2, 2 * T))
XoldCM = np.zeros((2, 2 * T))
XnewCM = np.zeros((2, 2 * T))                             
erCM = [0]*2*T
      
# ========================================= Exogenous variables ========================================= #
for i in range(T):
    for j in range(2):
        if i == 0:
            al[j, i] = a0[j]
            pbacktime[j, i] = pback[j]
            gsig[j, i] = gsigma1[j]/(1 + dsig[j])**tstep
            sigma[j, i] = sig0[j]
            etree[j, i] = eland0[j]
            l[j, i] = pop0[j]
        else:
            rr[j, i] = rr[j, i-1]/(1 + prstp[j])**tstep
            gsig[j, i] = gsig[j, i-1] * (1 + dsig[j])**tstep
            sigma[j, i] = sigma[j, i-1]*np.exp(gsig[j, i]*tstep)
            al[j, i] = al[j, i-1]/(1 - ga[j, i-1])**(tstep/5)
            pbacktime[j, i] = pbacktime[j, i-1]*(1 - gback[j])**(tstep/5)
            etree[j, i] = etree[j, i-1] * (1 - deland[j])**(tstep/5)
            l[j, i] = l[j, i-1]*(popasym[j]/l[j, i-1])**popadj[j]
        ga[j, i] = ga0[j] * np.exp(-dela[j] * tstep * i)
        cost1[j, i] = pbacktime[j, i] * sigma[j, i]/expcost2[j]/1000
    if i<(90/tstep):
        forcoth[i] = fex0 + (1/18)*(tstep/5) * (fex1 - fex0) * (tstep/5 + i - 1)
    else:
        forcoth[i] = fex1
RFCM[:, 0] = fco22x * (np.log(mat0/mateq))/np.log(2) + forcoth[0]
# ============================================ Output Variables ============================================ #
 # == DICE values == #
 
 #Output gross of abatement cost and climate damage ($trill)
Yt_CM = np.zeros((zx, 2, T))

#Output net of abatement cost and climate damage ($trill)
Qt_CM = np.zeros((zx, 2, T))

#Total carbon emissions (GTCO2 per year)                             
Et_CM = np.zeros((zx, 2, T))

#Industrial emissions (GTCO2 per year)                             
Eind_CM = np.zeros((zx, 2, T))

#World emissions intensity (sigma)                           
Eint_CM = np.zeros((zx, 2, T))

#Total damage (fraction of gross output)                           
Dm_CM = np.zeros((zx, 2, T))

#Climate damages (trillion $)                             
Dt_CM = np.zeros((zx, 2, T))

#Abatement cost (fraction of output)
Bm_CM = np.zeros((zx, 2, T))

#Abatement cost ($ trillion)
Bt_CM = np.zeros((zx, 2, T))

#Geoengineering cost (fraction of output)
Gm_CM = np.zeros((zx, 2, T))

#Geoengineering cost ($ trillion)
Gt_CM = np.zeros((zx, 2, T))

#Saving ($trill, 2005$)
Sv_CM = np.zeros((zx, 2, T)) 

#Consumption ($trill per year)
Cn_CM = np.zeros((zx, 2, T))

#Consumption per capita ($thous per year)
cn_CM = np.zeros((zx, 2, T)) 

#Total period utility
Ut_CM = np.zeros((zx, 2, T))
                             
#Utility of p. c. consumption
ut_CM = np.zeros((zx, 2, T))

#Carbon price (2005$ per ton of CO2)
Cp_CM = np.zeros((zx, 2, T))

#Optimal Abatement rate (DICE)
Abt_CM = np.zeros((zx, 2, T))

#Optimal SGE rate (DICE)
Geo_CM = np.zeros((zx, 2, T))

# ========================================= Transition Function ========================================= #

def state( STi, EX1, EX2, Act1, Act2, Inf1, Inf2 ):
    
    [Matx1, Mupx1, Mlox1, Tatx1, Tlox1, RFCx1] = STi

    [Ki1, A1, sig1, Eland1, L1, cost11, df1, pb1, Fexxi1, Fexxii1] = EX1
    [Ki2, A2, sig2, Eland2, L2, cost12, df2, pb2, Fexxi2, Fexxii2] = EX2
    
    [xa1, xg1] = Act1
    [xa2, xg2] = Act2
    
    [a2tempx1, a2oceanx1, a2atmosx1, nuGx1] = Inf1
    [a2tempx2, a2oceanx2, a2atmosx2, nuGx2] = Inf2
    
    RFCx1g = (fco22x * (np.log(Matx1/mateq))/np.log(2) + (Fexxi1 + Fexxi2)) - (xg1 + xg2)
    Tatx1 = Tatx1 + c1 * (RFCx1g - RFCx1)
    
    Yi1 = A1 * (L1/1000)**(1 - gama[0]) * Ki1**gama[0]
    Yi2 = A2 * (L2/1000)**(1 - gama[1]) * Ki2**gama[1]
    
    Eind1 = sig1 * Yi1 * (1 - xa1)
    Eind2 = sig2 * Yi2 * (1 - xa2)
    
    E1 = Eind1 + Eland1
    E2 = Eind2 + Eland2
    E0 = E1 + E2
    
    Eint1 = E1/Yi1
    Eint2 = E2/Yi2
    
#    Dm1 = (a1 * Tatx1 + a2 * Tatx1**a3) + (nuGx1 * xg1**2)
#    Dm2 = (a1 * Tatx1 + a2 * Tatx1**a3) + (nuGx2 * xg2**2)
    Dm1 = (a2tempx1 * Tatx1**a3 + a2oceanx1 * (Mupx1 - mueq)**2 + a2atmosx1 * (Matx1 - mateq)**2) + (nuGx1 * (xg1 + xg2)**2)
    Dm2 = (a2tempx2 * Tatx1**a3 + a2oceanx2 * (Mupx1 - mueq)**2 + a2atmosx2 * (Matx1 - mateq)**2) + (nuGx2 * (xg1 + xg2)**2)
    
    D1 = Yi1 * Dm1
    D2 = Yi2 * Dm2
    
    Bm1 = cost11 * xa1**expcost2[0]
    Bm2 = cost12 * xa2**expcost2[1]
    
    B1 = Bm1 * Yi1
    B2 = Bm2 * Yi2
    
    Gm1 = Gcost[0] * xg1**Gexpcost2[0]
    Gm2 = Gcost[1] * xg2**Gexpcost2[1]
    
    G1 = Gm1 * Yi1
    G2 = Gm2 * Yi2
    
    Q1 = Yi1 - D1 - B1 - G1
    Q2 = Yi2 - D2 - B2 - G2
    
    S1 = optlrsav[0] * Q1
    S2 = optlrsav[1] * Q2
    
    Con1 = Q1 - S1
    Con2 = Q2 - S2
    
    con1 = (Con1/L1) * 1000
    con2 = (Con2/L2) * 1000
    
    u1 = (con1**(1 - elasmu[0]) - 1)/ (1 - elasmu[0]) - 1
    u2 = (con2**(1 - elasmu[1]) - 1)/ (1 - elasmu[1]) - 1
    
    U1 = u1 * L1 * df1
    U2 = u2 * L2 * df2
    
    Cp1   = pb1 * (xa1)**(expcost2[0]-1)
    Cp2   = pb2 * (xa2)**(expcost2[1]-1)
    
    Matx2 = b11 * Matx1 + b21 * Mupx1 + (E0) * tstep/3.666
    Mupx2 = b12 * Matx1 + b22 * Mupx1 + b32 * Mlox1
    Mlox2 = b23 * Mupx1 + b33 * Mlox1

    RFCx2 = fco22x * (np.log(Matx2/mateq))/np.log(2) + (Fexxii1 + Fexxii2)
    Tatx2 = Tatx1 + c1 * (RFCx2 - (fco22x/t2xco2) * Tatx1 - c3 * (Tatx1 - Tlox1))
    Tlox2 = Tlox1 + c4 * (Tatx1 - Tlox1)
  
    Kii1 = (1 - dk[0])**tstep * Ki1 + tstep * S1
    Kii2 = (1 - dk[1])**tstep * Ki2 + tstep * S2

    STii = [Matx2, Mupx2, Mlox2, Tatx2, Tlox2, RFCx2]
    LEVi1 = [Kii1, Yi1, Q1, E1, Eind1, Eint1, Dm1, D1, Bm1, B1, Gm1, G1, S1, Con1, con1, U1, u1, Cp1]
    LEVi2 = [Kii2, Yi2, Q2, E2, Eind2, Eint2, Dm2, D2, Bm2, B2, Gm2, G2, S2, Con2, con2, U2, u2, Cp2]
    STGEi = [Tatx1, RFCx1g]
    
    return ( STii, LEVi1, LEVi2, STGEi )

# =========================================== Welfare Function =========================================== #

def SGE(v):
    W = 0
    for t in range(T-1):
        vxa = v[t]
        vxg = v[t + T]
        if region == 1:
            vxa1 = vxa
            vxg1 = vxg
            vxa2 = XnewCM[1, t]
            vxg2 = XnewCM[1, t + T]
        else:
            vxa1 = XnewCM[0, t]
            vxg1 = XnewCM[0, t + T]
            vxa2 = vxa
            vxg2 = vxg
            
        vSTi = [MatCM[stp, t], MupCM[stp, t], MloCM[stp, t], TatCM[stp, t], TloCM[stp, t], RFCM[stp, t]]
        vEX1 = [Kt_CM[stp, 0, t], al[0, t], sigma[0, t], etree[0, t], l[0, t], cost1[0, t], rr[0, t], pbacktime[0, t], forcoth[t]/2, forcoth[t + 1]/2]
        vEX2 = [Kt_CM[stp, 1, t], al[1, t], sigma[1, t], etree[1, t], l[1, t], cost1[1, t], rr[1, t], pbacktime[1, t], forcoth[t]/2, forcoth[t + 1]/2]
        vAct1 = [vxa1, vxg1]
        vAct2 = [vxa2, vxg2]
        vInf1 = [a2temp[0], a2ocean[0], a2atmos[0], nuG[0]]
        vInf2 = [a2temp[1], a2ocean[1], a2atmos[1], nuG[1]]
        
        ( vSTii, vLEVi1, vLEVi2, vSTGEi ) = state(vSTi, vEX1, vEX2, vAct1, vAct2, vInf1, vInf2)
        
        [MatCM[stp, t + 1], MupCM[stp, t + 1], MloCM[stp, t + 1], TatCM[stp, t + 1], TloCM[stp, t + 1], RFCM[stp, t + 1]] = vSTii
        [Kt_CM[stp, 0, t+1], Yt_CM[stp, 0, t], Qt_CM[stp, 0, t], Et_CM[stp, 0, t], Eind_CM[stp, 0, t], Eint_CM[stp, 0, t], Dm_CM[stp, 0, t], Dt_CM[stp, 0, t], Bm_CM[stp, 0, t], Bt_CM[stp, 0, t], Gm_CM[stp, 0, t], Gt_CM[stp, 0, t], Sv_CM[stp, 0, t], Cn_CM[stp, 0, t], cn_CM[stp, 0, t], Ut_CM[stp, 0, t], ut_CM[stp, 0, t], Cp_CM[stp, 0, t]] = vLEVi1
        [Kt_CM[stp, 1, t+1], Yt_CM[stp, 1, t], Qt_CM[stp, 1, t], Et_CM[stp, 1, t], Eind_CM[stp, 1, t], Eint_CM[stp, 1, t], Dm_CM[stp, 1, t], Dt_CM[stp, 1, t], Bm_CM[stp, 1, t], Bt_CM[stp, 1, t], Gm_CM[stp, 1, t], Gt_CM[stp, 1, t], Sv_CM[stp, 1, t], Cn_CM[stp, 1, t], cn_CM[stp, 1, t], Ut_CM[stp, 1, t], ut_CM[stp, 1, t], Cp_CM[stp, 1, t]] = vLEVi2
        [TatCM[stp, t], RFCM[stp, t]] = vSTGEi
        
        vU1 = Ut_CM[stp, 0, t]
        vU2 = Ut_CM[stp, 1, t]
        if region == 1:
            vU = vU1 + stp/3 * vU2
        else:
            vU = vU2 + stp/3 * vU1               
        W = W +  vU * tstep
    return -W

# ======================================== Optimization Algorithm ======================================== #
for stp in range(0, zx, 1):
    
    # == bounds == #
    bnds = 2 * T * [(0, 1.0)]
    bnds[0] = (0.039, 0.039)
    bnds[T:2*T] = T*[(0.0, 5.0)]
    bnds[T] = (0.0, 0.0)
    
    for i in range(niter):

        # == optimization == #
        ftol = 1e-12
        eps = 1e-6
        maxiter = 10000
    
        # == initial guess == #
        x0 = XnewCM[1, :]
        #2 * (T) * [0.01]
    
        for region in [1, 2]:
            res = minimize(SGE, x0, method='SLSQP', bounds=bnds, options={'ftol': ftol, 'eps': eps, 'disp': True, 'maxiter': maxiter})
            # == results == #
            result = res.fun
            x = res.x
            XoptCM[(region - 1), 0:T]= x[0:T]
            XoptCM[(region - 1), T:2*T] = x[T:2*T]
        
        for t in range(2*T):
            if i == 0:
                XnewCM[:, t] = XoptCM[:, t]
            XoldCM[:, t] = XnewCM[:, t]
            XnewCM[:, t] = XoptCM[:, t] * stepsize + XoldCM[:, t] * (1 - stepsize)
            erCM[t] = sum(abs(XoldCM[:, t] - XoptCM[:, t])) 

    for j in range(2):
        Abt_CM[stp, j, 0:T] = XnewCM[j, 0:T]
        Geo_CM[stp, j, 0:T] = XnewCM[j, T:2*T]

# ======================================== Comparative Analysis ======================================== #

T0 = T - 10
Xaxis = range(2010, 2010 + tstep * T0, tstep)

plt.plot(Xaxis[0:T0], Et_CM[0, 0, 0:T0], label = "Competition")
plt.plot(Xaxis[0:T0], Et_CM[1, 0, 0:T0], label = "Coordination-Competition")
plt.plot(Xaxis[0:T0], Et_CM[2, 0, 0:T0], label = "Coordination-Cooperation")
plt.plot(Xaxis[0:T0], Et_CM[3, 0, 0:T0], label = "Cooperation")
plt.xlabel('Time (years)')
plt.ylabel('Emissions (GtC)')
plt.title('Optimal Emissions')
plt.legend(loc=1, prop={'size':8})
plt.ylim(ymax = 200, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 5 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], Geo_CM[0, 0, 0:T0], label = "Competition")
plt.plot(Xaxis[0:T0], Geo_CM[1, 0, 0:T0], label = "Coordination-Competition")
plt.plot(Xaxis[0:T0], Geo_CM[2, 0, 0:T0], label = "Coordination-Cooperation")
plt.plot(Xaxis[0:T0], Geo_CM[3, 0, 0:T0], label = "Cooperation")
plt.xlabel('Time (years)')
plt.ylabel('Geoengineering (W/m2)')
plt.title('Optimal Geoengineering')
plt.legend(loc=1, prop={'size':8})
plt.ylim(ymax = 5, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 5 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], TatCM[0, 0:T0], label = "Competition")
plt.plot(Xaxis[0:T0], TatCM[1, 0:T0], label = "Coordination-Competition")
plt.plot(Xaxis[0:T0], TatCM[2, 0:T0], label = "Coordination-Cooperation")
plt.plot(Xaxis[0:T0], TatCM[3, 0:T0], label = "Cooperation")
plt.xlabel('Time (years)')
plt.ylabel('Atmospheric temperature (degree C)')
plt.title('Temperature')
plt.legend(loc=4, prop={'size':8})
plt.ylim(ymax = 4, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 5 * tstep))
plt.show()