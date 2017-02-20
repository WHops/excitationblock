## Continue in $I_\mathrm{dc}$ to illustrate frequency dependence ##
## Simulation must be run first (p["simfile"])                    ##
####################################################################
import matplotlib

import pdb
#matplotlib.use("Agg")
import pylab, brian2, brianutils, os, json, sympy, scipy, sys, \
    datetime, shelve, autoutils, auto, contextlib, inibrian
from sympy import S


########### Get initial condition through brian ###########

unames, pnames, ini, autobifpar = inibrian.find_ini_brian()
autobifpar['Cm'] = 1


################# CONT FP & LC ############################


calc_period = 0 #do we want to calculate the orbits?
refine_hb = 0 #do we want to get the HB position even more exactly? (Not needed and also not really working :))
directionlist = [1, -1]

r1_fwd= auto.run(ini, e='tm_new',
    c='tm_new', parnames= pnames, unames=unames,
    ICP=['I','Cm'], ISP=1,ILP=1, SP=['LP','HB','BP'],
    PAR=autobifpar, ITNW=17, NWTN=13, NMX = 200000,NPR=200000, #NMX=113500, NPR=5000,
    DS= directionlist[0] * 1e-3, DSMAX= directionlist[0] * 1e-2, STOP=['HB1'],
    UZSTOP= {})

r1_bwd= auto.run(ini, e='tm_new',
    c='tm_new', parnames= pnames, unames=unames,
    ICP=['I','Cm'], ISP=1,ILP=1, SP=['LP','HB','BP'],
    PAR=autobifpar, ITNW=17, NWTN=13, NMX = 3000, NPR=3000,#NMX=113500, NPR=5000,
    DS= directionlist[1] * 1e-3, DSMAX= directionlist[1] * 1e-2, #STOP=['HB1'],
    UZSTOP= {})

r1 = r1_fwd + r1_bwd

s1HB = r1_fwd.getLabel('HB')[0]
s1LP = r1.getLabel('LP')[0]
s1UZ = r1.getLabel('EP')[0]

################# Optional: continue orbit #######################

if calc_period:
    r1_period = auto.run(s1HB, e='tm_new', c='tm_new',
        parnames= pnames, unames= unames,
        ICP=['I','period'], ILP=1, ISW=1,IPS=2,
        ITNW=7, NWTN=3, NMX=1000, NPR=1000,
        DS=1e-3, DSMAX=1e-2, NTST= 300,
        SP=['BT','LP','HB','BP','CP'],
        UZSTOP={})
    r1_withperiod = r1 + r1_period
    s1HB = r1_withperiod.getLabel('HB')[0]

################# Optional: refine HB ############################

if refine_hb:
    try:
        r21= auto.run(s1HB, e='tm_new',
            c='tm_new', parnames= pnames, unames=unames,
            ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=10000, NPR=10000,
            DS=-1e-5, DSMAX=-1e-4, STOP=['HB1'],
            IID = 3)
        r22= auto.run(s1HB, e='tm_new',
            c='tm_new', parnames= pnames, unames=unames,
            ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=10000, NPR=10000,
            DS=1e-5, DSMAX=1e-4, STOP=['HB1'],
            IID = 3)
        r2 = r21 + r22
        s2HB = r2.getLabel('HB')[0]
        r31= auto.run(s2HB, e='tm_new',
            c='tm_new', parnames= pnames, unames=unames,
            ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=10000, NPR=10000,
            DS=-1e-5, DSMAX=-1e-4, STOP=['HB1'],
            IID = 3)
        r32= auto.run(s2HB, e='tm_new',
            c='tm_new', parnames= pnames, unames=unames,
            ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=10000, NPR=10000,
            DS=1e-5, DSMAX=1e-4, STOP=['HB1'],
            IID = 3)
        r3 = r31 + r32
        s3HB = r3.getLabel('HB')[0]
        s1HB = s3HB
    except:
        print("Could not refine position: Bifurcation point not found. Continuing with unrefined HB.")

################# CONT HB ############################


runHB = [0 for item in directionlist]
for counter, direction in enumerate(directionlist):
    runHB[counter] = auto.run(s1HB, e='tm_new', c='tm_new',
            parnames=pnames, unames=unames,
            ICP=['I','Cm'], #IPS = 1
            ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
            ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
            ITNW=17, NWTN=13, NMX=10000, NPR=1000, #Not sure if needed: PAR=autobifpar,
            DS= direction * 1e-2, DSMAX= direction * 1e-1, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
            UZSTOP= {}
            )

runHB = runHB[0] + runHB[1]
rhb = runHB + r1
auto.save(rhb, 'rhb')
