## Continue in $I_\mathrm{dc}$ to illustrate frequency dependence ##
## Simulation must be run first (p["simfile"])                    ##
####################################################################
import matplotlib

import pdb
#matplotlib.use("Agg")
import pylab, brian2, brianutils, os, json, sympy, scipy, sys, \
    datetime, shelve, autoutils, auto, contextlib, inibrian
from sympy import S


################# Get initial condition ####################

unames, pnames, ini, autobifpar = inibrian.find_ini_brian()
autobifpar['Cm'] = 1


################# CONT FP & LC ############################

calc_period = 0
refine_hb = 0
direction = [1, -1]

r1_fwd= auto.run(ini, e='tm_new',
    c='tm_new', parnames= pnames, unames=unames,
    ICP=['I'], ISP=1,ILP=1, SP=['LP','HB','BP'],
    PAR=autobifpar, ITNW=17, NWTN=13, NMX = 200000,NPR=200000, #NMX=113500, NPR=5000,
    DS= direction[0] * 1e-3, DSMAX= direction[0] * 1e-2, STOP=['HB1'],
    UZSTOP= {})

r1_bwd= auto.run(ini, e='tm_new',
    c='tm_new', parnames= pnames, unames=unames,
    ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'],
    PAR=autobifpar, ITNW=17, NWTN=13, NMX = 3000, NPR=3000,#NMX=113500, NPR=5000,
    DS= direction[1] * 1e-3, DSMAX= direction[1] * 1e-2, #STOP=['HB1'],
    IID = 3,
    UZSTOP= {})

r1 = r1_fwd + r1_bwd

s1HB = r1_fwd.getLabel('HB')[0]
s1LP = r1.getLabel('LP')[0]
s1UZ = r1.getLabel('EP')[0]

if calc_period:
    r1_period = auto.run(s1HB, e='tm_new', c='tm_new',
        parnames= pnames, unames= unames,
        ICP=['I','period'], ILP=1, ISW=1,IPS=2,
        ITNW=7, NWTN=3, NMX=1000, NPR=1000,
        DS=1e-2, DSMAX=1e1, NTST= 300,
        SP=['BT','LP','HB','BP','CP'],
        UZSTOP={'period':0.15})
    r1_withperiod = r1 + r1_period
    s1HB = r1_withperiod.getLabel('HB')[0]

#--------------Optional: refine the HB ------------------------

if refine_hb:
    r2= auto.run(s1HB_final, e='tm_new',
        c='tm_new', parnames= pnames, unames=unames,
        ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=1000000, NPR=1000000,
        DS=-1e-6, DSMAX=1e-5, STOP=['HB1'],
        IID = 3)
    s2HB = r2.getLabel('HB')[0]

    r3= auto.run(s2HB, e='tm_new',
        c='tm_new', parnames= pnames, unames=unames,
        ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=1000000, NPR=1000000,
        DS=1e-8, DSMAX=1e-7, STOP=['HB1'],
        IID = 3)
    s3HB = r3.getLabel('HB')[0]
    s1HB = s3HB


pdb.set_trace()
################# CONT HB ############################

runHB = [0 for item in direction]
for counter, dir_ in enumerate(direction):
    runHB[1] = auto.run(s1HB, e='tm_new', c='tm_new',
            parnames=pnames, unames=unames,
            ICP=['I','Cm'], #IPS = 1
            ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
            ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
            ITNW=17, NWTN=13, NMX=3000, NPR=3000, #Not sure if needed: PAR=autobifpar,
            DS= dir_ * 1e-3, DSMAX= dir_ * 1e-3, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
            UZSTOP= {}
            )

runHB = runHB[0] + runHB[1]
rhb = runHB + r1
auto.save(rhb, 'rhb')
pdb.set_trace()



#
# rb1EP = runBranch1.getLabel('EP')[0]
#
#
# runBranch2 = auto.run(rb1EP, e='tm_new', c='tm_new',
#                  parnames=pnames_wrap, unames=unames_wrap, TY='HB',
#                  ICP=['I','period'], #IPS = 1
#                  ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
#                  ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
#                  ITNW=17, NWTN=13, NMX=100, NPR=1, #Not sure if needed: PAR=autobifpar,
#                  DS=direction*1e-1, DSMAX=1e-1, IID = 5, #RL1=0.58, EPSL = 1e-1, EPSU = 1e-1, EPSS = 10, NTST = 500,
#                  UZSTOP= {}
#                  )
#
#
#
#
#
# ####Der Quatsch wird erstmal rausgelassen
# import numpy as np
#
#
# r2= auto.run(s1HB, e='tm_new', c='tm_new',
#     parnames= pnames_wrap, unames= unames_wrap,
#     ICP=['Cm','I'], ILP=1, ISW=2,IPS=1,
#     ITNW=17, NWTN=13, NMX=10, NPR=1,
#     DS=stepsize, DSMAX=stepsize, NTST= 3,
#     SP=['BT','LP','HB','BP','CP'],
#     UZSTOP={'period':0.15},
#     EPSL=1e-3, EPSU=1e-3, EPSS=1e-2, ITMX = 20)
#
#
# period_015 = r2.getLabel('UZ')[0]
# r3= auto.run(period_015, e='tm', c='tm',
#     parnames= pnames_wrap, unames= unames_wrap,
#     ICP=['I','period'], ILP=1, ISW=1,IPS=2,
#     ITNW=7, NWTN=3, NMX=1000, NPR=1,
#     DS=1e-2, DSMAX=1e1, NTST= 300,
#     SP=['BT','LP','HB','BP','CP'],
#     UZSTOP={})
#
# per_all = r1 + r3
# auto.save(per_all,'per_all')
# pdb.set_trace()
# # Jacobian at the saddle-node
# statepardic = dict(zip(s1LP.coordnames,s1LP.coordarray.flatten()))
# statepardic.update(s1LP.PAR)
# J0 = [[S(k).subs(statepardic) for k in j] for j in J]
#
# # fig = pylab.figure(1,figsize=(3,3))
# # fig.clf()
# # ax = fig.add_subplot(111)
# # ax.plot(r1["I"],r1["v"],color="k",lw=1)
# # ax.plot(r2["I"],r2["MAX v"],color="k",lw=2)
# # fig.savefig("bifdiag.svg")
#
#
#
#
#
# ##############
# ## CONT BVP ##
# sol = r2.getLabel('UZ')[0] # STARTING POINT
#
# ix = sol['v'].argmax()
# orbdat = sol.coordarray[:,:-1]
# orbdat = scipy.roll(orbdat,-ix,1)
# tdat = sol.indepvararray[:-1]
# tdat = scipy.mod(scipy.roll(tdat,-ix) - tdat[ix],1)
# dat= scipy.zeros((2*min(orbdat.shape)+1,max(orbdat.shape)))
# dat[0,:] = tdat
# dat[1:min(orbdat.shape)+1,:]= orbdat
# dat[min(orbdat.shape)+1:,:] = scipy.zeros((min(orbdat.shape),max(orbdat.shape)))
# dat[min(orbdat.shape)+1:,:] = -10.*(1-scipy.cos(2*scipy.pi*dat[0,:])[None,:].repeat(min(orbdat.shape),0))
# goldstone= scipy.array([scipy.gradient(k)/scipy.gradient(dat[0,:]) for k in dat[1:min(orbdat.shape)+1,:]])
# dotZF= scipy.trapz((goldstone*dat[min(orbdat.shape)+1:,:]).sum(0),dat[0,:])
#
# autobifpar['period'] = sol['period']
# autobifpar['lam'] = 0
# autobifpar['dotZF'] = dotZF
# autobifpar.update([(i,j) for i,j in sol.PAR.items()
#                    if i in autobifpar.keys()])
#
# unames, pnames= autoutils.writeBVP('tmBVP',
#     bifpar=autobifpar,
#     rhs=list(rhs)+adlinsys,
#     var=list(var)+advar,
#     bc=['{0}_left-{0}_right'.format(v) for v in list(var)+advar]\
#     +spikecriterion,
#     ic=[prcnorm])
#
# r3= auto.run(dat, e='tmBVP', c='tmBVP',
#     parnames= pnames, unames=unames, PAR=autobifpar,
#     ICP=['dotZF','lam','period','I'],
#     NTST= 200, ITNW=17, NWTN=13, NMX=200, NPR=200,
#     DS=1e-4, DSMAX=1e1, UZSTOP={'dotZF':1})
# pdb.set_trace()
#
# r4= auto.run(r3.getLabel('UZ')[0], e='tmBVP', c='tmBVP',
#     parnames= pnames, unames=unames, ICP=['period','I','lam'],
#     NTST=200, ITNW=7, NWTN=3, NMX=2000, NPR=2000, DS=5e-2, DSMAX=1e1,
#     UZSTOP={'period':.5})
# pdb.set_trace()
#
# Cm_list = [float(eval(k,units)) for k in p['bifpar']['Cm']]
#
# Cm_list = [.01,.011,.014,.0143,.0145,.0146,.01465]
#
# r5= auto.run(r4.getLabel("UZ")[0], e='tmBVP', c='tmBVP',
#        parnames= pnames, unames=unames, ICP=['Cm','I','lam'],
#        NTST=200, NCOL=2, ITNW=13, NWTN=10, NMX=8000,
#        NPR=8000, DS=1e2, DSMAX=1e12, UZR={"Cm":Cm_list[:-1]},
#        UZSTOP={"Cm":Cm_list[-1]},THL={"I":0,"lam":0})
# pdb.set_trace()
#
# with open(p['contfile'],'w') as f:
#   d = {}
#   d.update(p)
#   sols = r5.getLabel('UZ')
#   d['sol'] = [k.toArray() for k in sols]
#   d['par'] = [k.PAR.todict() for k in sols]
#   d['t'] = [k.indepvararray.tolist() for k in sols]
#   d['coordnames'] = k.coordnames
#   # d['J0'] = J0
#   d['saddle'] = {}
#   d['saddle']['coordarray'] = s1LP.toArray()
#   d['saddle']['coordnames'] = s1LP.coordnames
#   d['saddle']['par'] = s1LP.PAR.todict()
#   d['diffu'] = [(j,str(k)) for j,k in diffuterms.items()]
#   json.dump(d,f)
