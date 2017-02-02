## Continue in $I_\mathrm{dc}$ to illustrate frequency dependence ##
## Simulation must be run first (p["simfile"])                    ##
####################################################################
import matplotlib

import pdb
#matplotlib.use("Agg")
import pylab, brian2, brianutils, os, json, sympy, scipy, sys, \
    datetime, shelve, autoutils, auto, contextlib
from sympy import S
units= dict(brian2.units.__dict__.items()
            + brian2.units.allunits.__dict__.items()
            + brian2.__dict__.items())
unitlist=["mV","ms","cm2","uF","psiemens","um2","msiemens","cm"]
baseunits=[(k,float(eval(k,units).base)) for k in unitlist]



######################
# PARAMETER SETTINGS #
p = {"simfile"  : "sim_DNap_2015-12-14_16:55:14.603625.shv",
     "modfile"  : "cfg/wb01.json",      #This is our model
     "contfile" : "dat/cont_WB.json",   #This is what brian produces and hands to auto
     "workdir"  : os.getenv("HOME") + "/Uni/CNS/3.Semester/schreiberlab/project/sim1",
     "dt"       : "0.05*ms",
     "bifpar"   : {
       "I" : ["0.*uA/cm2"],
       "Cm" : ["1.0*uF/cm2","1.1*uF/cm2","1.4*uF/cm2","1.41*uF/cm2","1.42*uF/cm2"]
       }
     }

os.chdir(p["workdir"])
#Now, we define how to load in the model in a nice pythonic way and sort them
#in a way that brian will understand them.
def load_mod(modfile,bp):
    nakdic = json.load(open(modfile))
    fundic = dict([(j,k.split(":")[0]) for j,k in nakdic["fun"].items()])
    pardic = dict([(j,k) for j,k in nakdic["par"].items()])
    bifpar = [(k,pardic.pop(k)) for k in bp]
    sdelist = [[j,] + k.split(':') for j,k in nakdic['aux_odes'].items()]
    sdelist = [(i,":".join((str(S(j).subs(fundic).subs(fundic).subs(pardic)),k))) for i,j,k in sdelist]
    sdelist+= [('v', str(sympy.solve(nakdic['current_balance_eq'],'dv/dt')[0].subs(nakdic['currents']).subs(fundic).subs(fundic).subs(pardic))+":volt")]
    sde = brian2.Equations("d{}/dt = {}".format(*sdelist[0]))
    for i,j in sdelist[1:]:
        sde += brian2.Equations("d{}/dt = {}".format(i,j))
    return sde

## FIND INITIAL STEADY STATE ##
sde = load_mod(p["modfile"],p['bifpar'].keys())
ode = brianutils.sde2ode(sde)
diffuterms=dict([(k[3:],(S(j).coeff(k)**2).subs(baseunits)) for i,j in sde.eq_expressions for k in sde.stochastic_variables if S(j).coeff(k)!=0])

## ADD PAR ##
for j,k in p["bifpar"].items():
  ode += brian2.Equations("{} : {}".format(j,repr(eval(k[0],units).dim)))

brian2.defaultclock.dt = eval(p["dt"], units)
G = brian2.NeuronGroup(1, model=ode, method="rk4",
                       threshold='not_refractory and (v>5*mV)',
                       refractory='v>-40*mV')

# PAR INIT #
for j,k in p["bifpar"].items():
  setattr(G,j,eval(k[0],units))

# STATE INIT #
G.v= eval("-66 * mV", units)

states = brian2.StateMonitor(G, ode.eq_names, record=True)
spikes = brian2.SpikeMonitor(G)
net = brian2.Network(G,states,spikes)
duration = eval("500 * ms",units)
net.run(duration)

autobifpar = dict([(i,float(eval(j[0],units))) for i,j in p['bifpar'].items()])


## CREATE ADJOINT LINEAR SYSTEM ##
varrhs = [(i,sympy.S(j).subs(baseunits))
                for i,j in ode.eq_expressions]
varrhs.sort(cmp=lambda x,y:cmp(x[0],y[0]),reverse=True)
var,rhs = zip(*varrhs);
advar = sympy.S(["ad{}".format(k) for k in var])
J = [[S(i).diff(j) for j in var] for i in rhs]
J = [[j.subs(baseunits) for j in k] for k in J]
adlinsys = [str(k) for k in
            (sympy.S("lam")*sympy.eye(len(advar))-sympy.Matrix(J).T)*sympy.Matrix(advar)]
prcnorm=str((sympy.Matrix(sympy.S(advar)).T*sympy.Matrix(sympy.S(rhs)))[0,0] - sympy.S("dotZF/period"))

# Ipar = eval("1.2 * uA/cm2", units) # LC
spikecriterion = [str(S(k).subs([(i,"{}_left".format(i)) for i in var]))
                  for j,k in zip(var,rhs) if j=="v"]

#Hier kommt das Pythondings von Jan-Hendrik zum Einsatz!
if "A" in autobifpar:
    autobifpar.pop("A")

unames,pnames= autoutils.writeFP('tm',
    bifpar=autobifpar, rhs=rhs, var=var,
    bc=['{0}_left-{0}_right'.format(v) for v in var] + spikecriterion,
    ic=[])

#pdb.set_trace()
################
# CONT FP & LC #
r1= auto.run([float(getattr(states,j)[0][-1]) for j in var], e='tm',
    c='tm', parnames= pnames, unames=unames,
    ICP=['I', 'Cm'], ISP=1,ILP=1, SP=['LP','HB','BP'],
    PAR=autobifpar, ITNW=17, NWTN=13, NMX=500000, NPR=500000,
    DS=1e-6, DSMAX=1e-5, STOP=['HB1'],
    UZSTOP= {'I': 350.0})

s1HB = r1.getLabel('HB')[0]
s1LP = r1.getLabel('LP')[0]
pdb.set_trace()


##Now: continuate HB in I/Cm - direction.
direction = 1
runBranch = auto.run(s1HB, e='tm', c='tm',
                 parnames=pnames, unames=unames,
                 ICP=['I'],
                 ISP=2,ILP=1, SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
                 ISW=1, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
                 ITNW=17, NWTN=13, NMX=10000, NPR=2000, #Not sure if needed: PAR=autobifpar,
                 DS=direction*1e-6, DSMAX=1e-5,
                 UZSTOP= {'I':350.0}
                 )
pdb.set_trace()






####Der Quatsch wird erstmal rausgelassen

r2= auto.run(s1HB, e='tm', c='tm',
    parnames= pnames, unames= unames,
    ICP=['I','period'], ILP=1, ISW=1,IPS=2,
    ITNW=7, NWTN=3, NMX=1000, NPR=1000,
    DS=-1e-2, DSMAX=1e1, NTST= 300,
    SP=['BT','LP','HB','BP','CP'],
    UZSTOP={'period':0.15})
pdb.set_trace()
# Jacobian at the saddle-node
statepardic = dict(zip(s1LP.coordnames,s1LP.coordarray.flatten()))
statepardic.update(s1LP.PAR)
J0 = [[S(k).subs(statepardic) for k in j] for j in J]

# fig = pylab.figure(1,figsize=(3,3))
# fig.clf()
# ax = fig.add_subplot(111)
# ax.plot(r1["I"],r1["v"],color="k",lw=1)
# ax.plot(r2["I"],r2["MAX v"],color="k",lw=2)
# fig.savefig("bifdiag.svg")





##############
## CONT BVP ##
sol = r2.getLabel('UZ')[0] # STARTING POINT

ix = sol['v'].argmax()
orbdat = sol.coordarray[:,:-1]
orbdat = scipy.roll(orbdat,-ix,1)
tdat = sol.indepvararray[:-1]
tdat = scipy.mod(scipy.roll(tdat,-ix) - tdat[ix],1)
dat= scipy.zeros((2*min(orbdat.shape)+1,max(orbdat.shape)))
dat[0,:] = tdat
dat[1:min(orbdat.shape)+1,:]= orbdat
dat[min(orbdat.shape)+1:,:] = scipy.zeros((min(orbdat.shape),max(orbdat.shape)))
dat[min(orbdat.shape)+1:,:] = -10.*(1-scipy.cos(2*scipy.pi*dat[0,:])[None,:].repeat(min(orbdat.shape),0))
goldstone= scipy.array([scipy.gradient(k)/scipy.gradient(dat[0,:]) for k in dat[1:min(orbdat.shape)+1,:]])
dotZF= scipy.trapz((goldstone*dat[min(orbdat.shape)+1:,:]).sum(0),dat[0,:])

autobifpar['period'] = sol['period']
autobifpar['lam'] = 0
autobifpar['dotZF'] = dotZF
autobifpar.update([(i,j) for i,j in sol.PAR.items()
                   if i in autobifpar.keys()])

unames, pnames= autoutils.writeBVP('tmBVP',
    bifpar=autobifpar,
    rhs=list(rhs)+adlinsys,
    var=list(var)+advar,
    bc=['{0}_left-{0}_right'.format(v) for v in list(var)+advar]\
    +spikecriterion,
    ic=[prcnorm])

r3= auto.run(dat, e='tmBVP', c='tmBVP',
    parnames= pnames, unames=unames, PAR=autobifpar,
    ICP=['dotZF','lam','period','I'],
    NTST= 200, ITNW=17, NWTN=13, NMX=200, NPR=200,
    DS=1e-4, DSMAX=1e1, UZSTOP={'dotZF':1})
pdb.set_trace()

r4= auto.run(r3.getLabel('UZ')[0], e='tmBVP', c='tmBVP',
    parnames= pnames, unames=unames, ICP=['period','I','lam'],
    NTST=200, ITNW=7, NWTN=3, NMX=2000, NPR=2000, DS=5e-2, DSMAX=1e1,
    UZSTOP={'period':.5})
pdb.set_trace()

Cm_list = [float(eval(k,units)) for k in p['bifpar']['Cm']]

Cm_list = [.01,.011,.014,.0143,.0145,.0146,.01465]

r5= auto.run(r4.getLabel("UZ")[0], e='tmBVP', c='tmBVP',
       parnames= pnames, unames=unames, ICP=['Cm','I','lam'],
       NTST=200, NCOL=2, ITNW=13, NWTN=10, NMX=8000,
       NPR=8000, DS=1e2, DSMAX=1e12, UZR={"Cm":Cm_list[:-1]},
       UZSTOP={"Cm":Cm_list[-1]},THL={"I":0,"lam":0})
pdb.set_trace()

with open(p['contfile'],'w') as f:
  d = {}
  d.update(p)
  sols = r5.getLabel('UZ')
  d['sol'] = [k.toArray() for k in sols]
  d['par'] = [k.PAR.todict() for k in sols]
  d['t'] = [k.indepvararray.tolist() for k in sols]
  d['coordnames'] = k.coordnames
  # d['J0'] = J0
  d['saddle'] = {}
  d['saddle']['coordarray'] = s1LP.toArray()
  d['saddle']['coordnames'] = s1LP.coordnames
  d['saddle']['par'] = s1LP.PAR.todict()
  d['diffu'] = [(j,str(k)) for j,k in diffuterms.items()]
  json.dump(d,f)
