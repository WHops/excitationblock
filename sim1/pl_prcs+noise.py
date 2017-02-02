import matplotlib 
matplotlib.use("Agg")
import pylab, shelve, scipy, sympy, brian2, json, os, re
from scipy import stats
from sympy import S
units= dict(brian2.units.__dict__.items() + brian2.units.allunits.__dict__.items() + brian2.__dict__.items())

###################################################################
# PLOT PARAMETER (WILL BE SAVED)                                  #
p={"modfile"          : "cfg/wb01.json",
   "contfile"         : "dat/cont_WB.json",
   "workdir"          : os.getenv('HOME') + /Uni/CNS/3.Semester/labrot2/pythonwrapper#"/pj/channelprcs/sim/",
   "golden"           : 1.618,
   "figname"          : "fig/prcs+noise",
   "formats"          : ['svg','pdf','png'],
   "colordict"        : {
     "prc"   : "#20a386",
     "noise" : "#440154"
     },
   "font.size"        : 10,
   "svg.fonttype"     : "none",
   "axes.labelsize"   : 10,
   "axes.linewidth"   : 0.5,
   "axes.color_cycle" : ["#20a386","#440154"],
   "xtick.labelsize"  : 10,
   "ytick.labelsize"  : 10,
   "xtick.direction"  : "out",
   "ytick.direction"  : "out",
   "xtick.major.size" : 1.5,
   "ytick.major.size" : 1.5,
   "figure.figsize"   : [1.5,1.5]
  }
#############

os.chdir(p['workdir'])
matplotlib.rcParams.update([(j,k) for (j,k) in p.items() if j in matplotlib.rcParams.keys()])

pylab.close('all')
pylab.locator_params(nbins=4)
fig = pylab.figure(1)
fig.clf()

with open(p['contfile']) as f:
  d = json.load(f)

name2pl = 'Kdr1'
diffu = dict(d['diffu'])
diffufoo = sympy.lambdify(map(str,d["coordnames"]),diffu[name2pl],dict(units.items()+[('NKdr',1)]))
impact = []
for j,(par,s) in enumerate(zip(d['par'],d['sol'])):
  ix = d['coordnames'].index('ad'+name2pl)
  diffu_pl = scipy.array([diffufoo(*k[1:]) for k in s])
  s = zip(*s)
  s_pl = scipy.array(s[ix])
  ax = fig.add_axes([0,j*2,1*p['golden'],1])
  ax.plot(s[0],s_pl,color=p['colordict']['prc'])
  # ax.fill_between(s[0],scipy.ones(len(s_pl)),(s_pl*diffu_pl)**2,color=p['colordict']['prc'])
  ax.set_ylim([min(s[ix]),-min(s[ix])])
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)
  ax.xaxis.set_ticks_position("bottom")
  ax.spines["bottom"].set(**{"position":("outward",2)})
  ax.spines["left"].set(**{"position":("outward",2)})
  ax.set_xlabel(r"$\phi$")
  ax.set_ylabel(r"$\Delta\phi$")
  ax2 = ax.twinx()
  ax2.plot(s[0],diffu_pl,color=p['colordict']['noise'])
  ax2.set_ylim([-ax2.get_ylim()[1],ax2.get_ylim()[1]])
  ax2.spines["top"].set_visible(False)
  ax2.spines["right"].set(**{"position":("outward",2)})
  ax2.spines["bottom"].set_visible(False)
  ax2.spines["left"].set_visible(False)
  impact += scipy.trapz((diffu_pl*s_pl)**2,s[0])
  ax2.tick_params(axis="both",which="both",color=p['colordict']['noise'])

for k in p['formats']:
  fig.savefig(p['figname']+'.'+k,bbox_inches="tight")
