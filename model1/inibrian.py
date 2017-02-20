def find_ini_brian():

    import pylab, brian2, brianutils, os, json, sympy, scipy, sys, \
        datetime, shelve, autoutils, contextlib
    from sympy import S

    units= dict(brian2.units.__dict__.items()
                + brian2.units.allunits.__dict__.items()
                + brian2.__dict__.items())
    unitlist=["mV","ms","cm2","uF","psiemens","um2","msiemens","cm"]
    baseunits=[(k,float(eval(k,units).base)) for k in unitlist]

    p = {"simfile"  : "sim_DNap_2015-12-14_16:55:14.603625.shv",
         "modfile"  : "cfg/wangBuzsaki_brian.json",      #This is our model
         "contfile" : "dat/cont_WB.json",   #This is what brian produces and hands to auto
         "workdir"  : os.getenv("HOME") + "/Uni/CNS/3.Semester/schreiberlab/excitationblock/model1",
         "dt"       : "0.05*ms",
         "bifpar"   : {
           "I" : ["0.0* uA/cm2"],
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

    baseunits2 = [('mV', 1), ('ms', 1), ('cm2', 1), ('uF', 1), ('psiemens', 1), ('um2', 1), ('msiemens', 1), ('cm', 1)]
    varrhs = [(i,sympy.S(j).subs(baseunits2))
                    for i,j in ode.eq_expressions]

    # varrhs_2 = [(i,sympy.S(j).subs(baseunits))
    #                 for i,j in ode.eq_expressions]

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

    unames,pnames= autoutils.writeFP('tm_new',
        bifpar=autobifpar, rhs=rhs, var=var,
        bc=['{0}_left-{0}_right'.format(v) for v in var] + spikecriterion,
        ic=[])

    inivals = ([float(getattr(states,j)[0][-1]) for j in var])

    #convert first value (V) from V to mV.
    inivals[0] *= 1000
    return unames, pnames, inivals, autobifpar
