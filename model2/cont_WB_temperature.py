## Continue in $I_\mathrm{dc}$ to illustrate frequency dependence ##
## Simulation must be run first (p["simfile"])                    ##
####################################################################

def cont_WB_temperature_f(hb_continuation, outputname, gK, gNa):

    #Input parameters:
        # hb_continuation: What do we want to do?
        #        0: Cont Hopf in Cm
        #        1: Cont Hopf in T
        #        2: Cont Hopf in gK
        #        3: Cont Hopf in gK, followed by several conts in gNa
        #        4: Cont Hopf in gNa
        # outputname: Name of resultfiles that will be written

    import matplotlib
    import ipdb
    #matplotlib.use("Agg")
    import pylab, brian2, brianutils, os, json, sympy, scipy, sys, \
        datetime, shelve, autoutils, auto, contextlib, inibrian
    from sympy import S


    ########### Get initial condition through brian ###########

    bifpar = {
      "I" : ["0.0* uA/cm2"],
      "Cm" : ["1.0*uF/cm2","1.1*uF/cm2","1.4*uF/cm2","1.41*uF/cm2","1.42*uF/cm2"],
      "dT" : ["0. * kelvin"],
      "gK" : [str(gK) + " * msiemens/cm2"],
      "gNa": [str(gNa) + " * msiemens/cm2"]#,
      #"q_gk": ["1.2 * kelvin/kelvin"]
      }
    unames, pnames, ini, autobifpar = inibrian.find_ini_brian('wangBuzsaki_brian_temperature.json', 'model2', bifpar, 'wbtemp')

    autobifpar['Cm'] = 1
    autobifpar['gK'] = gK
    autobifpar['gNa'] = gNa
    ################# CONT FP & LC ############################


    calc_period = 0 #do we want to calculate orbits?
    refine_hb = 0 #do we want to get the HB position even more exactly? (Not needed and also not really working :))



    directionlist = [1, -1]

    r1_fwd= auto.run(ini, e='wbtemp',
        c='wbtemp', parnames= pnames, unames=unames,
        ICP=['I'], ISP=1,ILP=1, SP=['LP','HB','BP'],
        PAR=autobifpar, ITNW=17, NWTN=13, NMX = 100000,NPR=100000, #NMX=113500, NPR=5000,
        DS= directionlist[0] * 1e-3, DSMAX= directionlist[0] * 1e-2, STOP=['HB1'],
        UZSTOP= {})

    r1_bwd= auto.run(ini, e='wbtemp',
        c='wbtemp', parnames= pnames, unames=unames,
        ICP=['I'], ISP=1,ILP=1, SP=['LP','HB','BP'],
        PAR=autobifpar, ITNW=17, NWTN=13, NMX = 3000, NPR=3000,#NMX=113500, NPR=5000,
        DS= directionlist[1] * 1e-3, DSMAX= directionlist[1] * 1e-2, STOP=['HB1'],
        UZSTOP= {})

    r1 = r1_fwd + r1_bwd
    s1HB = r1_fwd.getLabel('HB')[0]
    #s1LP = r1.getLabel('LP')[0]
    #s1UZ = r1.getLabel('EP')[0]

    ################# Optional: continue orbit #######################

    if calc_period:
        r1_period = auto.run(s1HB, e='wbtemp', c='wbtemp',
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
            r21= auto.run(s1HB, e='wbtemp',
                c='wbtemp', parnames= pnames, unames=unames,
                ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=10000, NPR=10000,
                DS=-1e-5, DSMAX=-1e-4, STOP=['HB1'],
                IID = 3)
            r22= auto.run(s1HB, e='wbtemp',
                c='wbtemp', parnames= pnames, unames=unames,
                ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=10000, NPR=10000,
                DS=1e-5, DSMAX=1e-4, STOP=['HB1'],
                IID = 3)
            r2 = r21 + r22
            s2HB = r2.getLabel('HB')[0]
            r31= auto.run(s2HB, e='wbtemp',
                c='wbtemp', parnames= pnames, unames=unames,
                ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=10000, NPR=10000,
                DS=-1e-5, DSMAX=-1e-4, STOP=['HB1'],
                IID = 3)
            r32= auto.run(s2HB, e='wbtemp',
                c='wbtemp', parnames= pnames, unames=unames,
                ICP=['I','Cm','period'], ISP=1,ILP=1, SP=['LP','HB','BP'], ITNW=17, NWTN=13, NMX=10000, NPR=10000,
                DS=1e-5, DSMAX=1e-4, STOP=['HB1'],
                IID = 3)
            r3 = r31 + r32
            s3HB = r3.getLabel('HB')[0]
            s1HB = s3HB
        except:
            print("Could not refine position: Bifurcation point not found. Continuing with unrefined HB.")

    ################# CONT HB ############################


    if (hb_continuation == 0):
        runHB = [0 for item in directionlist]
        for counter, direction in enumerate(directionlist):
            runHB[counter] = auto.run(s1HB, e='wbtemp', c='wbtemp',
                    parnames=pnames, unames=unames,
                    ICP=['Cm', 'I'], #IPS = 1
                    ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
                    ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
                    ITNW=17, NWTN=13, NMX=10000, NPR=1000, #Not sure if needed: PAR=autobifpar,
                    DS= direction * 1e-2, DSMAX= direction * 1e-1, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
                    UZSTOP= {}
                    )

    elif (hb_continuation == 1):
        runHB = [0 for item in directionlist]
        for counter, direction in enumerate(directionlist):
            runHB[counter] = auto.run(s1HB, e='wbtemp', c='wbtemp',
                    parnames=pnames, unames=unames,
                    ICP=['dT', 'I'], #IPS = 1
                    ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
                    ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
                    ITNW=17, NWTN=13, NMX=10000, NPR=1000, #Not sure if needed: PAR=autobifpar,
                    DS= direction * 1e-2, DSMAX= direction * 1e-1, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
                    UZSTOP= {'dT': direction * 90}
                    )
        runHB_fwdbwd = runHB[0] + runHB[1]

    elif (hb_continuation == 2):
        runHB = [0 for item in directionlist]
        for counter, direction in enumerate(directionlist):
            runHB[counter] = auto.run(s1HB, e='wbtemp', c='wbtemp',
                    parnames=pnames, unames=unames,
                    ICP=['gK', 'I'], #IPS = 1
                    ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
                    ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
                    ITNW=17, NWTN=13, NMX=1000, NPR=100, #Not sure if needed: PAR=autobifpar,
                    DS= direction * 1e-2, DSMAX= direction * 1e-1, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
                    UZSTOP= {}
                    )
        runHB_fwdbwd = runHB[0] + runHB[1]

    elif (hb_continuation == 3):
        runHB_gk = [0 for item in directionlist]
        runHB_gna = [0 for item in directionlist]
        #continue in gk
        for counter, direction in enumerate(directionlist):
            runHB_gk[counter] = auto.run(s1HB, e='wbtemp', c='wbtemp',
                    parnames=pnames, unames=unames,
                    ICP=['gK', 'I', 'gNa'], #IPS = 1
                    ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
                    ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
                    ITNW=17, NWTN=13, NMX=1000, NPR=100, #Not sure if needed: PAR=autobifpar,
                    DS= direction * 1e-2, DSMAX= direction * 1e-1, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
                    UZSTOP= {'gK': direction * 75}
                    )
            runHB_gna[counter] = auto.run(s1HB, e='wbtemp', c='wbtemp',
                    parnames=pnames, unames=unames,
                    ICP=['gNa', 'I', 'gK'], #IPS = 1
                    ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
                    ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
                    ITNW=17, NWTN=13, NMX=1000, NPR=100, #Not sure if needed: PAR=autobifpar,
                    DS= direction * 1e-2, DSMAX= direction * 1e-1, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
                    UZSTOP= {'gNa': direction * 75}
                    )

        #take every gk point and continue it in gN
        runHB_fwdbwd_gk = runHB_gk[0] + runHB_gk[1]

        runHB_fwdbwd_gk = auto.relabel(runHB_fwdbwd_gk)
        start_idxs = runHB_fwdbwd_gk.getLabels()

        for pointcounter, startpoint in enumerate(start_idxs):
            for listcounter, direction in enumerate(directionlist):
                runHB_fwdbwd_gk += auto.run(runHB_fwdbwd_gk.getLabel(startpoint), e='wbtemp', c='wbtemp',
                        parnames=pnames, unames=unames,
                        ICP=['gNa', 'I', 'gK'], #IPS = 1
                        ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
                        ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
                        ITNW=17, NWTN=13, NMX=1000, NPR=100, #Not sure if needed: PAR=autobifpar,
                        DS= direction * 1e-2, DSMAX= direction * 1e-1, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
                        UZSTOP= {'gK': 0.1}
                        )
        #oppsite: take every gna point and cont it in gK
        runHB_fwdbwd_gna = runHB_gna[0] + runHB_gna[1]
        runHB_fwdbwd_gna = auto.relabel(runHB_fwdbwd_gna)
        start_idxs = runHB_fwdbwd_gna.getLabels()
        #run_gkna = [0 for i in range(len(start_idxs)*2)]
        for pointcounter, startpoint in enumerate(start_idxs):
            for listcounter, direction in enumerate(directionlist):
                runHB_fwdbwd_gna += auto.run(runHB_fwdbwd_gna.getLabel(startpoint), e='wbtemp', c='wbtemp',
                        parnames=pnames, unames=unames,
                        ICP=['gK', 'I', 'gNa'], #IPS = 1
                        ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
                        ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
                        ITNW=17, NWTN=13, NMX=1000, NPR=100, #Not sure if needed: PAR=autobifpar,
                        DS= direction * 1e-2, DSMAX= direction * 1e-1, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
                        UZSTOP= {'gK': 0.1}
                        )
        runHB_fwdbwd = runHB_fwdbwd_gk + runHB_fwdbwd_gna
        test0 = runHB_fwdbwd['I']
        test1 = runHB_fwdbwd['gNa']
        runHB_fwdbwd = auto.klb(runHB_fwdbwd)
        runHB_fwdbwd = auto.dlb(runHB_fwdbwd)
        #take every gna point and continue it in gK

    elif (hb_continuation == 4):
        runHB = [0 for item in directionlist]

        #continue in gk
        for counter, direction in enumerate(directionlist):
            runHB[counter] = auto.run(s1HB, e='wbtemp', c='wbtemp',
                    parnames=pnames, unames=unames,
                    ICP=['gNa', 'I', 'gK'], #IPS = 1
                    ISP=2,ILP=1,# SP=['LP','HB','BP'],   # ISP: Bifurcation detection; 0=off, 1=BP(FP), 3=BP(PO,BVP), 2=all | ILP: Fold detection; 1=on, 0=off
                    ISW=2, # ISW: Branch switching; 1=normal, -1=switch branch (BP, HB, PD), 2=switch to two-parameter continuation (LP, BP, HB, TR), 3=switch to three-parameter continuation (BP)
                    ITNW=17, NWTN=13, NMX=1000, NPR=100, #Not sure if needed: PAR=autobifpar,
                    DS= direction * 1e-2, DSMAX= direction * 1e-1, #IID = 5, RL1=0.58, EPSL = 1e+1, EPSU = 1e+1, EPSS = 10, NTST = 500,
                    UZSTOP= {'gNa': direction * 75}
                    )
        runHB_fwdbwd = runHB[0] + runHB[1]

    else:
        print("Please specify HB cont mode!")

    rhbt = runHB_fwdbwd + r1
    auto.save(rhbt, outputname) #for "Run HB temperature"
    ipdb.set_trace()



cont_WB_temperature_f(3, 'cont_grid', 9, 35)
#cont_WB_temperature_f(3, 'cont_grid2', 9, 35)
# hb_continuation: What do we want to do?
#        0: Cont Hopf in Cm
#        1: Cont Hopf in T
#        2: Cont Hopf in gK
#        3: Cont Hopf in gK and gNa followed by several conts in gK and gNa -> grid
#        4: Cont Hopf in gNa
