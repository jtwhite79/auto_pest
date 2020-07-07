import os
import shutil
import string
import platform

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
import pyemu

abet = string.ascii_uppercase

unit_dict = {"head":"sw-gw flux $\\frac{ft^3}{d}$",
                "tail": "sw-gw flux $\\frac{ft^3}{d}$",
                "trgw" : "gw level $ft$",
                "gage" : "sw flux $\\frac{ft^3}{d}$"}
label_dict = {"head": "headwater",
             "tail": "tailwater",
             "trgw_2_2_9": "gw_1",
              "trgw_2_33_7": "gw_2",
              "trgw_0_9_1" : "gw_3",
             "gage": "sw_1"}

def prep_mf6_model(org_ws):
    
    #if not os.path.exists(os.path.join(org_ws,"freyberg6.nam")):
    if True:
        #first mod the nam file and the dumbass last reach bottom - sigh
        m = flopy.modflow.Modflow.load("freyberg.nam",model_ws=org_ws,check=False, verbose=True)
        print(m.sfr.reach_data.dtype)
        last_reach_bot = m.sfr.reach_data["strtop"][-1] - m.sfr.reach_data["strthick"][-1]
        cell_bot = m.dis.botm[0].array[m.sfr.reach_data["i"][-1],m.sfr.reach_data["j"][-1]]
        print(cell_bot,last_reach_bot)
        if last_reach_bot <= cell_bot:
            m.sfr.reach_data["strtop"][-1] += 0.001
        #m.wel.stress_period_data[7]["flux"] = m.wel.stress_period_data[8]["flux"] * 1.01
        for kper in range(1,m.nper):
            m.wel.stress_period_data[kper]["flux"] = m.wel.stress_period_data[0]["flux"] * (1+ (np.random.random() * 0.001))
            #m.rch.rech[kper] = m.rch.rech[0].array * (1 + (np.random.random() * 0.001))
            #m.wel.stress_period_data[kper]["flux"] *= (1+ (np.random.random() * 0.001))
            m.rch.rech[kper] *= (1 + (np.random.random() * 0.001))

        m.external_path = "."

        m.write_input()
        lines = open(os.path.join(org_ws,"freyberg.nam"),'r').readlines()
        with open(os.path.join(org_ws,"freyberg_mod.nam"),'w') as f:
            for line in lines:
                if ".sfr.out" in line and not "REPLACE" in line:
                    line = line.strip() + " REPLACE\n"
                f.write(line)
        pyemu.os_utils.run("mf5to6 freyberg_mod.nam freyberg6",cwd=org_ws)

    new_ws = org_ws + "_mf6"
    if os.path.exists(new_ws):
        shutil.rmtree(new_ws)
    os.mkdir(new_ws)
    sim = flopy.mf6.MFSimulation.load(sim_ws=org_ws)
    sim.simulation_data.mfpath.set_sim_path(new_ws)
    sim.set_all_data_external()

    sim.name_file.continue_ = True
    m = sim.get_model("freyberg6")

    redis_fac = m.dis.nrow.data / 40 #in case this is a finely discret version

    obs_df = pd.read_csv("obs_loc.csv")
    obs_df.loc[:,"row"] = (obs_df.row * redis_fac).apply(np.int)
    obs_df.loc[:, "col"] = (obs_df.col * redis_fac).apply(np.int)


    obs_df.loc[:,"layer"] = 3
    obs_df.loc[:,"name"] = obs_df.apply(lambda x: "trgw_{0}_{1}_{2}".\
        format(x.layer-1,x.row-1,x.col-1),axis=1)
    obs_df.loc[:,"obstype"] = "HEAD"
    obs_df.loc[:,"col"] = obs_df.col.apply(np.int)
    obs_df2 = obs_df.copy()
    obs_df2.loc[:,"layer"] = 1
    obs_df2.loc[:,"name"] = obs_df2.apply(lambda x: "trgw_{0}_{1}_{2}".\
        format(x.layer-1,x.row-1,x.col-1),axis=1)

    obs_df = pd.concat([obs_df, obs_df2])
    obs_df['idxs'] = list(zip(obs_df.layer - 1, obs_df.row - 1, obs_df.col - 1))
    obs = flopy.mf6.ModflowUtlobs(
        m, pname='head_obs', digits=10, print_input=True,
        continuous={
            'heads.csv':
                obs_df.loc[:, ["name", "obstype", "idxs"]].to_records(
                    index=False)})


    m.sfr.boundnames = True
    m.sfr.obs.initialize(continuous={
        'sfr.csv': [('headwater', 'SFR', 'headwater'),
                    ('tailwater', 'SFR', 'tailwater'),
                    ('gage_1', 'inflow', 40*redis_fac)]}, filename='sfr.obs')
    pkd = m.sfr.packagedata.array.copy()
    pkd['boundname'] = 'headwater'
    pkd['boundname'][pkd.rno > int(m.dis.nrow.data / 2)] = 'tailwater'
    m.sfr.packagedata = pkd

    sim.set_all_data_external()
    sim.write_simulation()
    #shutil.copy2(os.path.join(org_ws,"mf6.exe"),os.path.join(new_ws,"mf6.exe"))
    pyemu.os_utils.run("mf6",cwd=new_ws)
    #ext_dict = {""}



    if "monthly" in org_ws:
        # sample the recharge arrays from the daily model to the monthly model
        daily_ws = "temp_daily_test"
        assert os.path.exists(daily_ws)
        rch_files = [f for f in os.listdir(daily_ws) if "rch_recharge_" in f and f.endswith(".txt")]
        rch_df = pd.DataFrame({"filename":rch_files})
        rch_df.loc[:,"sp"] = rch_df.filename.apply(lambda x: int(x.split('.')[1].split('_')[-1]))
        print(rch_df.sp)
        td = pd.to_timedelta(rch_df.sp.values,unit='d')
        rch_df.loc[:,"datetime"] = pd.to_datetime("12-31-2017") + td
        rch_df.sort_values(by="sp",inplace=True)
        stop = pd.to_datetime("1-1-2020")
        mn_range = pd.date_range(pd.to_datetime("12-31-2017"),stop,freq='m')
        print(mn_range)
        ss_rch_filename = rch_df.filename[0]
        shutil.copy2(os.path.join(daily_ws,ss_rch_filename),os.path.join(new_ws,ss_rch_filename))
        for kper,(start,end) in enumerate(zip(mn_range[:-1],mn_range[1:])):
            mn_rch_df = rch_df.loc[rch_df.datetime.apply(lambda x: x > start and x <= end),:]
            print(start,end,mn_rch_df)
            tot = 0.0
            for filename in mn_rch_df.filename:
                arr = np.loadtxt(os.path.join(daily_ws,filename))
                tot += arr.mean()
            mn_rch_mean = tot / mn_rch_df.shape[0]
            mn_arr = np.zeros((m.dis.nrow.data,m.dis.ncol.data)) + mn_rch_mean
            np.savetxt(os.path.join(new_ws,"freyberg6.rch_recharge_{0}.txt".format(kper+1)),mn_arr,fmt="%15.6E")
        pyemu.os_utils.run("mf6", cwd=new_ws)

def setup_interface(org_ws,num_reals=100):
    np.random.seed(123456)
    # load the mf6 model with flopy to get the spatial reference
    sim = flopy.mf6.MFSimulation.load(sim_ws=org_ws)
    m = sim.get_model("freyberg6")
    
    # work out the spatial rediscretization factor
    redis_fac = m.dis.nrow.data / 40

    # where the pest interface will be constructed
    template_ws = org_ws.split('_')[1] + "_template"
    
    
    # instantiate PstFrom object
    pf = pyemu.prototypes.PstFrom(original_d=org_ws, new_d=template_ws,
                 remove_existing=True,
                 longnames=True, spatial_reference=m.model_grid,
                 zero_based=False,start_datetime="1-1-2018")
    
    # add observations from the sfr observation output file
    df = pd.read_csv(os.path.join(org_ws, "sfr.csv"), index_col=0)
    pf.add_observations("sfr.csv", insfile="sfr.csv.ins", index_cols="time", 
                        use_cols=list(df.columns.values),
                        prefix="sfr")

    # add observations for the heads observation output file
    df = pd.read_csv(os.path.join(org_ws, "heads.csv"), index_col=0)
    pf.add_observations("heads.csv", insfile="heads.csv.ins", 
                        index_cols="time", use_cols=list(df.columns.values),
                        prefix="hds")

    # the geostruct object for grid-scale parameters
    grid_v = pyemu.geostats.ExpVario(contribution=1.0,a=500)
    grid_gs = pyemu.geostats.GeoStruct(variograms=grid_v)

    # the geostruct object for pilot-point-scale parameters
    pp_v = pyemu.geostats.ExpVario(contribution=1.0, a=2000)
    pp_gs = pyemu.geostats.GeoStruct(variograms=pp_v)

    # the geostruct for recharge grid-scale parameters
    rch_v = pyemu.geostats.ExpVario(contribution=1.0, a=1000)
    rch_gs = pyemu.geostats.GeoStruct(variograms=rch_v)

    # the geostruct for temporal correlation
    temporal_gs = pyemu.geostats.GeoStruct(variograms=pyemu.geostats.ExpVario(contribution=1.0,a=60))
    
    # import flopy as part of the forward run process
    pf.extra_py_imports.append('flopy')

    # use the idomain array for masking parameter locations
    ib = m.dis.idomain[0].array
    
    # define a dict that contains file name tags and lower/upper bound information
    tags = {"npf_k_":[0.1,10.],"npf_k33_":[.1,10],"sto_ss":[.1,10],"sto_sy":[.9,1.1],
            "rch_recharge":[.5,1.5]}
    dts = pd.to_datetime("1-1-2018") + \
          pd.to_timedelta(np.cumsum(sim.tdis.perioddata.array["perlen"]),unit="d")

    # loop over each tag, bound info pair
    for tag,bnd in tags.items():
        lb,ub = bnd[0],bnd[1]
        # find all array based files that have the tag in the name
        arr_files = [f for f in os.listdir(template_ws) if tag in f and f.endswith(".txt")]

        if len(arr_files) == 0:
            print("warning: no array files found for ",tag)
            continue
        
        # make sure each array file in nrow X ncol dimensions (not wrapped)
        for arr_file in arr_files:
            arr = np.loadtxt(os.path.join(template_ws,arr_file)).reshape(ib.shape)
            np.savetxt(os.path.join(template_ws,arr_file),arr,fmt="%15.6E")
        
        # if this is the recharge tag
        if "rch" in tag:
            # add one set of grid-scale parameters for all files
            pf.add_parameters(filenames=arr_files, par_type="grid", par_name_base="rch_gr",
                              pargp="rch_gr", zone_array=ib, upper_bound=ub, lower_bound=lb,
                              geostruct=rch_gs)

            # add one constant parameter for each array, and assign it a datetime
            # so we can work out the temporal correlation
            for arr_file in arr_files:
                kper = int(arr_file.split('.')[1].split('_')[-1]) - 1
                pf.add_parameters(filenames=arr_file,par_type="constant",par_name_base=arr_file.split('.')[1]+"_cn",
                                  pargp="rch_const",zone_array=ib,upper_bound=ub,lower_bound=lb,geostruct=temporal_gs,
                                  datetime=dts[kper])
        # otherwise...
        else:
            # for each array add both grid-scale and pilot-point scale parameters
            for arr_file in arr_files:
                pf.add_parameters(filenames=arr_file,par_type="grid",par_name_base=arr_file.split('.')[1]+"_gr",
                                  pargp=arr_file.split('.')[1]+"_gr",zone_array=ib,upper_bound=ub,lower_bound=lb,
                                  geostruct=grid_gs)
                pf.add_parameters(filenames=arr_file, par_type="pilotpoints", par_name_base=arr_file.split('.')[1]+"_pp",
                                  pargp=arr_file.split('.')[1]+"_pp", zone_array=ib,upper_bound=ub,lower_bound=lb,
                                  pp_space=int(5 * redis_fac),geostruct=pp_gs)

    
    # get all the list-type files associated with the wel package
    list_files = [f for f in os.listdir(org_ws) if "freyberg6.wel_stress_period_data_" in f and f.endswith(".txt")]
    # for each wel-package list-type file 
    for list_file in list_files:
        kper = int(list_file.split(".")[1].split('_')[-1]) - 1
        # add spatially constant, but temporally correlated parameter
        pf.add_parameters(filenames=list_file,par_type="constant",par_name_base="twel_mlt_{0}".format(kper),
                          pargp="twel_mlt".format(kper),index_cols=[0,1,2],use_cols=[3],
                          upper_bound=1.5,lower_bound=0.5, datetime=dts[kper], geostruct=temporal_gs)

        # add temporally indep, but spatially correlated grid-scale parameters, one per well
        pf.add_parameters(filenames=list_file, par_type="grid", par_name_base="wel_grid_{0}".format(kper),
                          pargp="wel_{0}".format(kper), index_cols=[0, 1, 2], use_cols=[3],
                          upper_bound=1.5, lower_bound=0.5)

    # add grid-scale parameters for SFR reach conductance.  Use layer, row, col and reach 
    # number in the parameter names
    pf.add_parameters(filenames="freyberg6.sfr_packagedata.txt", par_name_base="sfr_rhk",
                      pargp="sfr_rhk", index_cols=[0,1,2,3], use_cols=[9], upper_bound=10.,
                      lower_bound=0.1,
                      par_type="grid")

    # add model run command
    pf.mod_sys_cmds.append("mf6")
    
    # build pest control file
    pst = pf.build_pst('freyberg.pst')

    # draw from the prior and save the ensemble in binary format
    pe = pf.draw(num_reals, use_specsim=True)
    pe.to_binary(os.path.join(template_ws, "prior.jcb"))

    # set some algorithmic controls
    pst.control_data.noptmax = 0
    pst.pestpp_options["additional_ins_delimiters"] = ","

    # write the control file
    pst.write(os.path.join(pf.new_d, "freyberg.pst"))

    # run with noptmax = 0
    pyemu.os_utils.run("{0} freyberg.pst".format(
        os.path.join("pestpp-ies")), cwd=pf.new_d)

    # make sure it ran
    res_file = os.path.join(pf.new_d, "freyberg.base.rei")
    assert os.path.exists(res_file), res_file
    pst.set_res(res_file)
    print(pst.phi)

    # if successful, set noptmax = -1 for prior-based Monte Carlo
    pst.control_data.noptmax = -1
    
    # define what file has the prior parameter ensemble
    pst.pestpp_options["ies_par_en"] = "prior.jcb"

    # write the updated pest control file
    pst.write(os.path.join(pf.new_d, "freyberg.pst"))

    #assert pst.phi < 1.0e-5, pst.phi


def run_prior_mc(t_d):
    pyemu.os_utils.start_workers(t_d,"pestpp-ies","freyberg.pst",num_workers=5,worker_root=".",
                                 master_dir=t_d.replace("template","master"))

def make_kickass_figs():
    m_d_c = "monthly_master"
    m_d_f = "daily_master"

    sim = flopy.mf6.MFSimulation.load(sim_ws=m_d_f)
    m = sim.get_model("freyberg6")
    redis_fac = m.dis.nrow.data / 40

    pst_c = pyemu.Pst(os.path.join(m_d_c,"freyberg.pst"))
    pst_f = pyemu.Pst(os.path.join(m_d_f, "freyberg.pst"))

    oe_c = pd.read_csv(os.path.join(m_d_c,"freyberg.0.obs.csv"),index_col=0)
    oe_f = pd.read_csv(os.path.join(m_d_f,"freyberg.0.obs.csv"),index_col=0)

    obs_c = pst_c.observation_data
    hds_c = obs_c.loc[obs_c.obsnme.apply(lambda x: x.startswith("hds")),:].copy()
    hds_c.loc[:, "k"] = hds_c.obsnme.apply(lambda x: int(x.split('_')[2]))
    hds_c.loc[:, "time"] = hds_c.obsnme.apply(lambda x: float(x.split('_')[-1].split(':')[1]))

    #hds_c = hds_c.loc[hds_c.k == 2]

    obs_f = pst_f.observation_data
    hds_f = obs_f.loc[obs_f.obsnme.apply(lambda x: x.startswith("hds")), :].copy()
    hds_f.loc[:, "k"] = hds_f.obsnme.apply(lambda x: int(x.split('_')[2]))
    #hds_f = hds_f.loc[hds_f.k == 2]
    hds_f.loc[:,"i"] = hds_f.obsnme.apply(lambda x: int(x.split('_')[3]))
    hds_f.loc[:, "j"] = hds_f.obsnme.apply(lambda x: int(x.split('_')[4]))
    hds_f.loc[:, "time"] = hds_f.obsnme.apply(lambda x: float(x.split('_')[-1].split(':')[1]))

    hds_f.loc[:,"org_i"] = (hds_f.i / redis_fac).apply(np.int)
    hds_f.loc[:, "org_j"] = (hds_f.j / redis_fac).apply(np.int)
    hds_f.loc[:,"org_obgnme"] = hds_f.apply(lambda x: "hds_usecol:trgw_{0}_{1}_{2}".format(x.k,x.org_i,x.org_j),axis=1)

    gage_c = obs_c.loc[obs_c.obsnme.apply(lambda x: "gage" in x),:].copy()
    gage_c.loc[:,"time"] = gage_c.obsnme.apply(lambda x: float(x.split(':')[-1]))
    gage_f = obs_f.loc[obs_f.obsnme.apply(lambda x: "gage" in x), :].copy()
    gage_f.loc[:, "time"] = gage_f.obsnme.apply(lambda x: float(x.split(':')[-1]))
    gage_c.sort_values(by="time",inplace=True)
    gage_f.sort_values(by="time", inplace=True)

    grp_c = set(hds_c.obgnme.unique())

    grp_f = set(hds_f.org_obgnme.unique())

    assert len(grp_f.symmetric_difference(grp_c)) == 0

    grp_c = hds_c.obgnme.unique()
    grp_c.sort()
    print(grp_c)
    print(label_dict)
    grp_c = [g for g in grp_c if g.replace("hds_usecol:","") in label_dict]
    print(grp_c)

    fig,axes = plt.subplots(len(grp_c)+1,1,figsize=(8,8))
    #tags = ["trgw_2_2_9","trgw_2_33_7","trgw_0_9_1"]
    #labels = ["gw_1","gw_2","gw_3"]
    tags = ["trgw_0_9_1"]
    labels = ["gw_3"]
    np.random.seed(1111)
    for i,(tag,label) in enumerate(zip(tags,labels)):
    #for i,(grp,ax) in enumerate(zip(grp_c,axes)):
        grp = None
        for g in grp_c:
            if tag in g:
                grp = g
                break
        c = hds_c.loc[hds_c.obgnme==grp,:].copy()
        c.sort_values(by="time",inplace=True)
        f = hds_f.loc[hds_f.org_obgnme == grp,:].copy()
        f.sort_values(by="time",inplace=True)
        oe_c_g = oe_c.loc[:,c.obsnme].copy()
        oe_f_g = oe_f.loc[:, f.obsnme].copy()
        oe_c_g.values[oe_c_g.values < 30] = np.nan
        oe_f_g.values[oe_f_g.values < 30] = np.nan

        tr = oe_f_g.iloc[0,:].values.copy()
        tr += np.random.normal(0,0.15,tr.shape[0])

        # lab = None
        # for l,t in label_dict.items():
        #     if l in grp:
        #         lab = t
        ax = axes[i]
        [ax.plot(c.time.values,oe_c_g.loc[i,:].values,color='b',alpha=0.5,lw=0.2) for i in oe_c_g.index]
        ax.plot(f.time.values,tr,color="r",lw=1.5)

        up = oe_c_g.mean() + (3 * oe_c_g.std())
        ax.plot(c.time.values, up.values,color="b",lw=1.0,ls="--")


        #[print(f.time.values, oe_f_g.loc[i, :].values) for i in oe_f_g.index]
        ax.set_ylabel("simulated groundwater level ($L$)")
        ax.set_xlabel("simulation time ($T$)")
        ax.set_title("{0}) {1} lower resolution".format(abet[i], label), loc="left")
        ax = axes[i+1]
        [ax.plot(f.time.values, oe_f_g.loc[i, :].values, color='g', alpha=0.5, lw=0.2) for i in oe_f_g.index]
        ax.plot(f.time.values, tr, color="r", lw=1.5)
        up = oe_f_g.mean() + (3 * oe_f_g.std())
        ax.plot(f.time.values, up.values, color="g", lw=1.0, ls="--")
        ax.set_ylabel("simulated groundwater level ($L$)")
        ax.set_xlabel("simulation time ($T$)")
        ax.set_title("{0}) {1} higher resolution".format(abet[i],label),loc="left")
        #ax.set_ylim(34,40)

    tr = oe_f.loc[:,gage_f.obsnme].iloc[0,:].values.copy()
    tr += (tr * 0.05 * np.random.normal(0,1,tr.shape[0]))
    [axes[-2].plot(gage_c.time.values,oe_c.loc[i,gage_c.obsnme].values,color='b',alpha=0.5,lw=0.2) for i in oe_c.index]
    up = oe_c.loc[:,gage_c.obsnme].mean() + (3 * oe_c.loc[:,gage_c.obsnme].std())
    axes[-2].plot(gage_c.time.values, up.values, color="b", lw=1.0, ls="--")



    axes[-2].plot(gage_f.time.values,tr,color='r',lw=1.5)
    axes[-2].set_ylabel("simulated surface water flow ($\\frac{L^3}{T}$)")
    axes[-2].set_xlabel("simulation time ($T$)")
    axes[-2].set_title("C) sw_1 lower resolution", loc="left")
    axes[-2].set_ylim(1000,5500)

    [axes[-1].plot(gage_f.time.values, oe_f.loc[i, gage_f.obsnme].values, color='g', alpha=0.5, lw=0.2) for i in
     oe_f.index]
    up = oe_f.loc[:, gage_f.obsnme].mean() + (3 * oe_f.loc[:, gage_f.obsnme].std())
    axes[-1].plot(gage_f.time.values, up.values, color="g", lw=1.0, ls="--")

    axes[-1].plot(gage_f.time.values, tr, color='r', lw=1.5)
    axes[-1].set_ylabel("simulated surfacewaer flow ($\\frac{L^3}{T}$)")
    axes[-1].set_xlabel("simulation time ($T$)")
    axes[-1].set_title("D) sw_1 higher resolution", loc="left")
    axes[-1].set_ylim(1000, 5500)
    plt.tight_layout()
    plt.savefig("obs_prior.pdf")
    plt.close(fig)


    fig,axes = plt.subplots(1,2,figsize=(8,4))
    forecasts = ["sfr_usecol:headwater_time:610.0","sfr_usecol:tailwater_time:610.0"]
    labels = ["A) headwater flux","B) tailwater_flux"]
    xmin = min(oe_c.loc[:,forecasts].min().min(),oe_f.loc[:,forecasts].min().min())
    xmax = max(oe_c.loc[:, forecasts].max().max(),oe_f.loc[:, forecasts].max().max())

    for forecast,ax,label in zip(forecasts,axes,labels):
        oe_c.loc[:,forecast].hist(ax=ax,bins=15,edgecolor="none",facecolor="b", alpha=0.5,label="lower resolution",density=True)
        oe_f.loc[:, forecast].hist(ax=ax, bins=15, edgecolor="none", facecolor="g", alpha=0.5, label="higher resolution",density=True)
        ax.set_yticks([])
        ax.set_xlabel("surface-water/groundwater flux ($\\frac{L^3}{T}$)")
        ax.set_title(label,loc="left")
        ax.set_xlim(xmin,xmax)
        ax.legend(loc="upper right")
        ax.grid(False)
        ax.set_ylabel("increasing probability density")
    plt.tight_layout()
    plt.savefig("forecast_prior.pdf")
    plt.close(fig)

def write_par_sum(pst_file):
    pst = pyemu.Pst(pst_file)
    par = pst.parameter_data
    par.loc[par.pargp.str.startswith("wel"),"pargp"] = "wel"
    def group_namer(grp):
        name = ''
        if "layer" in grp:
            name += "Layer " + grp.split('layer')[1][0]
        if "gr" in grp:
            name += " grid-scale"
        elif "pp" in grp:
            name += " pilot points"
        elif "const" in grp:
            name += " constant"
        if "_k_" in grp:
            name += " HK"
        elif "k33" in grp:
            name += " VK"
        elif "ss" in grp:
            name += " SS"
        elif "sy" in grp:
            name += " SY"
        elif "rch" in grp:
            name += " recharge"
        if "sfr" in grp:
            name = " SFR stream-bed conductance"
        if "twel" in grp:
            name = "temporal wel flux constants"
        #elif "wel" in grp:
        #    name = "grid-scale wel flux for stress period " + grp.split('_')[1]
        elif "wel" in grp:
            name = "grid-scale wel flux"
        return name

    ugrps = pst.parameter_data.pargp.unique()
    name_dict = {ug:group_namer(ug) for ug in ugrps}
    pst.write_par_summary_table(os.path.split(pst_file)[0] + "_par_sum.tex",group_names=name_dict)


if __name__ == "__main__":
    #daily has to go first for recharge sampling in monthly
    prep_mf6_model("temp_daily")
    prep_mf6_model("temp_monthly")
    setup_interface("temp_monthly_test",num_reals=100)
    run_prior_mc("monthly_template")
    # #

    setup_interface("temp_daily_test",num_reals=100)
    run_prior_mc("daily_template")
    #
    make_kickass_figs()
    #write_par_sum(os.path.join("monthly_master","freyberg.pst"))
    #write_par_sum(os.path.join("daily_master", "freyberg.pst"))

