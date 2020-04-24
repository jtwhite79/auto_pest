import os
import shutil
import string
import numpy as np
import pandas as pd
import platform
import flopy
import pyemu


def prep_mf6_model(org_ws):
    
    if True:#not os.path.exists(os.path.join(org_ws,"freyberg6.nam")):

        #first mod the nam file and the dumbass last reach bottom - sigh
        m = flopy.modflow.Modflow.load("freyberg.nam",model_ws=org_ws,check=False)
        print(m.sfr.reach_data.dtype)
        last_reach_bot = m.sfr.reach_data["strtop"][-1] - m.sfr.reach_data["strthick"][-1]
        cell_bot = m.dis.botm[0].array[m.sfr.reach_data["i"][-1],m.sfr.reach_data["j"][-1]]
        print(cell_bot,last_reach_bot)
        if last_reach_bot <= cell_bot:
            m.sfr.reach_data["strtop"][-1] += 0.001
        m.wel.stress_period_data[7]["flux"] = m.wel.stress_period_data[8]["flux"] * 1.01
        m.external_path = "."

        m.write_input()
        lines = open(os.path.join(org_ws,"freyberg.nam"),'r').readlines()
        with open(os.path.join(org_ws,"freyberg_mod.nam"),'w') as f:
            for line in lines:
                if ".sfr.out" in line and not "REPLACE" in line:
                    line = line.strip() + " REPLACE\n"
                f.write(line)
        pyemu.os_utils.run("mf5to6 freyberg_mod.nam freyberg6",cwd=org_ws)


    new_ws = org_ws + "_test"
    if os.path.exists(new_ws):
        shutil.rmtree(new_ws)
    os.mkdir(new_ws)
    sim = flopy.mf6.MFSimulation.load(sim_ws=org_ws)
    sim.simulation_data.mfpath.set_sim_path("test")
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

    obs_df = pd.concat([obs_df,obs_df2])

    #head_obs = {"head_obs.csv":[("trgw_{0}_{1}".format(r,c),"NPF",(2,r-1,c-1)) for r,c in zip(obs_df.row,obs_df.col)]}
    with open(os.path.join(new_ws,"head.obs"),'w') as f:
        f.write("BEGIN CONTINUOUS FILEOUT heads.csv\n")
        obs_df.loc[:,["name","obstype","layer","row","col"]].to_csv(f,sep=' ',line_terminator='\n',
            index=False,header=False,mode="a")
        f.write("END CONTINUOUS\n")

    props_3d = [("npf","k"),("npf","k33"),("sto","ss"),("sto","sy")]
    props_trans = [("rch","recharge")]
    #print(dir(m.dis.nlay))
    #print(m.dis.nlay.data)
    print(type(m.rch.recharge))

    for pack,attr in props_3d:
        #print(m.get_package(pack).__getattribute__(attr))
        for k in range(m.dis.nlay.data):
            filename = "{0}_{1}_{2}.dat".format(pack,attr,k)
            m.get_package(pack).__getattribute__(attr).store_as_external_file(filename,layer=k)

    sim.write_simulation()
    lines = open(os.path.join(new_ws,"freyberg6.nam"),'r').readlines()
    new_lines = []
    for line in lines:
        if line.strip() == "END Packages":
            new_lines.append("   obs6 head.obs\n")
        new_lines.append(line)
    with open(os.path.join(new_ws,"freyberg6.nam"),'w') as f:
        [f.write(line) for line in new_lines]
        f.flush()

    # mod sfr
    # lines = open(os.path.join(new_ws,"freyberg6.sfr"),'r').readlines()
    # with open(os.path.join(new_ws,"freyberg6.sfr"),'w') as f:
    #     iline = 0
    #     while True:
    #         if iline >= len(lines):
    #             break
    #         line = lines[iline]
    #         f.write(line)
    #         if "begin options" in line.lower():
    #             f.write("  BOUNDNAMES\n")
    #             f.write("OBS6 FILEIN sfr.obs\n")
    #         if "begin packagedata" in line.lower():
    #             iline += 1
    #             #line = lines[iline]
    #             irch = 0
    #             while "end" not in line.lower():
    #                 print(iline)
    #                 line = lines[iline]
    #                 if "end" in line.lower():
    #                     break
    #                 bn = "headwater"
    #                 if irch > int(m.dis.nrow.data / 2):
    #                     bn = "tailwater"
    #                 line = line.strip() + " " + bn + "\n"
    #                 f.write(line)
    #                 irch += 1
    #                 iline += 1
    #             f.write(line)
    #         iline += 1

    lines = open(os.path.join(new_ws,'freyberg6.sfr_packagedata.txt'),'r').readlines()
    with open(os.path.join(new_ws,'freyberg6.sfr_packagedata.txt'),'w') as f:
        for i,line in enumerate(lines):
            bn = "headwater"
            if i > int(m.dis.nrow.data / 2):
                bn = "tailwater"
            f.write(line.strip() + " " + bn + "\n")


    with open(os.path.join(new_ws,"sfr.obs"),'w') as f:
        f.write("BEGIN CONTINUOUS FILEOUT sfr.csv\nheadwater sfr headwater\n")
        f.write("tailwater sfr tailwater\ngage_1 inflow 40\nEND CONTINUOUS")

    #shutil.copy2(os.path.join(org_ws,"mf6.exe"),os.path.join(new_ws,"mf6.exe"))
    pyemu.os_utils.run("mf6",cwd=new_ws)
    #ext_dict = {""}


if __name__ == "__main__":
    prep_mf6_model("temp_daily")
