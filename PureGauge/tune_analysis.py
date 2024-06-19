#!/usr/bin/env python3
# coding: utf-8

import re
import numpy as np
import pandas as pd
import os
import getopt
import sys
import pickle

def stepcounter(idx:int,idxmax:list):
    lidx=idx+1
    lst=[]

    for i in range(len(idxmax)):
        lst=[np.mod(lidx,idxmax[-i-1])] + lst
        lidx=lidx//idxmax[-i-1]
    return (idx+1,lst)
    


def read_outfile(linfile_name:str,nopickle:bool=False):
    run_info={}
    header=True

    if os.path.isfile(linfile_name+".pkl") and not nopickle:
        print(f"Found the pickle preproces {linfile_name+'.pkl'}")
        with open(linfile_name+".pkl", 'rb') as fp:
            run_info,data_df= pickle.load(fp)
        print(f"Memory footprint of the df in Mb={data_df.memory_usage(deep=True).sum()/10**6}")
        return run_info,data_df
    
    if not os.path.isfile(linfile_name):
        print(f"Unable to find {linfile_name}")
        return
        
    with open(linfile_name) as lfile:
        for lline in lfile:
            if("[INIT ML][0]number of MultiLevels=" in lline and header):
                run_info["n mlev"] = int(re.split("=",lline)[1])
                run_info["nit per ml"]=[None]* run_info["n mlev"]
                gbidx=toridx=-1
                run_info["frozen t"]={}
                for lid in range(run_info["n mlev"]+1):
                    lsize =  run_info["glbt"]//2**(lid+1)
                    run_info["frozen t"][lid] = [(i+1) * lsize -1 for i in range(2**(lid+1)-1) ]
            elif("[GEOMETRY_INIT][0]Global size is" in lline and header):
                stringsplit = re.split("x| ",lline)
                run_info["glbt"]=int(stringsplit[3])
                print(f'Global lattice temporal extent T={run_info["glbt"]}')
            elif("[INIT ML][0]ML tuning level=" in lline and header):
                tunelev = int(re.split("=",lline)[1]) 
                run_info["nit per ml"]=run_info["nit per ml"][:(tunelev+1)]
                print(f"Output file contain tuning on lev={tunelev}")
                print(f'Each conf contains {np.prod(run_info["nit per ml"])} different measures split along {tunelev+1} remaining update levels')
                store=[]
            elif("[INIT ML][0]lev " in lline and header):
                stringsplit = re.split("=| ",lline)
                llevid=int(stringsplit[2])
                lup=int(stringsplit[4])
                run_info["nit per ml"][llevid]=lup

            elif("[MAIN][0]Configuration " in lline):
                cnfgid=int(lline.split()[1])
                assert(np.mod(gbidx,np.prod(run_info["nit per ml"]))== np.prod(run_info["nit per ml"])-1)
                assert(np.mod(toridx,np.prod(run_info["nit per ml"]))== np.prod(run_info["nit per ml"])-1)
                header=False
            elif("[MEASURE_TUNE][0]Glueball operator tune measure" in lline):
                (gbidx,idxlst) = stepcounter(gbidx,run_info["nit per ml"])
            elif("[MEASURE_TUNE][0]Torellon operator tune measure" in lline):
                (toridx,idxlst) = stepcounter(toridx,run_info["nit per ml"])
            elif("[Measure ML][0]1pt function P=" in lline or "[Measure ML][0]1ptTor function P=" in lline):
                stringsplit = lline.split()
                if(stringsplit[1]=='ML][0]1pt'):
                    typestring="glueball"
                else:   
                    typestring="torellon"
                plst = list(map(lambda x:int(x), re.split("P=\\(|,|\\)", stringsplit[3] )[1:-1]))
                irrepstring = re.split("=", stringsplit[4] )[1]
                evstring = re.split("=", stringsplit[6] )[1]
                chargestring = re.split("=", stringsplit[7] )[1]
                data_dic={"cfgid":cnfgid,"type":typestring,"momenta": plst ,"irrep":irrepstring,"charge":chargestring,"ev":evstring}
                for i in range(tunelev+1):
                    data_dic["ML_index"+str(i)]=idxlst[i]
            elif("[Measure ML][0] t=" in lline):
                stringsplit = re.split("=|\\(|\\)| ",lline)[3:]
                data_dic["t"]=int(stringsplit[0])
                data_dic["values"]=[float(i) for i in stringsplit[1:-1] if i ]
                store.append(data_dic.copy())
    if(toridx!=-1):
        print(f'Found {(toridx+1)//np.prod(run_info["nit per ml"])} torellons configuration measurements')
    if(gbidx!=-1):
        print(f'Found {(gbidx+1)//np.prod(run_info["nit per ml"])} glueballs configuration measurements')
    data_df = pd.DataFrame(store)
    with open(linfile_name+'.pkl', 'wb') as fp:
        pickle.dump((run_info,data_df), fp)
    print(f"Memory footprint of the df in Mb={data_df.memory_usage(deep=True).sum()/10**6}")
    return run_info,data_df
 
def data_select(dataf,ltype,irrep,charge,ev):
    data_df2=dataf[ ( dataf["type"]==ltype ) & ( dataf["irrep"]==irrep ) & ( dataf["charge"]==charge ) & ( dataf["ev"]==ev ) ]
    assert(data_df2.size>0)
    print(f'Available t indices {data_df2["t"].unique()}')
    nop=data_df2['values'].apply(lambda x: len(x)).unique()
    assert(len(nop)==1)
    print(f'Number of available operators {nop[0]}')
    return data_df2
    
def main():
    try:
        opt, args = getopt.getopt(sys.argv[1:], 'f:m:', ('File to pickle',"list of pkl files to merge"))
    except:
        print('Something wrong in your options')
        sys.exit()
    
    if len(opt)!= 1:
        print('Select one option (either f or m)')
        sys.exit()
        
    for o, a in opt:
        if o == '-f':
            file = a
            action="pickle"
        elif o == "-m":
            listfile = a
            action="merge"

        
    if(action =="pickle"):
        print(f"Reading and pickling {file}")
        read_outfile(a)
    if(action =="merge"):
        print(f"Merging the files in {listfile}")
        counter = 1
        optype="glueball"
        irrep="A1plusOhP"
        charge="+"
        ev="1of1"
        outfilename = f"out_tune_measure_list_{optype}_{irrep}_{charge}_{ev}"
        with open(listfile, "r") as list_file:
            file=list_file.readline().strip()
            rn_inf,data_df1=read_outfile(file)
            data_df2=data_select(data_df1 ,  optype ,irrep , charge , ev.replace("of","/")  )
            for file in list_file:
                counter = counter +1
                del data_df1
                _,data_df1=read_outfile(file.strip())
                data_df2= pd.concat([data_df2,data_select(data_df1 , optype ,irrep , charge , ev.replace("of","/") )])
                if counter % 10 == 0 :
                    with open(f'{outfilename}_{counter}.pkl', 'wb') as fp:
                        pickle.dump((rn_inf,data_df2), fp)
        with open(f'{outfilename}.pkl', 'wb') as fp:
            pickle.dump((rn_inf,data_df2), fp)
        print(f"Memory footprint of the final df in Mb={data_df2.memory_usage(deep=True).sum()/10**6}")
        
main()