#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-

import argparse
import os.path
import sys
import re
import numpy as np
from operator import itemgetter


from pylab import *
from subprocess import *

bstat="./bs_stat "
trantab={ord('('): None,ord(')'): None}
trantab2={ord('/'):"x"}
trantab3={ord('\n'):None}
trantab4={ord('='):" "}

def correlator_def(l_inputfile,g_correlator_list):
    tmpfile=open(l_inputfile,'r')
    for line in tmpfile:
        if("[INIT ML][0] Cor Id=" in line):
            infostring=line.split("[INIT ML][0] Cor Id=",1)[1]
            mykey=re.split('\ |\=|\(|\)', infostring)
            if mykey[0] not in g_correlator_list.keys():
                g_correlator_list[mykey[0]]=[]
            g_correlator_list[mykey[0]].append([mykey[6],mykey[7]])
        if("[MAIN][0]ML Measure #" in line):
            break 

def correlator_meas(nmeas,corr_data,sub_data,l_correlator_list,file):
    for id in l_correlator_list.keys():
        sum=0.
        dt=abs(int(l_correlator_list[id][0][1])-int(l_correlator_list[id][0][0]))
        npoints=len(l_correlator_list[id])
        if not sub_data:
            for pair in l_correlator_list[id]:
                sum += 1.0/(2*npoints)*(np.outer(np.conjugate(corr_data[pair[0]]), corr_data[pair[1]])+
                    np.outer(np.conjugate(corr_data[pair[1]]), corr_data[pair[0]])).real
        else:
             for pair in l_correlator_list[id]:
 #               print("qui",corr_data[pair[1]],sub_data[pair[1]],corr_data[pair[0]],sub_data[pair[0]])
 #               sum += (1.0/(2*npoints)*(np.outer(np.conjugate(corr_data[pair[0]]), corr_data[pair[1]])+
 #                   np.outer(np.conjugate(corr_data[pair[1]]), corr_data[pair[0]])).real+
 #                   -0.0/(2*npoints)*(np.outer(sub_data[pair[0]],sub_data[pair[1]])
 #                   +np.outer(np.conjugate(sub_data[pair[1]]), sub_data[pair[0]])).real)
                sum += 1.0/(2.0*npoints)*(np.outer(np.conjugate(corr_data[pair[0]] - sub_data[pair[0]]), corr_data[pair[1]]- sub_data[pair[1]])+
                        np.outer(np.conjugate(corr_data[pair[1]]- sub_data[pair[1]]), corr_data[pair[0]]- sub_data[pair[0]])).real
        file.write(str(nmeas)+" "+str(dt)+"\n")
        np.savetxt(file, sum)


def onept_outfilename(instring):
    head=instring.split("[Measure ML][0]1pt function ",1)[1]
    head=separator.join(list(head.split(' ')[i] for i in [-6,-5,-2,-4,-3]))
    return prepend+("1pt_"+head+".dat").translate(trantab2)


def find_outfile_list(instring,linfile_name):
    lnmeas=0
    file_irreps=[]
    with open(linfile_name) as lfile:
        for lline in lfile:
            if("[MAIN][0]ML Measure #" in lline and "..." in lline):
                lnmeas+=1
            elif("[Measure ML][0]1pt function" in lline):
                file_irreps.append(onept_outfilename(lline))
            if(lnmeas==2):
                break
    req_irreps=[]
    if(instring):
        for subs in instring.split(","):
            if(not(bool(re.match(r"[0-9]+_[0-9]+_[0-9]+_[0-9A-Za-z]+_[\+|\-]_[0-9]+_[0-9]+", subs)))):
                print("The Irrep Pattern ",subs," is not matched")
                print("The list of irreps available is:")
                for name in file_irreps:
                    print(name)
                sys.exit()
            else:
                tmps=subs.split("_")
                req_irreps.append("1pt_P=("+tmps[0]+","+tmps[1]+","+tmps[2]+")_Irrep="+tmps[3]+"_Charge="+tmps[4]+"_Irrep_ev="+tmps[5]+"x"+tmps[6]+".dat")
        tmp_irreps=list(set(file_irreps).intersection(req_irreps))
        if(len(tmp_irreps)==0):
            print("No irreps matching the requests have been found in ",linfile_name)
            print("The list of irreps available is:")
            for name in file_irreps:
                print(name)
            sys.exit(0)
        if(len(tmp_irreps)!=len(req_irreps)):
            print("Not all the requested irreps have been found in ",linfile_name)
            print("The list of irreps available is:")
            for name in file_irreps:
                print(name)
            sys.exit(0)
        file_irreps=tmp_irreps
    return file_irreps
 

def find_timelist(linfile_name):
    lnmeas=0
    lset=set()
    with open(linfile_name) as lfile:
        for lline in lfile:
            if("[MAIN][0]ML Measure #" in lline and "..." in lline):
                lnmeas+=1
            elif("[Measure ML][0] t=" in lline):
                lset.add(lline.translate(trantab4).split()[3])
            if(lnmeas==2):
                break
    lret=sorted([ int(el) for el in lset ])
    return [ str(el) for el in lret ]

def find_n_up_x_lev(linfile_name):
    lnuplev=[]
    with open(linfile_name) as lfile:
        for lline in lfile:
            if("[INIT ML][0]lev " in lline):
                lnuplev.append(re.split(' |=', lline)[4]+" ")
            elif("[MAIN][0]ML Measure #" in lline and "..." in lline):
                break
    return lnuplev

def filter_set_OP(dataset, l):
    row = np.in1d(dataset[:, 1], l)
    return dataset[row, :]

def define_verify_operator_list(lin_op_string,linfile_name,lirrep):
    ltmp=set()
    for op in lin_op_string.split(","):
        try:
            ltmp.add(int(op))
            if(int(op)<0):
              raise()
        except:
            print("The operator cut string must be composed of comma separated positive distinct integers")
            sys.exit()
    if(len(ltmp)!=len(lin_op_string.split(","))):
        print("The operator cut string must be composed of comma separated positive distinct integers")
        sys.exit()
    ltmp=list(ltmp)
    ltmp.sort()
    lastop=ltmp[-1]
    lnmeas=0
    lp=lirrep.split("_")
    string="[Measure ML][0]1pt function P=("+lp[0]+","+lp[1]+","+lp[2]+") Irrep="+lp[3]+" Irrep ev="+lp[5]+"/"+lp[6]+" Charge="+lp[4]
    with open(linfile_name) as lfile:
        for lline in lfile:
            if(string in lline):
                if(int(lline.split("=")[-1])<=lastop):
                    print("The operator cut string contains a index larger than the number of available operators")
                    sys.exit()
                break
    lret=[0]
    for i in ltmp:
        lret.append(2*i+1)
        lret.append(2*i+2)
    return lret


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True,
                    help='input measure file')
    parser.add_argument('-p', type=str, required=False,
                    help='prefix directory for the output file')
    parser.add_argument('-b', type=str, required=False,
                    help='full path to the bs_stat analysis code')
    parser.add_argument('-s', required=False, action='store_true',
                    help='perform the striping of the 1pt functions')
    parser.add_argument('-q', required=False, type=str,
                    help='select the irreps to be analized (e.g. px_py_pz_Irrep_Charge_rho_irrepsize,... -> 0_0_0_A1plusOhP_+_1_1,0_0_0_EplusOhP_-_1_2 ) note that the format is fixed')
    parser.add_argument('-y', required=False, action='store_true',
                    help='evaluate the 1pt analysis')
    parser.add_argument('-c', required=False, action='store_true',
                    help='evaluate all the (subtracted) correlators defined in the simulation output file, force enable the y argument')
    parser.add_argument('-a', required=False, action='store_true',
                    help='evaluate the 2pt (non GEVP) analysis, force enable the c argument\n')
    parser.add_argument('-k', required=False, type=str,
                    help='select only the specified operators for striping (csv starting from 0), must always be used in conjuction to -q and only one irrep\n')
    parser.add_argument('-B', required=False, type=int,
                    help='bin the data, to give B final bins\n')

    args = parser.parse_args()
    inputdata_filename=args.i

    if(args.B):
        args.c=True

    if(args.c):
        args.y=True

    if(args.i and not(os.path.isfile(inputdata_filename))):
        sys.stderr.write("Missing file "+inputdata_filename+"\n")
        sys.exit()

    if(args.k):
        args.s=True

    if(args.k and not(args.q)):
        sys.stderr.write("When using the -k option you must use also the -q\n")
        sys.exit()

    if(args.k and len(args.q.split(","))!=1):
        sys.stderr.write("When using the -k and -q options you must specify only one irrep\n")
        sys.exit()

    if(args.k):
        args.s=True

    if(args.p):
        prepend=args.p
    else:
        prepend=""
    if(len(prepend)>0):
        if(not(os.path.isdir(prepend))):
            os.makedirs(prepend)
            print("Directory" , prepend ,  "created.") 
        if(prepend[-1]!="/"):
            prepend+="/"
    if(args.b or args.y or args.c):
        bstat="./bstat"
        if(args.b):
            bstat=args.b
        if(os.path.isfile(bstat)):
            if(not(os.access(bstat, os.X_OK))):
                sys.stderr.write("File "+bstat+" is not executable\n")
                sys.exit()  
        else:
            sys.stderr.write("Missing executable file for the statistica analisys (bstat)\n")
            sys.exit()
        

    if(args.i):
        print("Reading general setup from",inputdata_filename)
        separator='_'
        myopenfile={}
        tset="start"
        nuplev=[]
        tlist=[]
        outfilelist=[]
        activemeas={}
        irreplist=[]

        tlist=find_timelist(inputdata_filename)
        nuplev=find_n_up_x_lev(inputdata_filename)
        outfilelist=find_outfile_list(args.q,inputdata_filename)
        if(args.k):
            opcutstring=define_verify_operator_list(args.k,inputdata_filename,args.q)


    if(args.s):
        print("Striping the 1pt functions from",inputdata_filename)
        writedata=False
        countmeas=0
        with open(inputdata_filename) as file:
            for line in file:
                if("[MAIN][0]ML Measure #" in line and "..." in line):
                    nmeas=int((line.split("Measure #",1)[1]).split("...",1)[0])
                    countmeas+=1
                elif("[Measure ML][0]1pt function" in line):
                    if(not(line in activemeas)):
                        outfile=onept_outfilename(line)
                        if(outfile in outfilelist):
                            myopenfile[line]=open(outfile,'w')
                            activemeas[line]=True
                        else:
                            activemeas[line]=False
                    if(activemeas[line]):
                        filewrite=myopenfile[line]
                        writedata=True
                    else:
                        writedata= False
                elif( writedata and "[Measure ML][0] t=" in line):
                    if(args.k):
                        ldata=(line.split("[Measure ML][0] t=",1)[1]).translate(trantab).split()
                        separator = ' '
                        data=separator.join([ ldata[i] for i in opcutstring ])+"\n"
                    else:
                        data=(line.split("[Measure ML][0] t=",1)[1]).translate(trantab)
                    data=str(nmeas)+" "+data
                    filewrite.write(data)
        file.close()
        for string in myopenfile:
            myopenfile[string].close()
        print("Found",countmeas,"independent measures")
    if(args.y):
        print("Performing the 1pt statistical analyis")


        for filename in outfilelist: 
            if(not(os.path.isfile(filename))):
                print("Missing the file ",filename," for the analysis, you might want to run the striping (-s) also")
                sys.exit()
        
        command = bstat.split()
        separator = ' '
        nuplev=separator.join(nuplev)
    
        tmpfilename=prepend+"tmpfile"
        for datafile in outfilelist:
            if "P=(0,0,0)_Irrep=A1plusOhP_Charge=+" in datafile:
                subtract=True
            else:
                subtract=False
            with open(datafile) as file:
                datataufilename=datafile.replace(".dat", "_autocor.dat")
                dataavgfilename=datafile.replace(".dat", "_avg.dat")
                datavarfilename=datafile.replace(".dat", "_var.dat")
                datatunefilename=datafile.replace(".dat", "_tune.dat")

                datataufile=open(datataufilename,'w')
                dataavgfile=open(dataavgfilename,'w')
                datavarfile=open(datavarfilename,'w')
                datatunefile=open(datatunefilename,'w')


                for ti in tlist:
                    file.seek(0)         
                    tmpfile=open(tmpfilename,'w')
                    for line in file:
                        if(line.split()[1]==ti):
                            tmpfile.write(line)
                    tmpfile.close()
                    tmpfile=open(tmpfilename,'r')
                    nfields=len(tmpfile.readline().split())
                    outlinetau=[ti]
                    outlinevar=[ti]
                    outlineavr=[ti]
                    outlinetune=[ti,nuplev]
                    nrescale=1
                    for i in nuplev.split():
                        nrescale*=int(i)
                    size=0
                    for i in range(2, nfields):
                        tmpfile.seek(0)
                        process = Popen(command, text=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
                        for tmpline in tmpfile:
                            if(i==2):
                                size +=1
                            process.stdin.write(tmpline.split()[i]+ '\n')
                        outstring=re.split('\n|:', process.communicate()[0])
                        outlinetau.append(outstring[1])
                        outlineavr.append(outstring[3])
                        outlinevar.append(outstring[5])
                    size*=nrescale
                    for i in range(2, nfields, 2):
                        tmpfile.seek(0)
                        process = Popen(command, text=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
                        for tmpline in tmpfile:
                            if(subtract):
                                subs=float(outlineavr[i-1].split()[0])
                                cdata=complex(float(tmpline.split()[i])-subs,float(tmpline.split()[i+1]))
                            else:
                                cdata=complex(float(tmpline.split()[i]),float(tmpline.split()[i+1]))
                            rdata=size*cdata.conjugate()*cdata
                            process.stdin.write(str(rdata.real)+ '\n')
                        outstring=re.split('\n|:', process.communicate()[0])
                        outlinetune.append(outstring[3])
                        

                    datataufile.write(separator.join(outlinetau)+'\n')
                    dataavgfile.write(separator.join(outlineavr)+'\n')
                    datavarfile.write(separator.join(outlinevar)+'\n')
                    datatunefile.write(separator.join(outlinetune)+'\n')
                    tmpfile.close()
                    
                datataufile.close()
                dataavgfile.close()
                datavarfile.close()
                datatunefile.close()
        os.remove(tmpfilename)


    if(args.c):
        print("Creating the 2pt correlators")
        correlator_list={}
        correlator_def(inputdata_filename,correlator_list)
 
        for filename in outfilelist:
            subs={}
            if "P=(0,0,0)_Irrep=A1plusOhP_Charge=+" in filename:
                subfilename=filename.replace(".dat","_avg.dat")
                if(not(os.path.isfile(subfilename))):
                    sys.stderr.write("Missing the average 1pt function file "+subfilename+" \n")
                    sys.exit()
                with open(subfilename) as subfile:
                   for line in subfile:
                       tmp = line.replace("+-","").split()
                       tmpkey=tmp[0]
                       del tmp[0::2]
                       tmp = np.array(tmp).astype(np.float)
                       tmp1=tmp[0::2]+1j*tmp[1::2]
                       subs[tmpkey]=tmp1
                subfile.close
            if(not(os.path.isfile(filename))):
                sys.stderr.write("Missing the 1pt function file "+filename+" \n")
                sys.exit()
            blockcounter=0
            n_meas=""
            corrdata={}
            corfilename=filename.replace("1pt_P","2pt_P")
            outfile=open(corfilename,'w')
            with open(filename) as datafile:
                for line in datafile:
                    if(line[0]=="#"):
                        continue
                    
                    tmp=line.split()

                    tmp1 = np.array(tmp[2:]).astype(np.float)
                    corrdata[tmp[1]]=tmp1[0::2]+1j*tmp1[1::2]

                    if(blockcounter % len(tlist) == len(tlist)-1 ):
                        n_meas=tmp[0]
                        correlator_meas(n_meas,corrdata,subs,correlator_list,outfile)

                        corrdata={}
                    blockcounter+=1                    
            outfile.close()
            print(corfilename)

            if(args.B):
                bnmeas=blockcounter/len(tlist)
                if(args.B>bnmeas):
                    sys.stderr.write("The number of requested bins "+str(args.B)+" is larger than the number of measures "+str(bnmeas)+"\n")
                    sys.exit()
                n_elem_per_bin=int(bnmeas/args.B)
                print("Binning the data in",args.B,"bins, each bin contains",n_elem_per_bin,"elements")
                unbindatacorfilename=corfilename.replace(".dat", "_unbinned.dat")
                os.rename(corfilename,unbindatacorfilename)
                databin=np.zeros(shape=(len(tlist), 4))
                lcounter=0
                with open(unbindatacorfilename) as unbindatafile:
                    line=unbindatafile.readline()
                    nfields=len(unbindatafile.readline().split())
                    unbindatafile.seek(0)
                    line=unbindatafile.readline()
                    corrdata={}
                    measid=int(line.split()[0])
                    bindatafile=open(corfilename,'w')

                    while(line):
                        dt=int(line.split()[1])
                        if(int(line.split()[0])!=measid):
                            lcounter+=1
                            measid=int(line.split()[0])
                        tmp=np.loadtxt(unbindatafile,max_rows=nfields)
                        #print(tmp)
                        if dt in corrdata:
                            corrdata[dt]=np.add(corrdata[dt],tmp)
                        else:
                            corrdata[dt]=tmp

                        if(lcounter % n_elem_per_bin == n_elem_per_bin-1):
                            databin=np.divide(corrdata[dt], n_elem_per_bin) 
                            bindatafile.write(str(measid)+" "+str(dt)+"\n")
                            np.savetxt(bindatafile, databin)
                            corrdata.pop(dt)
                        line=unbindatafile.readline()
                    bindatafile.close()
      
    if(args.a):
        print("Performing the 2pt statistical analyis")

        correlator_list={}
        correlator_def(inputdata_filename,correlator_list)

        dtlist=[]
        for id in correlator_list.keys():
            dtlist.append(abs(int(correlator_list[id][0][1])-int(correlator_list[id][0][0])))


        for filename in outfilelist:
            corfilename=filename.replace("1pt_P","2pt_P") 
            if(not(os.path.isfile(corfilename))):
                print("Missing the file ",corfilename," for the analysis, you might want to run the cor creation (-c) also")
                sys.exit()
 
        for filename in outfilelist:
            corfilename=filename.replace("1pt_P","2pt_P")
            with open(corfilename) as datafile:
                datafile.readline()
                nfields=len(datafile.readline().split())
                print("Analysying ",corfilename)
                res=np.array([])
                for dti in dtlist:
                    mat = []
                    for i in range(nfields* nfields):
                        mat.append([])

                    datafile.seek(0)
                    line=datafile.readline()
                    while line:
                        if(int(line.split()[1])!=dti):
                            for _ in range(nfields):
                                next(datafile) 
                        else:
                            tmp=np.loadtxt(datafile,max_rows=nfields)
                            for j in range(0,nfields):
                                for i in range(j,nfields):
                                    mat[i+nfields*j].append(tmp.item((i,j)))
                        line=datafile.readline()
                    for j in range(0,nfields):
                        for i in range(j,nfields):
                            res=np.append(res,np.array((dti,i+nfields*j,np.mean(mat[i+nfields*j]),np.std(mat[i+nfields*j])/np.sqrt(len(mat[i+nfields*j])))))
                res=res.reshape(int(len(res)/4),4)
                for j in range(0,nfields):
                    for i in range(j,nfields):
                        corfilename_avg=corfilename.replace(".dat","_avg_C"+str(i)+"_"+str(j)+".dat")
                        outfile_corij=open(corfilename_avg,'w')
                        np.savetxt(outfile_corij,filter_set_OP(res, [i+nfields*j]))
                        outfile_corij.close()