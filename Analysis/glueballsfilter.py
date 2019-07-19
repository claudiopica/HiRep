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
                sum += (1.0/(2*npoints)*(np.outer(np.conjugate(corr_data[pair[0]]), corr_data[pair[1]])+
                    np.outer(np.conjugate(corr_data[pair[1]]), corr_data[pair[0]])).real+
                    -0.0/(2*npoints)*(np.outer(sub_data[pair[0]],sub_data[pair[1]])
                    +np.outer(np.conjugate(sub_data[pair[1]]), sub_data[pair[0]])).real)
           
        file.write(str(nmeas)+" "+str(dt)+"\n")
        np.savetxt(file, sum)


def filter_set_OP(dataset, l):
    row = np.in1d(dataset[:, 1], l)
    return dataset[row, :]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True,
                    help='input measure file')
    parser.add_argument('-p', type=str, required=False,
                    help='prefix directory for the output file')
    parser.add_argument('-b', type=str, required=False,
                    help='full path to the bs_stat analysis code')
    parser.add_argument('-y', required=False, action='store_true',
                    help='evaluate the 1pt analysis')
    parser.add_argument('-c', required=False, action='store_true',
                    help='evaluate all the (subtracted) correlators defined in the simulation output file, force enable the y argument')
    parser.add_argument('-a', required=False, action='store_true',
                    help='evaluate the 2pt (non GEVP) analysis, force enable the c argument\n')

    args = parser.parse_args()
    inputdata_filename=args.i

    if(args.a):
        args.c=True

    if(args.c):
        args.y=True

    if(not(os.path.isfile(inputdata_filename))):
        sys.stderr.write("Missing file "+inputdata_filename+"\n")
        sys.exit()   
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

    print("Striping 1pt from",inputdata_filename)
    separator='_'
    mydict={}
    tset="start"
    nuplev=[]
    tlist=[]
    outfilelist=[]
    with open(inputdata_filename) as file:
        for line in file:
            if("[MAIN][0]ML Measure #" in line and "..." in line):
                nmeas=int((line.split("Measure #",1)[1]).split("...",1)[0])
            elif("[Measure ML][0]1pt function" in line):
                if(not(line in mydict)):
                    head=line.split("[Measure ML][0]1pt function ",1)[1]
                    head=separator.join(list(head.split(' ')[i] for i in [-6,-5,-2,-4,-3]))
                    outfile=prepend+("1pt_"+head+".dat").translate(trantab2)
                    outfilelist.append(outfile)
                    mydict[line]=open(outfile,'w')
                    mydict[line].write("# "+head+"\n")
                filewrite=mydict[line]
                if(tset=="start"):
                    tset="firstset"
                else:
                    tset="done"
            elif("[INIT ML][0]lev " in line):
                nuplev.append(re.split(' |=', line)[4]+" ")
            elif("[Measure ML][0] t=" in line):
                data=str(nmeas)+" "+(line.split("[Measure ML][0] t=",1)[1]).translate(trantab)
                filewrite.write(data)
                if(tset=="firstset"):
                    tlist.append(data.split()[1])
    file.close()
    for string in mydict:
        mydict[string].close()

    
    if(args.y):
        print("Performing the 1pt statistical analyis")

        command = bstat.split()
        separator = ' '
        nuplev=separator.join(nuplev)
    
        tmpfilename=prepend+"tmpfile"
        for datafile in outfilelist:
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

                    for i in range(2, nfields, 2):
                        tmpfile.seek(0)
                        process = Popen(command, text=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
                        for tmpline in tmpfile:
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

                    blockcounter+=1                    
                    tmp1 = np.array(tmp[2:]).astype(np.float)
                    corrdata[tmp[1]]=tmp1[0::2]+1j*tmp1[1::2]

                    if(blockcounter==len(tlist)):
                        n_meas=tmp[0]
                        correlator_meas(n_meas,corrdata,subs,correlator_list,outfile)

                        corrdata={}
                        blockcounter=0
            outfile.close()

        
    if(args.a):
        print("Performing the 2pt statistical analyis")
        dtlist=[]
        for id in correlator_list.keys():
            dtlist.append(abs(int(correlator_list[id][0][1])-int(correlator_list[id][0][0])))

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
                            #res=np.append(res,np.array((dt,int(i+nfields*j)),dtype='int,int'))
                            res=np.append(res,np.array((dti,i+nfields*j,np.mean(mat[i+nfields*j]),np.std(mat[i+nfields*j])/np.sqrt(len(mat[i+nfields*j])))))
                res=res.reshape(int(len(res)/4),4)
                for j in range(0,nfields):
                    for i in range(j,nfields):
                        corfilename_avg=corfilename.replace(".dat","_avg_C"+str(i)+"_"+str(j)+".dat")
                        outfile_corij=open(corfilename_avg,'w')
                        np.savetxt(outfile_corij,filter_set_OP(res, [i+nfields*j]))
                        outfile_corij.close()