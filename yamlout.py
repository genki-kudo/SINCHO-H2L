#!/usr/bin/env python

from subprocess import run
import os
import yaml
import sys
from script_sincho.exec import *
from script_sincho.basic.basic_func import logo

bash=lambda x:run(x,shell=True)

if __name__ == '__main__':
    ###input args###
    conditions = str(sys.argv[1]) # condition.yaml

    ###load and preparation###
    with open(conditions,'r')as f:
        setting = yaml.safe_load(f)
    in_dir = setting['P2C_SINCHO']['working_directory']
    #####################################################
    ##### "TSUKUBA_OUTPUT" will be changed its name #####
    out_dir = setting['TSUKUBA_OUTPUT']['directory']
    #####################################################
    nums = setting['P2C_SINCHO']['num_of_parallel']


    for i in range(nums+1):
        n = str(i).zfill(3)
        if not os.path.exists(out_dir+"trajectory_"+n+"/"):
            os.makedirs(out_dir+"trajectory_"+n+"/")
        for k in ["prot_","lig_"]:
            lower_file = "/trajectory_"+n+"/"+k+n+".pdb"
            bash("cp "+in_dir+lower_file+" "+out_dir+lower_file)
    

    
    if list(setting['P2C_SINCHO']['output_method'].keys())[0]=="score_sort_evenly":
        num = int(list(setting['P2C_SINCHO']['output_method'].values())[0])
        for i in range(nums+1):
            n = str(i).zfill(3)
            yml = {}
            summary = {}
            flag = 0
            rank_index = 0
            count = 0
            for j in open(in_dir+"trajectory_"+n+"/sincho.log"):
                if "##### RESULTS #####" in j:
                    flag+=1
                if flag ==1 and "#" not in j and count<num:
                    count+=1
                    property = {}
                    property["atom_num"]=j.split()[0]
                    property["mw"]=float(j.split()[2])
                    ## will be updated
                    property["logp"]=float(j.split()[4])
                    property["acceptor"]={"min":0, "max":3}
                    property["donor"]={"min":0, "max":3}
                    ## will be updated
                    rank_index +=1
                    summary["rank"+str(rank_index).zfill(3)]=property
                
            yml['SINCHO_result']=summary
            with open(out_dir+"trajectory_"+n+"/sincho_result.yaml","w")as y:
                yaml.dump(yml, y, encoding='utf-8', allow_unicode=True)

        
        
        

