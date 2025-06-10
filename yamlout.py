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
    in_dir = os.path.join(setting['OUTPUT']['directory'], setting['SINCHO']['working_directory'])
    out_dir = setting['OUTPUT']['directory']
    nums = setting['SINCHO']['num_of_parallel']
    order_scale = int(len(str(int(nums)))+1)

    if list(setting['SINCHO']['output_method'].keys())[0]=="score_sort_evenly":
        num = int(list(setting['SINCHO']['output_method'].values())[0])
        rank_order_scale = int(len(str(int(num)))+1)
        for i in range(nums+1):
            n = str(i).zfill(rank_order_scale)
            yml = {}
            summary = {}
            flag = 0
            rank_index = 0
            count = 0
            print(os.path.join(in_dir, 'trajectory_'+n, 'sincho.log'))
            for j in open(os.path.join(in_dir,"trajectory_"+n, "sincho.log")):
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
                    summary["rank_"+str(rank_index).zfill(rank_order_scale)]=property
                    rank_index +=1
                
            yml['SINCHO_result']=summary
            with open(os.path.join(in_dir,"trajectory_"+n,"sincho_result.yaml"),"w")as y:
                yaml.dump(yml, y, encoding='utf-8', allow_unicode=True)

        
        
        

