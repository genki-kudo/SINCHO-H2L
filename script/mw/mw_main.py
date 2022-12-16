#!/usr/bin/env python

from script.basic.basic_func import *
from script.atom_select.min_theta_dist import *
from script.comp_stat_dataset.stat_correlation import *
from script.mw.volume_term import *
from script.mw.distance_term import *
import logging

def mw_calc(ligand, poc_and_atom, clusterdir, outputdir, mwsd, logger):
    pocket = clusterdir+poc_and_atom[0]

    ###correct naoshite!
    correct_term = 0

    #stat-correlation-set
    mw_from_vol, mw_from_dist = correlation_value()
    #print(mw_from_vol, mw_from_dist)

    estimate_pocvol = vol_calc(outputdir, poc_and_atom[0], clusterdir)
    estimate_dist = dist_calc(ligand, poc_and_atom, clusterdir)

    idealmw = (mw_from_vol[0]*estimate_pocvol+mw_from_vol[1])+(mw_from_dist[0]*estimate_dist+mw_from_dist[1])+(float(mwsd))
    logger.info(poc_and_atom[0]+"_"+poc_and_atom[1]+" estimate-mw: "+ str(idealmw))

    

