import json
import os
import pickle
import joblib

def correlation_value():
    mw_from_vol = json.loads(open(os.path.dirname(os.path.abspath(__file__))+'/make_model/vol_mw.json').read())
    mw_from_dist = json.loads(open(os.path.dirname(os.path.abspath(__file__))+'/make_model/dist_passmw.json').read())

    return mw_from_vol, mw_from_dist

    
def logp_rftree_load():
    #rf_tree = joblib.load(os.path.dirname(os.path.abspath(__file__))+'/logp_rf_tree.pkl')
    rf_tree = joblib.load(os.path.dirname(os.path.abspath(__file__))+'/make_model/logp_rf.pkl')
    return rf_tree