import json
import os

def correlation_value():
    mw_from_vol = json.loads(open(os.path.dirname(os.path.abspath(__file__))+'/vol_mw_fitting.json').read())
    mw_from_dist = json.loads(open(os.path.dirname(os.path.abspath(__file__))+'/dist_passmw_fitting.json').read())

    return mw_from_vol, mw_from_dist