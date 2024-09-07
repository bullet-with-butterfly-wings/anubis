import importlib
import sys
from tqdm import tqdm
import  os
import pickle
import glob
from datetime import datetime
import matplotlib.pyplot as plt
# Add the directories to the sys.path
dir_path = "C://Users//jony//Programming//Python//Anubis//anubis//" # insert your directory path
sys.path.append(dir_path + "Osiris//processing//python")
sys.path.append(dir_path + "tools")

import Analysis_tools as ATools
import Atlas_tools as AtTools
import proAnubis_Analysis_Tools
import Reconstruction_tools as RTools
import numpy as np
import scipy
import mplhep as hep
import Timing_tools as TTools
import pickle
import overview
import rawFileReader
import Visual_tools as VTools
from os.path import exists
importlib.reload(AtTools)

hep.style.use([hep.style.ATLAS])

dir = "C://Users//jony//Programming//Python//Anubis//anubis//data//"

def generate_bcr_histogram(data, atlas = False,  storage_name = "bcr_histogram.pkl"):
    res = AtTools.BCRHistogram(data, plot = False, atlas = atlas) #0,...,3653
    with open(dir + storage_name, "wb") as f:
        pickle.dump(res, f)
    return res

def calculate_background(anubis_data = None):
    period = 3564
    result = []
    
    if exists(dir + "bcr_histogram.pkl"):
        with open(dir + "bcr_histogram.pkl", "rb") as f:
            a_org = pickle.load(f)
    else:
        if anubis_data is None:
            raise ValueError("No data and no file provided")
        a_org = generate_bcr_histogram(anubis_data)
    
    bins = list(range(0, period))
    with tqdm(total=period) as pbar:    
        for diff in range(0,period):
            a_shift = np.roll(a_org, diff)
            correction = scipy.integrate.simpson((a_org*a_shift)[:3564-diff], x = bins[:3564-diff])
            result.append(correction)
            pbar.update(1)
    return result

def subtract_background(data, background, storage_name = "background_subtracted.pkl"):
    result = data - background
    with open(dir + storage_name, "wb") as f:
        pickle.dump(result, f)
    plt.plot(result)
    plt.show()
    return result