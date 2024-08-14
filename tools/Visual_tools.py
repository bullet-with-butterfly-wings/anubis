import importlib
import sys
from tqdm import tqdm
import numpy as np
import  os
import glob
# Add the directories to the sys.path
dir_path = "C://Users//jony//Programming//Python//Anubis//anubis//" # insert your directory path
sys.path.append(dir_path + "Osiris//processing//python")
sys.path.append(dir_path + "Osiris//monitoring//python")
sys.path.append(dir_path + "tools")

import Analysis_tools as ATools
import proAnubis_Analysis_Tools
import Reconstruction_tools as RTools
import mplhep as hep
import Timing_tools as TTools
import rawFileReader

hep.style.use([hep.style.ATLAS])

# Specify the directory
data_list = sorted([f for f in os.listdir("data") if os.path.isfile(os.path.join("data", f))], reverse=True) ##all files in data directory sorted from the newest to the oldest

file_path = [dir_path+"//data//"+"proAnubis_240731_0217.raw"] #list(map(lambda p: dir_path+"//data//"+p, data_list))[:1] # insert your file
def time_events(eventChunk):
    for event in eventChunk:
        print(event.tdcEvents[0].time)

def abs_count(eventChunk):
    total = [0 for tdc in range(5)]
    for event in eventChunk:
        for tdc in range(5):
            for word in event.tdcEvents[tdc].words:
                bad_channels = [[32],[0,96],[64],[31,32],[0,1,2,3]]
                if word & 0xfffff < 300 and (word >> 24) & 0x7f not in bad_channels[tdc]:                    
                    total[tdc] += 1
    return total

def all_hits_event(event):
    buffer = []
    for tdc in range(5):
        for word in event.tdcEvents[tdc].words:
            _, hit = ATools.tdcChanToRPCHit(word, tdc, 1)
            buffer.append(hit)
    return buffer

def all_hits(eventChunk, tdc1, tdc2, event_num):
    for i, event in enumerate(eventChunk[:event_num]):
        print("Event: ", i)
        print("--------------------")
        print(f"TDC{tdc1}")
        for word in event.tdcEvents[tdc1].words:
            _, hit = ATools.tdcChanToRPCHit(word, tdc1, 1)
            print(hit)
        print(f"TDC{tdc2}")
        for word in event.tdcEvents[tdc2].words:
            _, hit = ATools.tdcChanToRPCHit(word, tdc2, 1)
            print(hit)

def metric_possible(eventChunk, tdc1, tdc2):
    total = 0
    good_1 = 0
    good_2 = 0
    for event in eventChunk:
        checker = [False, False, False, False]
        for word in event.tdcEvents[tdc1].words:
            _, hit1 = ATools.tdcChanToRPCHit(word, tdc1, 1)
            if hit1.time < 300:
                if hit1.eta:
                    checker[0] = True
                else:
                    checker[1] = True
                good_1 += 1
        for word in event.tdcEvents[tdc2].words:
            _, hit2 = ATools.tdcChanToRPCHit(word, tdc2, 1)
            if hit2.time < 300:
                if hit2.eta:
                    checker[2] = True
                else:
                    checker[3] = True
                good_2 += 1
        if all(checker):
            total += 1
    return total, good_1, good_2

def hitHeatMap(event): #actually returns 6 heatmaps, one for each rpc
    heatMaps = [np.zeros((32,64)) for rpc in range(6)]
    for tdc in range(5):
        for word in event.tdcEvents[tdc].words:
            _, hit = ATools.tdcChanToRPCHit(word, tdc, 0)
            if hit.time < 300:
                if hit.eta:
                    heatMaps[hit.rpc][hit.channel,:] += np.ones(64)
                else:
                    heatMaps[hit.rpc][:,hit.channel] += np.ones(32)
    return heatMaps
"""
def hitHeatMap(eventChunk, evt_num):
    heatMap = np.zeros((32,64))
    for evt_num, event in enumerate(eventChunk):
        for tdc in range(5):
            for word in event.tdcEvents[tdc].words:
                _, hit = ATools.tdcChanToRPCHit(word, tdc, evt_num)
""" 
                
    # Finally, recontruction is done using cluster information
"""
    for tdc in range(5):
        print(eventChunk[2].tdcEvents[tdc].time)
        for word in eventChunk[2].tdcEvents[tdc].words:
            _, hit = ATools.tdcChanToRPCHit(word, tdc, 1)
            print(hit)
    """
   # times_words = [(word & 0xfffff, word) for word in words if (word >> 24) & 0x7f not in bad_channels[tdc]]
            
    

