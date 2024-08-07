import importlib
import sys
from tqdm import tqdm
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

def hitHeatMap(eventChunk, tdc):
    print(len(eventChunk))
    print(eventChunk[0].tdcEvents[tdc].words)
    last_rpc = -1
    for tdc in range(5):
        for word in eventChunk[2].tdcEvents[tdc].words:
            _, hit = ATools.tdcChanToRPCHit(word, tdc, 1)
            if last_rpc != hit.rpc:
                print(f"RPC{hit.rpc}")
                last_rpc = hit.rpc
            print(hit)
   # times_words = [(word & 0xfffff, word) for word in words if (word >> 24) & 0x7f not in bad_channels[tdc]]
            
    

