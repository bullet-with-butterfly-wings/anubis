import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'GTK3Agg', etc.
import mplhep as hep
hep.style.use([hep.style.ATLAS])
import sys
# import ANUBIS_triggered_functions as ANT
import matplotlib.backends.backend_pdf
from itertools import product
import numpy as np
# from scipy.stats import normpip install pillow
sys.path.insert(1, 'Osiris Temp\processing\python')
import Analysis_tools as ATools
import Reconstruction_tools as RTools
from scipy.optimize import minimize 
import math
from itertools import product
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

time_range = (150,350)
RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.
tot_speed_factor = 0.15204446322001586 #(25/32) ns / channel
tof_offsets = [[0, 0], [np.float64(0.25790911824685864), np.float64(0.5818488944295358)], [np.float64(0.06942782469835455), np.float64(0.5634127422432911)], [np.float64(2.17698945272121), np.float64(2.5044894876283283)], [np.float64(6.619902168271083), np.float64(2.91238011053787)], [np.float64(6.199631794023697), np.float64(2.617768586336547)]]
std_proportional_to_cluster_size = False
distance_per_phi_channel = 2.7625 #cm
distance_per_eta_channel = 2.9844 #cm
#there was an idea that the uncertainty is proportional to the cluster size
#Maybe it is right, maybe it is not
     
dir_path = "C://Users//jony//Programming//Python//Anubis//anubis//data//"
with open(dir_path + "tot_mean.pkl", "rb") as f:
    tot_mean = pickle.load(f)
with open(dir_path + "tot_std.pkl", "rb") as f:
    tot_std = pickle.load(f)

class Tanalyser():
    def __init__(self):
        self.tot_hits = np.empty((6, 64, 32), dtype=object)
        # Fill each element with a new empty list
        for index, _ in np.ndenumerate(self.tot_hits):
            self.tot_hits[index] = np.array([13])
        self.tot_mean = np.full((6, 64, 32), 13.0)
        self.tot_std = np.full((6, 64, 32), 1.0)
        
    def calculate_tot(self, chunks):
        # filter out the shady tracks
        for chun_number, chunk in enumerate(chunks):
            reconstructor = RTools.Reconstructor()
            reconstructor.update_chunk(chunk)
            chunk_clusters = reconstructor.cluster()
            for evt_clusters in chunk_clusters:
                if sum([1 for rpc_clust in evt_clusters if len(rpc_clust) == 1]) == 6: #1 cluster per rpc
                    track = RTools.Track([rpc_clust[0] for rpc_clust in evt_clusters])
                    if track.fit():
                        for rpc in range(6):
                            cluster = evt_clusters[rpc][0]
                            phi_time = min([hit.time for hit in cluster.hits[0]])
                            eta_time = min([hit.time for hit in cluster.hits[1]])
                            histogram = self.tot_hits[rpc][cluster.channel[0]][cluster.channel[1]]
                            histogram = np.append(histogram, eta_time - phi_time)
        self.tot_mean = np.mean(self.tot_hits, axis=2)
        self.tot_std = np.std(self.tot_hits, axis=2)
        return self.tot_mean, self.tot_std