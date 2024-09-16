import importlib
import sys
from tqdm import tqdm
import numpy as np
import  os
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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

def event_3d_plot(proAnubis_event, reconstructor, title, save=False):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Coordinates for the horizontal plane at z = 60
    RPC_heights = [0.6, 1.8, 3.0, 61.8, 121.8, 123]
    distance_per_phi_channel = 2.7625 #cm
    distance_per_eta_channel = 2.9844 #cm
    dimensions = [distance_per_phi_channel*64, distance_per_eta_channel*32]
    x = np.array([0, dimensions[0]])
    y = np.array([0, dimensions[1]])
    X, Y = np.meshgrid(x, y)
    for rpc_h in RPC_heights:
        z_plane = np.ones((1, 1))*rpc_h  # Horizontal plane at z = 60
        ax.plot_surface(X, Y, z_plane, color='blue', alpha=0.1)

    # Plot the full horizontal plane (green)

    # make the hits
    hits = []
    for tdc in range(5):
        for word in proAnubis_event.tdcEvents[tdc].words:
            _, hit = ATools.tdcChanToRPCHit(word, tdc, 0)
            hits.append(hit)

    colours = ['orange', "red"]
    for hit in hits:
        col = colours[int(150 < hit.time < 350)]
        if hit.eta: #eta
            stripe_coord = [hit.channel*distance_per_eta_channel,(hit.channel+1)*distance_per_eta_channel+1]
            x_stripe = np.array([[0, 0], [dimensions[0], dimensions[0]]])    # Spans the y-axis fully
            y_stripe = np.array([stripe_coord, stripe_coord])  # Only a portion of the x-axis
            z_stripe = np.ones((1, 1)) * RPC_heights[hit.rpc]              # Same height (z=60)
            ax.plot_surface(x_stripe, y_stripe, z_stripe, color=col, alpha=0.9)
        else:
            stripe_coord = [hit.channel*distance_per_phi_channel,(hit.channel+1)*distance_per_phi_channel]
            x_stripe = np.array([stripe_coord, stripe_coord])  # Only a portion of the x-axis
            y_stripe = np.array([[0, 0], [dimensions[1], dimensions[1]]])    # Spans the y-axis fully
            z_stripe = np.ones((1, 1)) * RPC_heights[hit.rpc]              # Same height (z=60)
            ax.plot_surface(x_stripe, y_stripe, z_stripe, color=col, alpha=0.9)

    #reconstructor = RTools.Reconstructor([proAnubis_event], 0)
    
    reconstructor.update_event([proAnubis_event], 0)
    cluster = reconstructor.cluster()
    print("RPC:",sum([1 for x in cluster[0] if x]))
    tracks = reconstructor.reconstruct_track(cluster)
    print("Tracks:", len(tracks[0]))
    for track in tracks[0]:
        optimised_centroid = track.centroid[1:]
        optimised_d = track.direction[1:]
        x = []
        y = []
        z = []
        for rpc, height in enumerate(RPC_heights):
            t = (height-optimised_centroid[2])/optimised_d[2]
            x = np.append(x, optimised_centroid[0] + optimised_d[0]*t)
            y = np.append(y, optimised_centroid[1] + optimised_d[1]*t)
            z = np.append(z, height)
        y = [dimensions[1] - i for i in y]
        # Plot points with a red connecting line
        ax.scatter(x, y, z, color='blue', s=20)
        ax.plot(x, y, z, color='purple', linewidth=2)
    

    # Set axis labels
    ax.set_xlabel('phi channels')
    ax.set_ylabel('eta channels')
    ax.set_zlabel('Z / cm')


    ax.set_xticks([i*distance_per_phi_channel for i in range(0, 65, 8)], labels=[str(i) for i in range(0, 65, 8)])
    ax.set_yticks([i*distance_per_eta_channel for i in range(0, 33, 8)], labels=[str(i) for i in range(0, 33, 8)])
    # Set axis limits
    ax.set_xlim(0, dimensions[0])
    ax.set_ylim(0, dimensions[1])
    ax.set_title(title)
    ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1, dimensions[1]/dimensions[0], 1, 1]))
    ax.autoscale()
    ax.view_init(elev=30, azim=-130)
    ax.invert_yaxis()
    ax.set_zlim(0, 130)
    
    if save:
        plt.savefig("video//images_video//"+title+".png")
    else:
        plt.show()
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
            
    

