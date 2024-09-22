import importlib
import sys
from tqdm import tqdm
import numpy as np
import  os
import glob
import random
import importlib
import sys
from tqdm import tqdm
import  os
import cv2
import glob
from datetime import datetime 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.rc('text', usetex=True)
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
from itertools import chain
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

def produce_video_images(chunk, tag = ""):
    with tqdm(total=len(chunk)) as pbar:
        for evt_num, evt in enumerate(chunk):
        #[print(hit) for hit in VTools.all_hits_event(evt)]
            maps = hitHeatMap(evt)
            for rpc in range(6):
                plt.imshow(maps[rpc], interpolation='nearest')
                plt.title(f'{tag}: Event {evt_num} RPC {rpc}')
                plt.savefig(f'video/images_video/{tag+str(10*evt_num+rpc)}.png')
            pbar.update(1)

def compose_video(title = "video", rpc = [0,1,2,3,4,5]):
    image_folder = 'video/images_video'
    images = [img for img in os.listdir(image_folder) if img.endswith(tuple([str(i)+".png" for i in rpc]))]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape


    video_name = f'video/{title}.mp4'
    video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), 8, (width,height))

    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

    cv2.destroyAllWindows()
    video.release()

def event_3d_plot(proAnubis_event, title, save=False):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Configuration and constants
    RPC_heights = [0.6, 1.8, 3.0, 61.8, 121.8, 123]
    distance_per_phi_channel = 2.7625  # cm
    distance_per_eta_channel = 2.9844  # cm
    dimensions = [distance_per_phi_channel * 64, distance_per_eta_channel * 32]

    # Plot horizontal planes
    x = np.array([0, dimensions[0]])
    y = np.array([0, dimensions[1]])
    X, Y = np.meshgrid(x, y)
    for rpc_h in RPC_heights:
        z_plane = np.full_like(X, rpc_h)  # Horizontal plane at the current z
        ax.plot_surface(X, Y, z_plane, color='blue', alpha=0.1)

    # Collect hits
    hits = [ATools.tdcChanToRPCHit(word, tdc, 0)[1] 
            for tdc in range(5) 
            for word in proAnubis_event.tdcEvents[tdc].words]

    # Plot hits
    colours = ['orange', 'red']
    print_hits = ""
    for hit in hits:
        print_hits += str(hit) + "\n"
        col = colours[int(150 < hit.time < 350)]
        stripe_coord = [hit.channel * (distance_per_eta_channel if hit.eta else distance_per_phi_channel), 
                        (hit.channel + 1) * (distance_per_eta_channel if hit.eta else distance_per_phi_channel)]
        
        if hit.eta:
            # Eta stripe spans y-axis
            x_stripe = np.array([[0, 0], [dimensions[0], dimensions[0]]])
            y_stripe = np.array([stripe_coord, stripe_coord])
        else:
            # Phi stripe spans x-axis
            x_stripe = np.array([stripe_coord, stripe_coord])
            y_stripe = np.array([[0, 0], [dimensions[1], dimensions[1]]])

        z_stripe = np.full_like(x_stripe, RPC_heights[hit.rpc])
        ax.plot_surface(x_stripe, y_stripe, z_stripe, color=col, alpha=0.9)
    #ax.text(0, 0, 0, print_hits, color='black', fontsize=12, ha='right', va='top')
    print(print_hits)
    # Reconstruct tracks
    reconstructor = RTools.Reconstructor()
    reconstructor.update_event([proAnubis_event])
    all_possible_clusters = [(cluster.time[0], cluster.coords) for cluster in list(chain(*reconstructor.cluster()[0]))]
    tracks = reconstructor.reconstruct_tracks([proAnubis_event])[0]
    tracks_to_plot = []
    while tracks:
        track = tracks.pop(0)
        if isinstance(track, RTools.Vertex):
            if track.inside: #inside
                tracks_to_plot.append(([RPC_heights[0], track.point[2]], [track.initial.centroid, track.initial.direction], track.initial.clusters))
                tracks_to_plot.append(([track.point[2], RPC_heights[-1]], [track.final[0].centroid, track.final[0].direction], track.final[0].clusters))
                tracks_to_plot.append(([track.point[2], RPC_heights[-1] ], [track.final[1].centroid, track.final[1].direction], track.final[1].clusters))
            else: #double
                tracks_to_plot = [([track.point[2], RPC_heights[-1]], [track.final[0].centroid, track.final[0].direction], track.final[0].clusters)]
                tracks_to_plot.append(([track.point[2], RPC_heights[-1]], [track.final[1].centroid, track.final[1].direction], track.final[1].clusters))
        else: #single
            tracks_to_plot.append(([RPC_heights[0],RPC_heights[-1]], [track.centroid, track.direction], track.clusters))
        #plot clusters as well
    # Plot tracks
    for track_data in tracks_to_plot:
        x, y, z = [], [], []
        centroid, direction = track_data[1]
        bottom, top = track_data[0]
        height_to_plot = [bottom, *[h for h in RPC_heights if bottom < h < top], top]
        for h in height_to_plot:
            t = (h - centroid[3]) / direction[3]
            x.append(centroid[1] +  direction[1] * t)
            y.append(dimensions[1] - (centroid[2] + direction[2] * t))
            z.append(h)
        ax.scatter(x, y, z, color='blue', s=20)
        ax.plot(x, y, z, color='purple', linewidth=2)
        
        triplet = []
        singlet = []
        doublet = []

        for cluster in track_data[2]:
            if cluster.rpc < 3:
                triplet.append(cluster)
            elif cluster.rpc == 3:
                singlet.append(cluster)
            elif cluster.rpc > 3:
                doublet.append(cluster)

        for l in [triplet, singlet, doublet]:
            times = []
            for cluster in l:
                all_possible_clusters.remove((cluster.time[0], cluster.coords))
                times.append(round(float(cluster.time[0]), 1))
                ax.scatter(cluster.coords[0], dimensions[1] - cluster.coords[1], cluster.coords[2], color='green', s=50, marker = 'x')
            if times:
              ax.text(cluster.coords[0], dimensions[1] - cluster.coords[1], cluster.coords[2] + 4,f"{times}", size=10, zorder=1, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
    # Plot remaining clusters
    if all_possible_clusters:
        for cluster in all_possible_clusters:
            time, coords = cluster
            ax.scatter(coords[0], dimensions[1] - coords[1], coords[2], color='orange', s=50, marker = 'x')
            #ax.text(coords[0], dimensions[1] - coords[1], coords[2],f"{time}", size=10, zorder=1, bbox=dict(facecolor='orange', alpha=0.5, edgecolor='none'))
    
    # Set axis labels, ticks, and limits
    ax.set_xlabel('$\phi$ channels')
    ax.set_ylabel('$\eta$ channels')
    ax.set_zlabel('Z [cm]')
    ax.set_xticks([i * distance_per_phi_channel for i in range(0, 65, 8)], labels=[str(i) for i in range(0, 65, 8)])
    ax.set_yticks([i * distance_per_eta_channel for i in range(0, 33, 8)], labels=[str(i) for i in range(0, 33, 8)])
    ax.set_xlim(0, dimensions[0])
    ax.set_ylim(0, dimensions[1])
    ax.set_zlim(0, 130)
    ax.set_title(title)

    # Adjust projection for axis scaling
    ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1, dimensions[1] / dimensions[0], 1, 1]))
    ax.invert_yaxis()
    ax.autoscale()
    ax.view_init(elev=30, azim=-130)

    # Save or show the plot
    if save:
        plt.savefig(f"video//images_video//{title}.png")
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
            
    

