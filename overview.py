import importlib
import sys
from tqdm import tqdm
import  os
import pickle
import glob

dir_path = "C://Users//jony//Programming//Python//Anubis//anubis//" # insert your directory path
sys.path.append(dir_path + "Osiris//processing//python")
sys.path.append(dir_path + "tools")

from matplotlib.backends.backend_pdf import PdfPages
import Analysis_tools as ATools
import proAnubis_Analysis_Tools
import Reconstruction_tools as RTools
import mplhep as hep
import Timing_tools as TTools
import rawFileReader
from datetime import datetime
hep.style.use([hep.style.ATLAS])
import matplotlib.pyplot as plt
import numpy as np

interval = 100 # Set your monitoring chunck size
order = [[0,1], [1,2], [2,3], [3,4]] # Order what you want to align


def get_chunks(file_name, max_process_event = 20_000, fReader = None, start = None, end = None, tdc5 = False, saves = False):
    if not fReader:
        fReader = rawFileReader.fileReader(dir_path+"data//"+file_name) # load in the classs object    
    processedEvents = 0 # Initialisation
    max_process_event_chunk = max_process_event//interval
    last_reset = 0
    time = []
    chunks = []
    #tdc_mets = [[] for tdc in range(5)]
    #Tot_TDC_info = [[] for tdc in range(5)]#
    initial_chunk = fReader.get_aligned_events(order=order, interval=interval, extract_tdc_mets = False)
    initial_time = max([initial_chunk[0].tdcEvents[tdc].time for tdc in range(5) if initial_chunk[0].tdcEvents[tdc].time])
    print("Initial time", initial_time, datetime.timestamp(initial_time))
    if start:#if it is string, its a date
        if type(start) == str:
            goal = datetime.strptime(start, '%Y-%m-%d %H:%M:%S')
            event_time = initial_time
            with tqdm(total=round((goal-event_time).total_seconds()), desc=f"Skipping Events {file_name}", unit='Events') as pbar:        
                while event_time < datetime.strptime(start, '%Y-%m-%d %H:%M:%S'):
                    fReader.skip_events(2_000)
                    try:
                        chunk = fReader.get_aligned_events(order=order, interval=interval, extract_tdc_mets = False)
                    except Exception as e:
                        print(e)
                        return chunks, time, fReader
                        
                    if chunk:
                        event_time = max([chunk[0].tdcEvents[tdc].time for tdc in range(5) if chunk[0].tdcEvents[tdc].time])
                        pbar.update(round(event_time.timestamp()-pbar.n - initial_time.timestamp()))
        else:  
            with tqdm(total=start, desc=f"Skipping Events {file_name}", unit='Events') as pbar:
                event_counter = 0
                while event_counter < start:
                    fReader.skip_events(2_000)
                    event_counter += 2_000
                    pbar.update(2_000)

        print("Skip time", event_time)          
    if end:
        total_limit = (datetime.strptime(end, '%Y-%m-%d %H:%M:%S') - datetime.strptime(start, '%Y-%m-%d %H:%M:%S')).total_seconds()
        unit = "seconds"
    else:
        total_limit = max_process_event_chunk
        unit = 'Chunks'
        
    running = True
    if tdc5:
        fReader.evtBuilder.tdcFiveBuffer = []
        data = []
    counter = 0
    with tqdm(total=total_limit, desc=f"Processing Chunks {file_name}", unit=unit) as pbar:
            while running:
                processedEvents += 1
                try:
                    event_chunk, tdc_met, TDC_info = fReader.get_aligned_events(order=order, interval=interval, extract_tdc_mets = True) # get the aligned events
                except Exception as e:
                    print(e)
                    max_process_event_chunk = processedEvents
                    break

                event_time = max([event_chunk[0].tdcEvents[tdc].time for tdc in range(5) if event_chunk[0].tdcEvents[tdc].time])
                
                if last_reset:
                    if (event_time - last_reset).total_seconds() > 60:
                        last_reset = event_time
                        fReader.reload_event_builder()
                else:
                    last_reset = event_time
                    originTime = event_time
                    
                #[tdc_mets[i].append(tdc_met[i]) for i in range(5) if tdc_met[i] != 0]
                #[Tot_TDC_info[i].extend(TDC_info[i]) for i in range(5) if TDC_info[i]]
                time.append((event_time-originTime).total_seconds())
                chunks.append(event_chunk)
                counter += 1
                if tdc5:
                    tdc5_events = fReader.getTDCFiveEvents()
                    data.append([event_chunk, tdc5_events])

                if saves:
                    if counter % 5_000 == 0:
                        print("Event time:", event_time)
                        storage = "recovery.pkl"
                        with open(storage, "wb") as f:
                            pickle.dump([chunks, time], f)
                        print("Saved")

               
                if end:
                    pbar.update(round((event_time - datetime.strptime(start, '%Y-%m-%d %H:%M:%S')).total_seconds() - pbar.n, 2))
                else:     
                    pbar.update(1)
                if end:
                    running = event_time < datetime.strptime(end, '%Y-%m-%d %H:%M:%S')
                else:
                    running = processedEvents < max_process_event_chunk
    pbar.close()
    print("Ending time:" , event_time)
    print("Number of chunks:", len(chunks))
    if tdc5:
        return data, time, fReader
    return chunks, time, fReader

def alignment(chunks, times = None, pdf = None):
    mets = [[] for tdc in range(5)]
    for chunk in chunks:
        for idx, (i, j) in enumerate(order):
            x, y, l, m = ATools.find_tdc_alignment_metric(i, j) # determining which RPCs to use for aignment metric
            alignMet = ATools.calcAvgAlign(chunk, 0, x, y, l, m, i, j, 0) # i do not think the processed events is needed here
            
            mets[idx].append(alignMet) # write to memory
    fig, ax = plt.subplots(figsize=(10, 8))
    for idx, item in enumerate(order):
        met = mets[idx]
        i, j = item
        if times:
            binsx = times[:len(met)]
        else:
            binsx = [x*interval for x in range(len(met))]
        ax.plot(binsx, met, label=f'TDC{i} and TDC{j}, offset 0')


    ax.set_xlim(0, len(chunks)*interval)
    ax.set_ylim(-1, 40)
    ax.legend()
    ax.set_title('Alignment graph')
    ax.set_ylabel('Average $\sqrt{d\eta^2+d\phi^2}$')
    if times:
        ax.set_xlabel('Time (s)')
    else:
        ax.set_xlabel('Processed Events')
    
    if pdf:
        pdf.savefig()
    else:
        plt.show()
    plt.close()
    print("Alignment Done")
    return mets
    
def bvg(mets, window_size = 100, times = None, pdf = None):
    fig, ax = plt.subplots(figsize=(10, 8))    
    mets = np.convolve(mets, np.ones(np.ones(window_size)/window_size), mode="valid")
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    for tdc in range(5):
        met = mets[tdc]
        if times:
            ax.set_xlim(0,times[-1])
            ax.set_xlabel('Time (s)')

        else:
            binsx = [x for x in range(len(met))] #it might be wrong
            ax.set_xlim(0,binsx[-1])
            ax.set_xlabel('Processed Events')

        ax.plot(binsx, met, label = f'tdc {tdc}', color = colors[tdc])
    
    ax.legend()
    ax.set_title(f'TDC monitoring metric')
    ax.set_ylabel('off-peak times fraction')
    if pdf:
        pdf.savefig()
    else:
        plt.show()
    print("BVG Done")
    plt.close()

def noise_time(mets, time = None):
    for idx, alg in enumerate(mets):
        for chunk_number, alignment_metric in enumerate(alg):
            if alignment_metric > 40:
                if time:
                    return (chunk_number, time[chunk_number])
                return chunk_number

    print("No noise found")
    return None

def time_vs_chunks(times, pdf = None):
    plt.figure(figsize=(10, 6))
    plt.plot(times, [i for i in range(len(times))], label="Time between events")
    plt.xlabel("Time (s)")
    plt.ylabel("Event Chunk Number")
    plt.title(f"Time vs events")
    if pdf:
        pdf.savefig()
    else:
        plt.show()
    plt.close()
    print("Time vs Chunk Done")

def abs_bvg_hits(chunks, times = None,  pdf = None, per_rpc = False): #actually better
    if per_rpc:
        hit_counts = [[[],[]] for rpc in range(6)] #(good, bad)
    else:
        hit_counts = [[[],[]] for tdc in range(5)] #(good, bad)

    for event_chunk in chunks:
        if per_rpc:
            temp = [[0,0] for rpc in range(6)]
        else:
            temp = [[0,0] for tdc in range(5)]
        for tdc in range(5):
            for event in event_chunk:
                for hit in event.tdcEvents[tdc].words:
                    time = (hit & 0xfffff) #*(25/32)
                    if per_rpc:
                        rpc, _ = ATools.tdcChanToRPCHit(hit, tdc, 0)
                        if 150 < time < 370:
                            temp[rpc][0] += 1
                        else:
                            temp[rpc][1] += 1 
                    else:
                        if 150 < time < 370:
                            temp[tdc][0] += 1
                        else:
                            temp[tdc][1] += 1 

        for i, (good, bad) in enumerate(temp):
            hit_counts[i][0].append(good)
            hit_counts[i][1].append(bad)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    if times:
        binsx = times[:len(hit_counts[0][0])]
        ax.set_xlabel('Time (s)')
    else:
        binsx = [x*interval for x in range(len(hit_counts[0][0]))]
        ax.set_xlabel('Processed Events')
    
    if per_rpc:
        for rpc in range(6):
            ax.plot(binsx, hit_counts[rpc][0], label=f'RPC{rpc} in-peak')
            ax.plot(binsx, hit_counts[rpc][1], label=f'RPC{rpc} off-peak')
    else:
        for tdc in range(5):
            ax.plot(binsx, hit_counts[tdc][0], label=f'TDC{tdc} in-peak')
            ax.plot(binsx, hit_counts[tdc][1], label=f'TDC{tdc} off-peak')
    ax.set_title('Hit counts')
    ax.set_ylabel('Absolute hit count in chunk')
    ax.legend()
    if pdf:
        pdf.savefig()
    else:
        plt.show()
    plt.close()
    print("Absolute BVG Done")
    return hit_counts

def efficiency(chunks, residual = False, pdf = None):
    TAnalyser = proAnubis_Analysis_Tools.Timing_Analyser(chunks[0], 0)
    if residual:
        for processedEvents, event_chunk in enumerate(chunks[1:]):
            TAnalyser.update_event(event_chunk, processedEvents)
            TAnalyser.readTDCTimeDiffs()
        
        residEta, residPhi = TAnalyser.Calculate_Residual_and_plot_TDC_Time_Diffs( 
                                                     pdf_filename='output/TDC_time_diffs.pdf', 
                                                     max_itr = 5)
    reconstructor = proAnubis_Analysis_Tools.Reconstructor(chunks[0], 0)
    for processedEvents, event_chunk in enumerate(chunks[1:]):

            reconstructor.update_event(event_chunk, processedEvents)
            # populate_hits turns TDC bit wise information into their corresponding strips
            reconstructor.populate_hits()
            # This is optionnal, and requires the residual of eta and phi
            if residual:
                reconstructor.apply_systematic_correction(residEta, residPhi)
            # make_cluster does temporal and spatial coincidence between the stips, and reconstruction is done
            cluster = reconstructor.make_cluster()
            print(cluster)
            reconstructor.reconstruct_and_extrapolate(cluster, only_second = True)
    plt.figure(figsize=(10, 6))
    max_eff = []
    for RPC in range(6):
        if reconstructor.possible_reconstructions[RPC] == 0:
            efficiency = [0 for x in reconstructor.successful_reconstructions[RPC]]
        else:
            efficiency = [x / reconstructor.possible_reconstructions[RPC] for x in reconstructor.successful_reconstructions[RPC]]
        max_eff.append(efficiency[-1])
        plt.plot(reconstructor.tol, efficiency, label=f'RPC {RPC}')

    plt.xlabel('Tolerance')
    plt.ylabel('Efficiency')
    plt.title('The efficiency of the reconstruction')
    plt.legend()
    plt.grid(True)
    if pdf:
        pdf.savefig()
    else:
        plt.show()
    plt.close()
    print("Efficiency Done")
    print("Possible reconstructions", reconstructor.possible_reconstructions)
    return max_eff

def hit_time_hist(chunks, ranges = None, per_rpc = False, pdf = None):
    result = []
    std = []
    if not ranges:
        ranges = [(0, len(chunks))]
    for i, (start, end) in enumerate(ranges):
        if per_rpc:
            histograms = [[] for rpc in range(6)]
        else:
            histograms = [[] for tdc in range(5)]

        over = 0
        total = 0
        maximum = 0
        for event_chunk in chunks[start:end]:
            if per_rpc:
                temp = [[] for rpc in range(6)]
            else:
                temp = [[] for tdc in range(5)]
            for tdc in range(5):
                for event in event_chunk:
                    for hit in event.tdcEvents[tdc].words:
                        time = (hit & 0xfffff) #*(25/32)
                        total += 1
                        if time > 1200:
                            over += 1
                            if maximum < time:
                                maximum = time
                        if per_rpc:
                            rpc, _ = ATools.tdcChanToRPCHit(hit, tdc, 0)
                            histograms[rpc].append(time) 
                        else:
                            histograms[tdc].append(time)
    
        fig, ax = plt.subplots(figsize=(10, 8))
        new_std = []
        new_hist = []
        if per_rpc:
            for rpc in range(6):
                new_std.append(np.std(histograms[rpc]))
                hist, binsx = np.histogram(histograms[rpc], bins=100, density=True)
                ax.plot(binsx[:-1], hist, label=f"RPC{rpc}" )
                new_hist.append((binsx, hist))
        else:
            for tdc in range(5):
                new_std.append(np.std(histograms[tdc]))
                hist, binsx = np.histogram(histograms[tdc], bins=100, density=True)
                ax.plot(binsx[:-1], hist, label=f"TDC{tdc}" )
                new_hist.append((binsx, hist))
        std.append(new_std)
        result.append(new_hist)
        ax.set_title("Hit time histogram")
        ax.set_ylabel('Normalized hit count')
        ax.set_xlabel('Time (25/32 ns)')
        ax.legend()
        if pdf:
            pdf.savefig()
        else:
            plt.show()
    plt.close()
    print("Percentage of hits >1200:",over/total)
    print("Hit Time Histogram Done")
    return result, std


def hit_channel_hist(time_hits, ranges = None, tdcs_to_plot=[0], pdf=None): 
    if not ranges:
        ranges = [(0, len(time_hits[0]))]
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    for tdc in tdcs_to_plot:
        plt.figure(figsize=(24, 12))
        for i, (range_start, range_end) in enumerate(ranges):
            chunk_index = range_start
            good_channels = [0 for channel in range(128)]
            bad_channels = [0 for channel in range(128)]
            for chunk_index in range(range_start, range_end):
                for (hit_time, hit_word) in time_hits[tdc][chunk_index]: #TDC_error_time[tdc] should send me all the hits
                    channel = (hit_word >> 24) & 0x7f
                    if 370 > hit_time > 200:
                        good_channels[channel] += 1
                    else:
                        bad_channels[channel] += 1
                        
            length = (range_end - range_start)
            
            # Plotting good channels
            plt.subplot(2, len(ranges), i + 1)
            x_channels = [channel for channel in range(128)]
            good_channels = [channel/length for channel in good_channels]
            bad_channels = [channel/length for channel in bad_channels] # normalizing per event
            
            if good_channels:
                if tdc == 4: #different number of channels
                    plt.step(x_channels, good_channels,
                            label=f'TDC {tdc} ({range_start}-{range_end}) In-peak', linewidth=2, color=colors[tdc])
                else:
                    plt.step(x_channels, good_channels,
                            label=f'TDC {tdc} ({range_start}-{range_end}) In-peak', linewidth=2, color=colors[tdc])
            plt.xlabel('Channel')
            plt.ylabel('Events')
            plt.title(f'In-peak Times Channels Histogram for TDC {tdc} ({range_start}-{range_end})')
            plt.legend()
            plt.grid(True)

            # Plotting bad channels
            plt.subplot(2, len(ranges), len(ranges) + i + 1)
            if bad_channels:
                if tdc == 4: #different number of channels
                    plt.step(x_channels, bad_channels,
                            label=f'TDC {tdc} ({range_start}-{range_end}) Off-peak', linewidth=2, linestyle='--', color=colors[tdc])
                else:
                    plt.step(x_channels, bad_channels,
                            label=f'TDC {tdc} ({range_start}-{range_end}) Off-peak', linewidth=2, linestyle='--', color=colors[tdc])
            plt.xlabel('Channel')
            plt.ylabel('Events')
            plt.title(f'Off-peak Times Channels Histogram for TDC {tdc} ({range_start}-{range_end})')
            plt.legend()
            plt.grid(True)
    plt.tight_layout()
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()
    print("Hit Channel Histogram Done")

def possible_alignment():
        pass
        """
        for i, (tdc1, tdc2) in enumerate(order):
            mets[i].append(VTools.metric_possible(event_chunk, tdc1, tdc2)[0])
        """
def tdc_times_stats(chunks, only_min = False):
    bad_channels = [[32],[0,96],[64],[31,32],[0,1,2,3]] 
    bvg = [[] for tdc in range(5)] #bad vs good
    time_hits = [[] for tdc in range(5)] #[(time, word),...] 
    for chunk in chunks:
        for tdc in range(5):
            bad_time_count = 0
            total_time_count = 0
            for i, event in enumerate(chunk):
                words = event.tdcEvents[tdc].words
                times_words = [(word & 0xfffff, word) for word in words if (word >> 24) & 0x7f not in bad_channels[tdc]]
                
                if times_words and only_min:
                        times_words = [min(times_words, key=lambda x: x[0])]
                
                time_hits[tdc].append(times_words)
                for hit in times_words:
                    total_time_count += 1
                    if not 200 < hit[0] <= 370:
                        bad_time_count += 1                                                
                        if only_min:
                            event.tdcEvents[tdc].qual = 0x10 #raising flag
            ratio = bad_time_count / total_time_count 
            bvg[tdc].append(ratio)

    return time_hits, bvg
    
def cluster_size(chunks, residual = False, seperate = False):
    importlib.reload(proAnubis_Analysis_Tools)
   
    phi_histograms = [[] for rpc in range(6)]
    eta_histograms = [[] for rpc in range(6)]

    TAnalyser = proAnubis_Analysis_Tools.Timing_Analyser(chunks[0], 0)
    if residual:
        for processedEvents, event_chunk in enumerate(chunks[1:]):
            TAnalyser.update_event(event_chunk, processedEvents)
            TAnalyser.readTDCTimeDiffs()
        
        residEta, residPhi = TAnalyser.Calculate_Residual_and_plot_TDC_Time_Diffs( 
                                                     pdf_filename='output/TDC_time_diffs.pdf', 
                                                     max_itr = 5)
    reconstructor = proAnubis_Analysis_Tools.Reconstructor(chunks[0], 0)
    for processedEvents, event_chunk in enumerate(chunks[1:]):
            
            reconstructor.update_event(event_chunk, processedEvents)
            # populate_hits turns TDC bit wise information into their corresponding strips
            reconstructor.populate_hits()
            for h in reconstructor.etaHits:
                h = sorted(h,key = lambda x: x.event_num)
                for j in h:
                    pass
                    #print(j)
            # This is optionnal, and requires the residual of eta and phi
            if residual:
                reconstructor.apply_systematic_correction(residEta, residPhi)
            # make_cluster does temporal and spatial coincidence between the stips, and reconstruction is done
            clusters = reconstructor.make_cluster()
            for evt in clusters:
                    for rpc in range(6):
                        for phi_cluster in evt[2][rpc][0]:
                            phi_histograms[rpc].append(len(phi_cluster)) 
                        for eta_cluster in evt[2][rpc][1]:
                            eta_histograms[rpc].append(len(eta_cluster))
                        
            #print("Clust")
    if not seperate:
        histograms = []
        for rpc in range(6):
            histograms.append([phi_histograms[rpc]]+[eta_histograms[rpc]]) # only 6 together
    else:
        histograms = phi_histograms + eta_histograms #first 6 phi, next 6 eta
    histograms = np.array(histograms, dtype=object)
    avg_cluster_size = [np.mean(rpc_clusters) for rpc_clusters in histograms]
    errors = []
    for rpc_clusters in histograms:
        if rpc_clusters:
            errors.append(np.std(rpc_clusters)/np.sqrt(len(rpc_clusters)))
        else:
            errors.append(100)
    #over the whole population
    print("Cluster Size Done")
    if seperate:    
        return (avg_cluster_size[:6],avg_cluster_size[6:] ), (errors[:6],errors[6:] ), (histograms[:6],histograms[6:])
    else:
        return avg_cluster_size, errors, histograms