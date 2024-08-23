import importlib
import sys
from tqdm import tqdm
import  os
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
import datetime
hep.style.use([hep.style.ATLAS])
import matplotlib.pyplot as plt
import numpy as np


def complete(file_name, start = None, end = None):
    interval = 100 # Set your monitoring chunck size
    order = [[0,1], [1,2], [2,3], [3,4]] # Order what you want to align
    max_process_event_chunk = 20_000 # End the loop early
    fReader = rawFileReader.fileReader(dir_path+"output"+file_name) # load in the classs object    
    processedEvents = 0 # Initialisation
    noiseStart = 0
    firstNoise = False
    initial_event_chunk = fReader.get_aligned_events(order=order, interval=interval)
    last_reset = 0
    mets = [[] for tdc in range(5)]
    time = []
    tdc_mets = [[] for tdc in range(5)]
    Tot_TDC_info = [[] for tdc in range(5)]

    event_time = max([initial_event_chunk[0].tdcEvents[tdc].time for tdc in range(5) if initial_event_chunk[0].tdcEvents[tdc].time])
    print(event_time)
    if start:
        while event_time < datetime.strptime(start, '%Y-%m-%d %H:%M:%S'):
            fReader.skip_events(2_000)
            initial_event_chunk = fReader.get_aligned_events(order=order, interval=interval)
            if initial_event_chunk:
                event_time = max([initial_event_chunk[0].tdcEvents[tdc].time for tdc in range(5) if initial_event_chunk[0].tdcEvents[tdc].time])
        print(event_time)    

    with PdfPages(file_name) as pdf:
        print("Alignment Metric")
        with tqdm(total=max_process_event_chunk, desc=f"Processing Events {file_name}", unit='Events') as pbar:
            while processedEvents < max_process_event_chunk:
                processedEvents += 1
                try:
                    event_chunk, tdc_met, TDC_info = fReader.get_aligned_events(order=order, interval=interval, extract_tdc_mets = True) # get the aligned events
                except Exception as e:
                    print("Exception:", e)
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
                    
                [tdc_mets[i].append(tdc_met[i]) for i in range(5) if tdc_met[i] != 0]
                [Tot_TDC_info[i].extend(TDC_info[i]) for i in range(5) if TDC_info[i]]
                time.append((event_time-originTime).total_seconds())
                for idx, (i, j) in enumerate(order):
                    x, y, l, m = ATools.find_tdc_alignment_metric(i, j) # determining which RPCs to use for aignment metric
                    alignMet = ATools.calcAvgAlign(event_chunk, 0, x, y, l, m, i, j, processedEvents) # determine the alignment metric
                    #alignMet_off = ATools.calcAvgAlign(event_chunk, 15, x, y, l, m, i, j, processedEvents) # determine the alignment metric
                    
                    #mets_off[idx].append(alignMet_off) # write to memory
                    mets[idx].append(alignMet) # write to memory
                    if alignMet > 40 and not firstNoise:
                        for tdc in range(5):
                            if event_chunk[0].tdcEvents[tdc].time != 0:
                                noiseTime = event_chunk[0].tdcEvents[tdc].time
                                break
                        if noiseTime == 0:
                            print("No hits in the event???")
                        startTime =  noiseTime - originTime
                        print(startTime)
                        noiseStart = processedEvents*interval
                        print(noiseStart)
                        firstNoise = True
                pbar.update(1)

        fig, ax = plt.subplots(figsize=(10, 8))
        for idx, item in enumerate(order):
            met = mets[idx]
            i, j = item
            binsx = [x for x in range(len(met))]
            ax.plot(binsx, met, label=f'TDC{i} and TDC{j}, offset 0')


        ax.set_xlim(0, max_process_event_chunk)
        ax.set_ylim(-1, 40)
        ax.legend()
        ax.set_title('Alignment graph')
        ax.set_ylabel('Average $\sqrt{d\eta^2+d\phi^2}$')
        ax.set_xlabel('Processed Event chunks')
        pdf.savefig()  # Save the current figure to the PDF
        plt.close()
    
        colors = ['blue', 'green', 'red', 'purple', 'orange']
        fig, ax = plt.subplots(figsize=(10, 8))
        for tdc in range(5):
            met = tdc_mets[tdc]
            binsx = [x*26.5 for x in range(len(met))]
            ax.plot(time[::25], met, label = f'tdc {tdc}', color = colors[tdc])

        ax.set_xlim(0,time[-1])
        ax.legend()
        ax.set_title(f'TDC monitoring metric {file_name}')
        ax.set_ylabel('bad time behavior / nominal time behavior')
        ax.set_xlabel('Time (s)')
        pdf.savefig()  # Save the current figure to the PDF
        plt.close()
        if noiseStart:
            TTools.plot_tdc_error_channels_custom_ranges(Tot_TDC_info, [(0, noiseStart), (noiseStart, max_process_event_chunk*interval)], tdcs_to_plot=None, output_pdf=f'output//TDC_channel_perEv_{file_name}.pdf')
        else:
            TTools.plot_tdc_error_channels_custom_ranges(Tot_TDC_info, [(0, max_process_event_chunk*interval)], tdcs_to_plot=None, output_pdf=f'output//TDC_channel_perEv_{file_name}.pdf')
            