import Analysis_tools as ATools
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


def tdcEventToRPCData(event,activeTDCs=[0,1,2,3,4], event_num = 0):
    rpcHits = [[] for rpc in range(6)]
    for tdc in activeTDCs:
        for word in event.tdcEvents[tdc].words:
            rpc, thisHit = ATools.tdcChanToRPCHit(word,tdc, event_num=event_num)
            rpcHits[rpc].append(thisHit)
    return rpcHits

def rpcHitToTdcChan(rpc, rpcChan, eta):
    tdc = -1
    tdcChannel = -1
    if rpc == 0:
        if eta:
            if rpcChan < 32:
                tdc = 0
                tdcChannel = rpcChan
            else:
                return None  # Invalid rpcChan for this combination
        else:
            if rpcChan < 64:
                tdc = 0
                tdcChannel = rpcChan + 32
            else:
                return None  # Invalid rpcChan for this combination
    elif rpc == 1:
        if eta:
            if rpcChan < 32:
                tdc = 0
                tdcChannel = rpcChan + 96
            else:
                return None  # Invalid rpcChan for this combination
        else:
            if rpcChan < 64:
                tdc = 1
                tdcChannel = rpcChan
            else:
                return None  # Invalid rpcChan for this combination
    elif rpc == 2:
        if eta:
            if rpcChan < 32:
                tdc = 1
                tdcChannel = rpcChan + 64
            else:
                return None  # Invalid rpcChan for this combination
        else:
            if rpcChan < 32:
                tdc = 2
                tdcChannel = rpcChan + 32
            elif rpcChan < 64:
                tdc = 1
                tdcChannel = rpcChan + 96
            else:
                return None  # Invalid rpcChan for this combination
    elif rpc == 3:
        if eta:
            if rpcChan < 32:
                tdc = 2
                tdcChannel = rpcChan + 32
            else:
                return None  # Invalid rpcChan for this combination
        else:
            if rpcChan < 64:
                tdc = 2
                tdcChannel = rpcChan + 64
            else:
                return None  # Invalid rpcChan for this combination
    elif rpc == 4:
        if eta:
            if rpcChan < 32:
                tdc = 3
                tdcChannel = rpcChan
            else:
                return None  # Invalid rpcChan for this combination
        else:
            if rpcChan < 64:
                tdc = 3
                tdcChannel = rpcChan + 32
            else:
                return None  # Invalid rpcChan for this combination
    elif rpc == 5:
        if eta:
            if rpcChan < 32:
                tdc = 3
                tdcChannel = rpcChan + 96
            else:
                return None  # Invalid rpcChan for this combination
        else:
            tdc = 4
            tdcChannel = rpcChan

    if tdc != -1 and tdcChannel != -1:
        # word = (tdcChannel << 24) & 0x7f000000
        # word |= rpcChan & 0xfffff
        return tdc, tdcChannel
    else:
        return -None

def plot_tdc_error_times(TDC_error_time):
    plt.figure(figsize=(12, 8))
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    bins = list(range(0, 1251, 50)) + [float('inf')]  
    text_offset = 0.1  

    for tdc in range(5):
        min_times = [entry[0] for entry in TDC_error_time[tdc]]
        counts, bin_edges = np.histogram(min_times, bins=bins)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        bin_centers[-1] = 1251
        plt.plot(bin_centers, counts, label=f'TDC {tdc}', linestyle='-', marker='o', color=colors[tdc])

        overflow_count = counts[-1]
        if overflow_count > 0:
            plt.text(1200, 100 * text_offset, f'overflow count == {overflow_count}', fontsize=12, 
                    ha='center', va='bottom', color=colors[tdc], bbox=dict(facecolor='white', alpha=0.6))
            text_offset *= 2

    plt.yscale('log')
    plt.xlabel('min_time')
    plt.ylabel('Events')
    plt.title('Histograms of min_time for all TDCs')
    plt.xlim(0, 1300)
    plt.legend()
    plt.grid(True)
    plt.show()
    
    
def plot_tdc_error_times_custom_ranges(TDC_error_time, ranges, output_pdf='TDC_error_times.pdf'):
    """
    Plot histograms of min_time for different TDCs over custom process time ranges and save the plots to a PDF.

    Parameters:
    TDC_error_time (list of list of tuples): The data for TDC error times. Each sublist corresponds to a TDC and contains tuples of (min_time, process_time).
    ranges (list of tuples): A list of tuples where each tuple defines the start and end of a range for process time (e.g., [(0, 15000), (25000, 40000)]).
    output_pdf (str): The name of the output PDF file where the plots will be saved.

    Usage Example:
    --------------
    # Define TDC_error_time as needed
    TDC_error_time = [
        [(100, 5000), (200, 12000), ...],  # TDC 0 data
        [(50, 3000), (150, 26000), ...],   # TDC 1 data
        ...
    ]

    # Define the custom ranges you want to plot
    ranges = [(0, 15000), (25000, 40000), (45000, 60000)]

    # Generate and save plots to a PDF
    plot_tdc_error_times_custom_ranges(TDC_error_time, ranges, output_pdf='TDC_error_times.pdf')
    """
    
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    bins = list(range(0, 1251, 50)) + [float('inf')]
    text_offset_base = 0.1

    with PdfPages(output_pdf) as pdf:
        for tdc in range(5):
            plt.figure(figsize=(12, 8))
            text_offset = text_offset_base
            for i, (range_start, range_end) in enumerate(ranges):
                min_times = []
                event_count = 0
                last_process = -1
                for entry, process in TDC_error_time[tdc]:
                    if last_process > process:
                        event_count += 2500
                    if range_start <= process+event_count < range_end:
                        min_times.append(entry[0])

                counts, bin_edges = np.histogram(min_times, bins=bins)
                bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
                bin_centers[-1] = 1251

                plt.plot(bin_centers, counts, label=f'TDC {tdc} ({range_start}-{range_end})', linestyle='-' if i % 2 == 0 else '--', marker='o' if i % 2 == 0 else 'x', color=colors[tdc])

                overflow_count = counts[-1]
                if overflow_count > 0:
                    plt.text(1200, 100 * text_offset, f'overflow count ({range_start}-{range_end}) == {overflow_count}', fontsize=12, 
                            ha='center', va='bottom', color=colors[tdc], bbox=dict(facecolor='white', alpha=0.6))
                    text_offset *= 2

            plt.yscale('log')
            plt.xlabel('min_time')
            plt.ylabel('Events')
            plt.title(f'Histograms of min_time for TDC {tdc}')
            plt.xlim(0, 1300)
            plt.legend()
            plt.grid(True)
            pdf.savefig()  # Save the current figure to the PDF
            plt.close()




def plot_tdc_error_channels(TDC_error_time, tdcs_to_plot=None):
    if tdcs_to_plot is None:
        tdcs_to_plot = range(5) 

    colors = ['blue', 'green', 'red', 'purple', 'orange']

    good_channels_dict = {tdc: [] for tdc in tdcs_to_plot}
    bad_channels_dict = {tdc: [] for tdc in tdcs_to_plot}

    for tdc in tdcs_to_plot:
        for (min_time, min_word), processedEvents in TDC_error_time[tdc]:
            channel = (min_word >> 24) & 0x7f
            if min_time > 300:
                bad_channels_dict[tdc].append(channel)
            elif min_time <= 300:
                good_channels_dict[tdc].append(channel)

    plt.figure(figsize=(24, 12))
    for tdc in tdcs_to_plot:
        if good_channels_dict[tdc]:
            plt.hist(good_channels_dict[tdc], bins=128, alpha=0.7, edgecolor=colors[tdc],
                    label=f'TDC {tdc}', histtype='step', linewidth=2, color=colors[tdc])

    plt.xlabel('Channel')
    plt.ylabel('Events')
    plt.title('Good Times Channels Histogram for Selected TDCs')
    plt.legend()
    plt.grid(True)
    plt.show()
    """
    plt.figure(figsize=(24, 12))
    for tdc in tdcs_to_plot:
        if bad_channels_dict[tdc]:
            plt.hist(bad_channels_dict[tdc], bins=128, alpha=0.7, edgecolor=colors[tdc],
                    label=f'TDC {tdc}', histtype='step', linewidth=2, color=colors[tdc])
    # plt.yscale('log')
    plt.xlabel('Channel')
    plt.ylabel('Events')
    plt.title('Bad Times Channels Histogram for Selected TDCs')
    plt.legend()
    plt.grid(True)
    plt.show()
    """
    
    
def plot_tdc_error_channels_custom_ranges(TDC_error_time, ranges, tdcs_to_plot=None, output_pdf='TDC_error_channels.pdf'):
    """
    Plot histograms of good and bad time channels for different TDCs over custom process time ranges and save the plots to a PDF.

    Parameters:
    TDC_error_time (list of list of tuples): The data for TDC error times. Each sublist corresponds to a TDC and contains tuples of (min_time, min_word), process_time.
    ranges (list of tuples): A list of tuples where each tuple defines the start and end of a range for process time (e.g., [(0, 15000), (25000, 40000)]).
    tdcs_to_plot (list of int, optional): A list of TDC indices to plot. Defaults to plotting all TDCs (0-4).
    output_pdf (str): The name of the output PDF file where the plots will be saved.

    Usage Example:
    --------------
    # Define TDC_error_time as needed
    TDC_error_time = [
        [((100, 0x1800000), 5000), ((200, 0x2400000), 12000)],  # TDC 0 data
        [((50, 0x3000000), 3000), ((150, 0x1C00000), 26000)],   # TDC 1 data
        ...
    ]

    # Define the custom ranges you want to plot
    ranges = [(0, 15000), (25000, 40000), (45000, 60000)]

    # Define the TDCs to plot
    tdcs_to_plot = [0, 1]

    # Generate and save plots to a PDF
    plot_tdc_error_channels_custom_ranges(TDC_error_time, ranges, tdcs_to_plot, output_pdf='TDC_error_channels.pdf')
    """
    if tdcs_to_plot is None:
        tdcs_to_plot = range(5)  # Default to plotting all TDCs

    colors = ['blue', 'green', 'red', 'purple', 'orange']

    with PdfPages(output_pdf) as pdf:
        for tdc in tdcs_to_plot:
            plt.figure(figsize=(24, 12))
            for i, (range_start, range_end) in enumerate(ranges):
                event_count = 0
                good_channels = [0 for channel in range(128)]
                bad_channels = [0 for channel in range(128)]
                last_process = -1
                for (hit_time, hit_word), process in TDC_error_time[tdc]: #TDC_error_time[tdc] should send me all the hits
                    if last_process > process:
                        event_count += 2500
                    if range_start <= process+event_count < range_end:
                        channel = (hit_word >> 24) & 0x7f
                        if 350 > hit_time > 150:
                            good_channels[channel] += 1
                        else:
                            bad_channels[channel] += 1
                    last_process = process
                length = (range_end - range_start)
                #print(good_channels_dict)
                # Plotting good channels
                plt.subplot(2, len(ranges), i + 1)
                x_channels = [channel for channel in range(128)]
                good_channels = [channel/length for channel in good_channels]
                bad_channels = [channel/length for channel in bad_channels] # normalizing per event
                
                if good_channels:
                    if tdc == 4:
                        plt.step(x_channels, good_channels,
                                label=f'TDC {tdc} ({range_start}-{range_end}) Good', linewidth=2, color=colors[tdc])
                    else:
                        plt.step(x_channels, good_channels,
                             label=f'TDC {tdc} ({range_start}-{range_end}) Good', linewidth=2, color=colors[tdc])
                plt.xlabel('Channel')
                plt.ylabel('Events')
                plt.title(f'Good Times Channels Histogram for TDC {tdc} ({range_start}-{range_end})')
                plt.legend()
                plt.grid(True)

                # Plotting bad channels
                plt.subplot(2, len(ranges), len(ranges) + i + 1)
                if bad_channels:
                    if tdc == 4:
                        plt.step(x_channels, bad_channels,
                             label=f'TDC {tdc} ({range_start}-{range_end}) Bad', linewidth=2, linestyle='--', color=colors[tdc])
                    else:
                        plt.step(x_channels, bad_channels,
                             label=f'TDC {tdc} ({range_start}-{range_end}) Bad', linewidth=2, linestyle='--', color=colors[tdc])
                plt.xlabel('Channel')
                plt.ylabel('Events')
                plt.title(f'Bad Times Channels Histogram for TDC {tdc} ({range_start}-{range_end})')
                plt.legend()
                plt.grid(True)

            plt.tight_layout()
            pdf.savefig()  # Save the current figure to the PDF
            plt.close()
            
              
        
            
            
        
def plot_tdc_alignment_channels_custom_ranges(TDC_alignment_time, ranges, tdcs_to_plot=None, output_pdf='TDC_alignment_channels.pdf'):
    """
    Plot histograms of alignment channels for different TDCs over custom process time ranges and save the plots to a PDF.

    Parameters:
    TDC_alignment_time (list of list of tuples): The data for TDC alignment times. Each sublist corresponds to a TDC and contains tuples of (min_time, min_word, processed_events).
    ranges (list of tuples): A list of tuples where each tuple defines the start and end of a range for processed events (e.g., [(0, 15000), (25000, 40000)]).
    tdcs_to_plot (list of int, optional): A list of TDC indices to plot. Defaults to plotting all TDCs (0-4).
    output_pdf (str): The name of the output PDF file where the plots will be saved.

    Usage Example:
    --------------
    # Define TDC_alignment_time as needed
    TDC_alignment_time = [
        [(100, (0, 45), 5000), (200, (0, 67), 12000)],  # TDC 0 data
        [(50, (0, 34), 3000), (150, (0, 89), 26000)],   # TDC 1 data
        ...
    ]

    # Define the custom ranges you want to plot
    ranges = [(0, 15000), (25000, 40000), (45000, 60000)]

    # Define the TDCs to plot
    tdcs_to_plot = [0, 1]

    # Generate and save plots to a PDF
    plot_tdc_alignment_channels_custom_ranges(TDC_alignment_time, ranges, tdcs_to_plot, output_pdf='TDC_alignment_channels.pdf')
    """
    if tdcs_to_plot is None:
        tdcs_to_plot = range(5)  # Default to plotting all TDCs

    colors = ['blue', 'green', 'red', 'purple', 'orange']
    bins = [i - 0.5 for i in range(129)]  # Creating bin edges to center bins on integer values

    with PdfPages(output_pdf) as pdf:
        for tdc in tdcs_to_plot:
            plt.figure(figsize=(24, 12))
            for i, (range_start, range_end) in enumerate(ranges):
                good_channels_dict = []
                bad_channels_dict = []

                for min_time, min_word, processed_events in TDC_alignment_time[tdc]:
                    if range_start <= processed_events < range_end:
                        channel = min_word[1]
                        if min_time > 300:
                            bad_channels_dict.append(channel)
                        else:
                            good_channels_dict.append(channel)

                plt.subplot(len(ranges), 1, i + 1)
                if good_channels_dict:
                    plt.hist(good_channels_dict, bins=bins, alpha=0.7, edgecolor=colors[tdc],
                             label=f'TDC {tdc} ({range_start}-{range_end})', histtype='step', linewidth=2, color=colors[tdc])
                if bad_channels_dict:
                    plt.hist(bad_channels_dict, bins=bins, alpha=0.7, edgecolor=colors[tdc],
                             label=f'TDC {tdc} ({range_start}-{range_end})', histtype='step', linewidth=2, linestyle='--', color=colors[tdc])
                plt.xlabel('Channel')
                plt.ylabel('Events')
                plt.title(f'Alignment Channels Histogram for TDC {tdc} ({range_start}-{range_end})')
                plt.legend()
                plt.xlim(0, 128)
                plt.grid(True)

            plt.tight_layout()
            pdf.savefig()
            plt.close()
            
            
            
def plot_tdc_alignment_times_custom_ranges(TDC_alignment_time, ranges, output_pdf='TDC_alignment_times.pdf'):
    """
    Plot histograms of min_time for different TDCs over custom process time ranges and save the plots to a PDF.

    Parameters:
    TDC_alignment_time (list of list of tuples): The data for TDC alignment times. Each sublist corresponds to a TDC and contains tuples of (min_time, min_word, processed_events).
    ranges (list of tuples): A list of tuples where each tuple defines the start and end of a range for processed events (e.g., [(0, 15000), (25000, 40000)]).
    output_pdf (str): The name of the output PDF file where the plots will be saved.

    Usage Example:
    --------------
    # Define TDC_alignment_time as needed
    TDC_alignment_time = [
        [(100, (0, 45), 5000), (200, (0, 67), 12000)],  # TDC 0 data
        [(50, (0, 34), 3000), (150, (0, 89), 26000)],   # TDC 1 data
        ...
    ]

    # Define the custom ranges you want to plot
    ranges = [(0, 15000), (25000, 40000), (45000, 60000)]

    # Generate and save plots to a PDF
    plot_tdc_alignment_times_custom_ranges(TDC_alignment_time, ranges, output_pdf='TDC_alignment_times.pdf')
    """
    
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    bins = list(range(0, 1251, 50)) + [float('inf')]
    text_offset_base = 0.1

    with PdfPages(output_pdf) as pdf:
        for tdc in range(5):
            plt.figure(figsize=(12, 8))
            text_offset = text_offset_base
            for i, (range_start, range_end) in enumerate(ranges):
                min_times = [entry[0] for entry in TDC_alignment_time[tdc] if range_start <= entry[2] < range_end]

                counts, bin_edges = np.histogram(min_times, bins=bins)
                bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
                bin_centers[-1] = 1251

                plt.plot(bin_centers, counts, label=f'TDC {tdc} ({range_start}-{range_end})', linestyle='-' if i % 2 == 0 else '--', marker='o' if i % 2 == 0 else 'x', color=colors[tdc])

                overflow_count = counts[-1]
                if overflow_count > 0:
                    plt.text(1200, 100 * text_offset, f'overflow count ({range_start}-{range_end}) == {overflow_count}', fontsize=12, 
                            ha='center', va='bottom', color=colors[tdc], bbox=dict(facecolor='white', alpha=0.6))
                    text_offset *= 2

            plt.yscale('log')
            plt.xlabel('min_time')
            plt.ylabel('Events')
            plt.title(f'Histograms of min_time for TDC {tdc}')
            plt.xlim(0, 1300)
            plt.legend()
            plt.grid(True)
            pdf.savefig()  # Save the current figure to the PDF
            plt.close()
            
            

def show_strip_time_info(outDict, ph, et, rpc):            
    fig, ax = plt.subplots(figsize=(10, 8))

    rpcTitles = {0:"Triplet Low",1: "Triplet Mid", 2:"Triplet Top", 3:"Singlet",4:"Doublet Low",5:"Doublet Top"}
    outDict['diffHists'][rpc][ph][et].plot1d(ax=ax, color="teal", lw=3, label='proANUBIS Data')

    ax.set_xlabel('$\eta$ Hit Time - $\phi$ Hit Time (ns)')
    yrange = ax.get_ylim()
    ax.text(-2, 0.8*yrange[1], rpcTitles[rpc]+", $\eta$="+str(et)+" $\phi$="+str(ph), fontsize=14,
                    verticalalignment='top')
    #ax.set_yscale('log')

    def gaus(x,a,x0,sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    # popt,pcov = curve_fit(gaus,readDict['diffHists'][rpc][ph][et].axes.centers[0],readDict['diffHists'][rpc][ph][et].values(),p0=[1,13,5])
    # fitX = [0.1*x-10 for x in range(400)]
    # plt.plot(fitX,gaus(fitX,*popt),'r:',label='Gaussian Fit')
    ax.set_xlim([-100,100])
    # ax.text(-2, 0.7*yrange[1], "Fit mean:"+str(round(popt[1],2))+", $\sigma$: "+str(round(popt[2],2)), fontsize=14,
    #                 verticalalignment='top')
    plt.legend()
    plt.show()
    
