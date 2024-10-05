import numpy as np
import mplhep as hep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('TkAgg')  # or 'Qt5Agg', 'GTK3Agg', etc.
hep.style.use([hep.style.ATLAS])
import sys
# import ANUBIS_triggered_functions as ANT
# from scipy.stats import normpip install pillow
sys.path.insert(1, 'Osiris Temp\processing\python')
from scipy.optimize import curve_fit
mpl.rcParams['text.usetex'] = False
import seaborn as sns
import pickle
import Reconstruction_tools as RTools
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
            self.tot_hits[index] = np.array([])
        self.tot_mean = np.full((6, 64, 32), 13.0)
        self.tot_std = np.full((6, 64, 32), 1.0)
        self.offsets = np.full((6, 32), 13.0)
        self.speed = 0.0
        
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
                            self.tot_hits[rpc][cluster.channel[0]][cluster.channel[1]] = np.append(self.tot_hits[rpc][cluster.channel[0]][cluster.channel[1]], eta_time - phi_time)
        self.calculate_mean_std()
        return self.tot_mean, self.tot_std
    
    def calculate_mean_std(self):
        for rpc in range(6):
            for phi in range(64):
                for eta in range(32):
                    if len(self.tot_hits[rpc][phi][eta]) > 3:
                        #try to fit a gaussian
                        hist, bins = np.histogram(self.tot_hits[rpc][phi][eta], bins=np.arange(-5, 20, 0.8), range=(-5, 20))
                        bins = bins[:-1]
                        try:
                            popt, pcov = curve_fit(gaus, bins, hist, p0=[1, 13, 5])
                            self.tot_mean[rpc][phi][eta] = popt[1]
                            self.tot_std[rpc][phi][eta] = popt[2]
                        except RuntimeError: #if not, just mean and std                            
                            self.tot_mean[rpc][phi][eta] = np.mean(self.tot_hits[rpc][phi][eta])
                            self.tot_std[rpc][phi][eta] = np.std(self.tot_hits[rpc][phi][eta])

    def caclulate_offset_speed(self):
        #calculate the offset and speed of the TOT
        speed_data = []
        for rpc in range(6):
            for eta in range(32):
                eta_stripe = self.tot_mean[rpc, :, eta]
                x = []
                y = []
                for phi in range(64):
                    if eta_stripe[phi] != 13:
                        y.append(eta_stripe[phi])
                        x.append(phi)
                #fit a line
                if len(y) > 10:
                    popt, pcor = curve_fit(line, x, y, p0=[1, 13])
                    speed_data.append(popt[0])
                #i do not need intercept
        self.speed = np.mean(speed_data)
        for rpc in range(6):
            for eta in range(32):
                eta_stripe = self.tot_mean[rpc, :, eta]
                x = []
                y = []
                for phi in range(64):
                    if eta_stripe[phi] != 13:
                        y.append(eta_stripe[phi])
                        x.append(phi)
                if len(y) > 10:
                    popt, pcor = curve_fit(lambda x, c: line(x, self.speed, c), x, y, p0=[13])
                    self.offsets[rpc, eta] = popt[0] + eta*self.speed
                else:
                    print(rpc, eta)
        hist, bins = np.histogram(speed_data)
        bins = bins[:-1]
        plt.plot(bins, hist, drawstyle='steps-mid')
        plt.show()
        return self.offsets, self.speed

    def plot_tot(self):
        fig_ratio = (20,12)  # Width: 10, Height: 20 to achieve 1:2 ratio
        # Create a PDF to save the plots
        pdf_pages = PdfPages("tot_mean.pdf")
        global_vmin = 0
        global_vmax = 20

        # Loop over the data and create two heatmaps per page
        for i in range(0, 6, 2):
            fig, axes = plt.subplots(2, 1, figsize=fig_ratio)

            # Plot the first heatmap of the pair
            sns.heatmap(np.transpose(self.tot_mean[i]), ax=axes[0], cmap="viridis", cbar=True, vmin=global_vmin, vmax=global_vmax)
            axes[0].set_title(f'RPC {i}')
            axes[0].set_xticks(np.arange(0, 64, 8))
            axes[0].set_yticks(np.arange(0, 32, 8))

            # Plot the second heatmap if it exists
            if i+1 < 6:
                sns.heatmap(np.transpose(self.tot_mean[i+1]), ax=axes[1], cmap="viridis", cbar=True, vmin=global_vmin, vmax=global_vmax)
                axes[1].set_title(f'RPC {i+1}')
                axes[1].set_xticks(np.arange(0, 64, 8))
                axes[1].set_yticks(np.arange(0, 32, 8))

            else:
                axes[1].axis('off')  # Turn off the second axis if there is no second plot

            # Save the current figure to the PDF
            pdf_pages.savefig(fig)
            plt.close(fig)

        # Close the PDF
        pdf_pages.close()
    
        pdf_pages = PdfPages("tot_std.pdf")
        global_vmin = 0
        global_vmax = 10

        # Loop over the data and create two heatmaps per page
        for i in range(0, 6, 2):
            fig, axes = plt.subplots(2, 1, figsize=fig_ratio)

            # Plot the first heatmap of the pair
            sns.heatmap(np.transpose(self.tot_std[i]), ax=axes[0], cmap="viridis", cbar=True, vmin=global_vmin, vmax=global_vmax)
            axes[0].set_title(f'RPC {i}')
            axes[0].set_xticks(np.arange(0, 64, 8))
            axes[0].set_yticks(np.arange(0, 32, 8))

            # Plot the second heatmap if it exists
            if i+1 < 6:
                sns.heatmap(np.transpose(self.tot_std[i+1]), ax=axes[1], cmap="viridis", cbar=True, vmin=global_vmin, vmax=global_vmax)
                axes[1].set_title(f'RPC {i+1}')
                axes[1].set_xticks(np.arange(0, 64, 8))
                axes[1].set_yticks(np.arange(0, 32, 8))

            else:
                axes[1].axis('off')  # Turn off the second axis if there is no second plot

            # Save the current figure to the PDF
            pdf_pages.savefig(fig)
            plt.close(fig)

        # Close the PDF
        pdf_pages.close()

    def plot_point(self, rpc,phi_channel, eta_channel):
        fig, ax = plt.subplots(figsize=(10, 8))
        data = self.tot_hits[rpc][phi_channel][eta_channel]
        
        try:
            hist, bins = np.histogram(data, bins=np.arange(-5, 20, 0.8), range=(-5, 20))
            bins = bins[:-1]
            popt, pcov = curve_fit(gaus, bins, hist, p0=[1, 13, 5])
            fitX = np.linspace(min(bins), max(bins), 400)
            yrange = ax.get_ylim()
            ax.plot(fitX, gaus(fitX, *popt), 'r:', label='Gaussian Fit')
            ax.text(-2, 0.7*yrange[1], f"Fit mean: {round(popt[1], 2)}, $\sigma$: {round(popt[2], 2)}", 
                fontsize=14, verticalalignment='top')

        except RuntimeError:
            print("No Gaussian fit possible")
        # Curve fitting with the new list structure

        # Plot the data as a line plot
        ax.plot(bins, hist, color="teal", lw=3, label='proANUBIS Data', drawstyle='steps-mid')

        ax.set_xlabel('$\eta$ Hit Time - $\phi$ Hit Time (ns)')
      
        
        # Set axis limits and display fit parameters
        ax.set_xlim([-5, 20])
      
        plt.legend()
        plt.show()

def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def line(x, m, c):
    return m*x + c