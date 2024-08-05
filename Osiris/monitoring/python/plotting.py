#-------------------------------------------------------------------#
# Script to produce a set of plots based on input data.
#   By default load the most recent raw file from eos/proj-anubis
#-------------------------------------------------------------------#

# Imports:
import sys, os
import anubisPlotUtils as anPlot
import anubisDataUtils as anData
from glob import glob
import re

    
#-----------------------------------------------
# The main function for running standalone
if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('-i','--inFiles', default=[], nargs="+", type=str, help="Input Files", required=False)
    parser.add_argument('--maxRate', default=False, action="store_true", help="Create a plot of the maximum rates within the input files", required=False)
    parser.add_argument('--heatMap', default=False, action="store_true", help="Create a plot of the heatmaps for the input files", required=False)
    parser.add_argument('--channel', default=False, action="store_true", help="Print the per channel counts/rates", required=False)
    parser.add_argument('--plotname', default="", help="Assign a name to the produced plot", required=False)
    args=parser.parse_args()

    if len(args.inFiles)==0:
        print("Using the most recent data file...")
        args.inFiles = [anData.getMostRecentData()]

    if args.maxRate:
        anPlot.plotMaxRates(args.inFiles, time=args.duration)

    if args.heatMap:
        for file in args.inFiles:
            duration = anData.getDuration(file)
            if args.plotname=="":
                plotname ="proANUBISMonitorHeatMap" 
            else:
                plotname = args.plotname

            anPlot.heatFromRawfile(file, time=duration, name=plotname)


    if args.channel:
        for file in args.inFiles:
            print(f"---------- {file} ------------")
            anData.printTotalPerChannel(file, doRate = False)
