# OSIRIS
The OperationS and monItoRIng Software framework is a collection of tools to assist in monitoring the online data-taking of proANUBIS and performing low-level data analysis. 

This follows the naming scheme of ANUBIS being an egyptian god, with Osiris being the god of the dead, resurrection, and life. Through this software we should be able to tell if proANUBIS is 'live' and taking data, or 'dead' due to an issue so that we can 'resurrect' it and continue taking our data.  


# Setup
Run the setup.sh script in the top-level directory of OSIRIS, this will automatically create a set of environment variables that may be useful. Such as `MON_PLOT_DIR` which gives the path to a plot directory in the monitoring folder, or `ANUBIS_LOGO` which gives the path to the ANUBIS logo image in the images folder. The script also uses setupATLAS and lsetup to load `LCG_105 x86_64-el9-gcc12-opt` which should contain all needed software. 
As such this framework must be run on lxplus. 
