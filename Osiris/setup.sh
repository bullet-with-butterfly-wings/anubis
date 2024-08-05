#!/bin/bash

currentDir=$PWD

for subDir in 'monitoring' 'images' ; do
	if [[ ${currentDir} ==  *"${subDir}"* ]]; then
		echo "This setup script must be run in the top level of OSIRIS"
		return 1
	fi
done

if [[ ${currentDir,,} != *"osiris"* ]]; then
        echo "You are currently in ${currentDir}."
	echo "This setup script must be run in the top level of OSIRIS."
        read -p  "Enter y if this is the top level and you wish to continue with the setup." confirmation
        
        if [[ $confirmation == [yY] || $confirmation == [yY][eE][sS] ]]; then 
            echo "Setting up..."
        else
            echo "Cancelling set up..."
            return 1
        fi
fi

# Setup directory folder bash variables 
export OSIRIS_DIR=$PWD
export MON_DIR=${OSIRIS_DIR}/monitoring
export MON_PLOT_DIR=${OSIRIS_DIR}/monitoring/plots
export IMG_DIR=${OSIRIS_DIR}/images
export PROCESS_DIR=${OSIRIS_DIR}/processing

# Export the logo filepath
export ANUBIS_LOGO=${IMG_DIR}/ANUBISLogo.png


# Setup ATLAS LCG environment
setupATLAS
lsetup "views LCG_105 x86_64-el9-gcc12-opt"
