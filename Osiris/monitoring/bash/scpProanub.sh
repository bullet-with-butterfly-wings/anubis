#!/bin/bash
MostRecentFile=$(ls /eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/*/*/*.raw -Art | tail -n 1)
#MostRecentFile=$(ls /eos/user/m/mireveri/anubis/data/*/*/*.raw -Art | tail -n 1)
MostRecentCksum=$(cksum $MostRecentFile | cut -d ' ' -f 1)
#rsync -avz proanubis@pcatlanubisbb501:~/michael-test/data/ /eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/
rsync -rvz --rsync-path=C:/cygwin64/bin/rsync.exe proanubis@pcatlanubisbb3:/cygdrive/c/Users/proanubis.PCATLANUBISBB3/Documents/ANUBIS/DataTaking/data/ /eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/
#rsync -avz proanubis@pcatlanubisbb501:~/michael-test/data/ /eos/user/m/mireveri/anubis/data/

export PYTHONPATH="${PYTHONPATH}:/eos/user/m/mireveri/SWAN_projects/anubisPlotting/"
python3 ~/makeMonitorHeat.py
rsync -avz /eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/monitor/*.raw /eos/user/m/mireveri/anubis/data/monitor/
rsync -avz /eos/user/m/mireveri/anubis/data/monitor/ /eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/monitor/

NewMostRecentFile=$(ls /eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/*/*/*.raw -Art | tail -n 1)
#NewMostRecentFile=$(ls /eos/user/m/mireveri/anubis/data/*/*/*.raw -Art | tail -n 1)
NewMostRecentCksum=$(cksum $NewMostRecentFile | cut -d ' ' -f 1)
if [ $MostRecentCksum -eq $NewMostRecentCksum ]
then
   sendEmail=1
   while IFS= read -r line
   do
      sendEmail=$line
   done < "/eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/monitor/emailTrigger.txt"
   if [ $sendEmail -eq 1 ]
   then
      echo "The proANUBIS rsync didn't see any changes in the most-recent datafile over the last 15 minutes. Please check that the data-taking is running and that the server hard drive hasn't filled up." | mail -s "*proANUBIS AUTOMAIL* No Recent Data" "mr2025@cam.ac.uk" "aashaq.shah@cern.ch" "paul.nathaniel.swallow@cern.ch"
      echo 0 > /eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/monitor/emailTrigger.txt
   fi
#   echo "No New Data"
fi
