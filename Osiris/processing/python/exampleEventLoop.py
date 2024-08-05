from rawFileReader import fileReader

fReader = fileReader("my_raw_file.raw")

while not fReader.doneReading():
    #Read a single TDC read block, double-check that it was OK
    goodBlock = fReader.readBlock()
    if not goodBlock:
        break;
    #Only make events when they have data from all TDCs
    if fReader.hasEvents():
        #getEvents empties the event buffer to keep the RAM usage down
        for event in fReader.getEvents():
            #Do some analysis
