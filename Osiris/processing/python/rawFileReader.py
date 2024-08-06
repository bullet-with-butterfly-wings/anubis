#from rawEventBuilder import eventBuilder
import rawEventBuilder
import struct
import datetime
import os
import importlib
import Analysis_tools as aTools

def readTimeStamp(tsData):
    ts = struct.iter_unpack("I", tsData)
    tStamp = []
    for tsWord in ts:
        tStamp.append(tsWord[0]) 
    return datetime.datetime.fromtimestamp(tStamp[0]+tStamp[1]/1000000000)

class fileReader():
    def __init__(self,filename,activeTDCs=[tdc for tdc in range(5)]):
        self.fname = filename
        importlib.reload(rawEventBuilder)
        self.evtBuilder = rawEventBuilder.eventBuilder(activeTDCs=activeTDCs)
        self.wordSize = 4
        self.data = open(self.fname,'rb')
        self.data.seek(0, os.SEEK_END)
        self.fsizeBytes = self.data.tell()
        self.data.seek(0)
        self.bytesRead = 0
        stDat = self.data.read(2*self.wordSize)
        self.st = readTimeStamp(stDat)
        self.leadWords = []
        self.adjustment = 0
        self.lastWasBad = False
        self.global_alignment = True
        self.tdcstatus = [True for tdc in range(5)]
        self.tdc_monitoring_event_buffer = []
        self.tdc_monitoring_counter = 0
        self.lasttdcStatus = [True for tdc in range(5)]

    def doneReading(self):
        return self.bytesRead==(self.fsizeBytes-2*self.wordSize)
        
    def hasEvents(self):
        return len(self.evtBuilder.events)>0

    def getEvents(self):
        evts = []
        for event in range(len(self.evtBuilder.events)):
            evts.append(self.evtBuilder.events.pop(0))
        return evts
    def skip_to(self):
        print(len(self.evtBuilder.events))

    def get_aligned_events(self, order = [(0,1), (1,2), (2,3), (3,4)], interval = 100, extract_tdc_mets = False):
        evts_chunk = []
        tdc_mets = [0 for tdc in range(5)]
        TDC_error_time = [[] for tdc in range(5)]
        i = 0
        while i < interval:
            if not self.readBlock():
                print("Bad Block Read")
                break
            if(self.hasEvents()):
                print(len(self.evtBuilder.events))
                for event in range(len(self.evtBuilder.events)):
                    evts_chunk.append(self.evtBuilder.events.pop(0))
                    i += 1
                    self.tdc_monitoring_counter += 1
        aligned, realigned = self.doRealign(evts_chunk, order)
        self.check_alignment_status(aligned, realigned) 
        self.update_adjustment_window(realigned)
        #self.global_alignment = True
        self.tdc_monitoring_event_buffer.extend(evts_chunk)
        if self.tdc_monitoring_counter >= 2500:
            TDC_error_time, tdc_mets = self.monitor_tdc_state(recordtimes=True)
            self.tdc_monitoring_event_buffer.clear()
            self.tdc_monitoring_counter = 0
        if extract_tdc_mets:
            return evts_chunk, tdc_mets, TDC_error_time
        elif self.global_alignment == True and not extract_tdc_mets:
            return evts_chunk
        else:
            return None
                                
    def doRealign(self, event_chunk, order, skipChans=[0]):
        aligned = True
        realigned = False
        offsetlist = [p for o in range(1, (4 + self.adjustment)) for p in (o, -o)]
        updates = [0 for _ in range(len(order))]
        for idx, item in enumerate(order):
            i, j = item
            x, y, l, m = aTools.find_tdc_alignment_metric(i, j)
            alignMet = aTools.calcAvgAlign(event_chunk, offSet=0, i=x, j=y, k=l, l=m, tdc1=i, tdc0=j, processedEvents=0, skipChans=skipChans)#ProcessedEvents is required for class object RPCHit. I know its redundency, but hard to fix
            if alignMet > 15 and alignMet < 100:
                aligned = False
                for testOffset in offsetlist:
                    testAlignMet = aTools.calcAvgAlign(event_chunk, offSet=testOffset, i=x, j=y, k=l, l=m, tdc1=i, tdc0=j, processedEvents=0, skipChans=skipChans)
                    if testAlignMet < 15:
                        updates[idx] += (testOffset)
                        realigned = True
                        break
        insertion_list = aTools.ConstructEventInsertionList(updates, order)
        if not all(x == 0 for x in insertion_list):
            self.InsertFakeEvents(insertion_list)

        return aligned, realigned
    
    
    def check_alignment_status(self, aligned, realigned):
        if not aligned and not realigned:
            self.lastWasBad = True
            self.global_alignment = False
        elif not aligned and realigned:
            self.lastWasBad = False
            self.global_alignment = False
        elif aligned and not realigned:
            self.lastWasBad = False
            self.global_alignment = True
        else:
            print(f'alignment error, realigned events already aligned')
    
    
    def update_adjustment_window(self, realigned):
        if self.lastWasBad and not realigned:
            if self.adjustment < 20:
                self.adjustment += 1
        else:
            self.adjustment = 0
            
    def monitor_tdc_state(self, recordtimes=False):
        TDC_error_time = [[] for tdc in range(5)]
        tdc_mets = [[] for tdc in range(5)]
        bad_channels = [[32],[0,96],[64],[31,32],[0,1,2,3]]
        for tdc in range(5):
            poor_time_count = 0
            good_time_count = 0
            for i, event in enumerate(self.tdc_monitoring_event_buffer[-(2500):]):
                words = event.tdcEvents[tdc].words
                times_words = [(word & 0xfffff, word) for word in words if (word >> 24) & 0x7f not in bad_channels[tdc]]
                
                if times_words:
                    min_time, min_word = min(times_words, key=lambda x: x[0])
                    if recordtimes:
                        TDC_error_time[tdc].append([(min_time, min_word), i])
                    if min_time > 300:
                        poor_time_count += 1
                    elif 200 < min_time <= 300:
                        good_time_count += 1

            if good_time_count == 0:
                tdc_mets[tdc].append(-1)
                if self.lasttdcStatus[tdc]:
                    print(f'tdc{tdc} error state through no good time')
                    self.lasttdcStatus[tdc] = False
            else:
                ratio = poor_time_count / good_time_count
                tdc_mets[tdc].append(ratio)

                if ratio > 0.3:
                    self.tdcstatus[tdc] = False
                    if self.lasttdcStatus[tdc]:
                        print(f'tdc{tdc} enters error state through metric')
                        self.lasttdcStatus[tdc] = False
                else:
                    self.tdcstatus[tdc] = True
                    if not self.lasttdcStatus[tdc]:
                        print(f'tdc{tdc} enters nominal state through metric')
                        self.lasttdcStatus[tdc] = True
                        buffer = self.reload_event_builder()
                        print(f'event builder reloaded, proof {buffer}')
        return TDC_error_time, tdc_mets
        

    def InsertFakeEvents(self, insertion_list):
        for tdc, insertion in enumerate(insertion_list):
            for fakeEvent in range(insertion):
                self.evtBuilder.insertFakeEvent(tdc = tdc)
            
    
    
    
    def readBlock(self,p=False):
       tsDat = self.data.read(2*self.wordSize)
       thisTime = readTimeStamp(tsDat)
       headDat = self.data.read(self.wordSize)
       headWord = 0
       try:
           headWord = struct.unpack("I",headDat)
       except:
           return False
       thisTDC = headWord[0]>>24
       nWords = headWord[0]&0xffffff
       if thisTDC>5:
           print("Bad TDC - number is",thisTDC,", header word was:",hex(leadWords[2]))
           return False
       if self.bytesRead+nWords>(self.fsizeBytes-2*self.wordSize):
           print("Going to over-read the file. Corrupted number of bytes? Header is", hex(leadWords[2]))
           return False
       tdcReadData = self.data.read(nWords*self.wordSize)
       self.bytesRead = self.bytesRead+self.wordSize*(nWords+3)
       #Skip the fifth TDC for now
       if thisTDC<5:
           self.evtBuilder.addTDCRead(thisTDC, thisTime, tdcReadData, p)
       return True
   
    def reload_event_builder(self):
        importlib.reload(rawEventBuilder)
        self.evtBuilder = rawEventBuilder.eventBuilder()
        return self.evtBuilder.tdcEventBuffer
