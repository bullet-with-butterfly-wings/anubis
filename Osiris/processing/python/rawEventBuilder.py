
import struct

class tdcEvent():
    def __init__(self, header, time, words = [], eob=None, qual=0):
        self.header = header
        self.words = words
        self.time = time
        self.EOB = eob
        self.qual = qual

class proAnubEvent():
    def __init__(self, tdcEvents):
        self.tdcEvents = tdcEvents

class eventBuilder():
    def __init__(self, activeTDCs = [0,1,2,3,4]):
        self.events = []
        self.tdcEventBuffer = [[] for tdc in activeTDCs]
        self.tdcFiveBuffer = []
        self.activeTDCs = activeTDCs
        self.missItr = [[0,0] for tdc in activeTDCs]
        self.eventCounts = [-1 for tdc in activeTDCs]
        self.headerMask = 0x400000
        self.EOBMask = 0x200000
        self.HEADER = 0
        self.EOB = 1
        self.DATA = 2
        self.CORRUPTHEADER = 3
        
    def addTDCRead(self, thisTDC, thisTime, tdcReadData, p=False):
        if thisTDC not in self.activeTDCs and thisTDC!=5:
            return
        readData = struct.iter_unpack("I", tdcReadData)
        thisRead = []
        splitEvtBuf = []
        lastWordType = None
        missedEvts = 0
        maybeMissed = []
        if thisTDC<5:
            startCount = self.eventCounts[thisTDC]
        lastHead = 0
        for word in readData:
            thisWord = word[0]
            if p: #and thisTDC!=5:
                print("TDC",thisTDC,hex(thisWord), (thisWord>>24)&0x7f, thisWord&0xfffff)
            wordType = self.getWordType(thisWord)
            if thisTDC==5:
                thisRead.append(thisWord)
            elif wordType==self.HEADER:
                #If we get a header, check the previous word to throw out extra headers
                if lastWordType==self.HEADER or lastWordType==self.CORRUPTHEADER:
                    thisRead[-1].qual = thisRead[-1].qual|0x10
                    lastWordType=self.EOB
                    continue
                else:
                    #Otherwise, create a new event.
                    #If we haven't had any events yet, start the counter with the first header in case they don't begin with zero
                    if self.eventCounts[thisTDC]<0:
                        self.eventCounts[thisTDC]=(thisWord&0x3ff)+1
                    else:
                        self.eventCounts[thisTDC]=self.eventCounts[thisTDC]+1
                    #Store the difference in TDC header number and the expected event number to check for missed events
                    evtDiff = int(thisWord&0x3ff)-int(self.eventCounts[thisTDC]%(0x3ff+1))+1
                    if p:
                        print("Header from TDC",thisTDC,"evtDiff:",evtDiff,"Saw:",thisWord&0x3ff,"Expected:", self.eventCounts[thisTDC]%(0x3ff+1)-1)
                    if evtDiff<0:
                        evtDiff = evtDiff+0x3ff
                    maybeMissed.append(evtDiff) 
                    thisRead.append(tdcEvent(thisWord, thisTime,[]))
            elif wordType==self.CORRUPTHEADER:
                #Can treat corrupted headers the same way, but as we know the event number is bad we don't check for missed events
                if lastWordType==self.HEADER or lastWordType==self.CORRUPTHEADER:
                    continue
                else:
                    self.eventCounts[thisTDC]=self.eventCounts[thisTDC]+1
                    thisRead.append(tdcEvent(thisWord, thisTime, [], qual=0x1))
            #For non-headers, need to check if we've started a word yet to see if we've got an event that's split over two reads.
            #Need to put split events into a temp buffer, as we don't yet know whether there were missed events in between.
            elif len(thisRead)==0:
                splitEvtBuf.append(thisWord)
            elif wordType==self.EOB:
                #End the current TDC event if we get an EOB, flag the event if we get two in a row.
                if lastWordType==self.EOB:
                    thisRead[-1].qual = thisRead[-1].qual+0x2
                else:
                    thisRead[-1].EOB = thisWord
            else:
                if len(thisRead)==0:
                    thisRead.append(tdcEvent(0, thisTime, [], qual=0x1))
                #Have a data word. Should add them into the event unless we've already got an EOB in the last event, suggesting that we missed a header. 
                if thisRead[-1].EOB is not None:
                    #If we did miss a header, make a fake one with a quality flag set to show we made one up
                    thisRead.append(tdcEvent((self.eventCounts[thisTDC]+1)%(0x3ff+1), thisTime, [], qual=0x5))
                    self.eventCounts[thisTDC]=self.eventCounts[thisTDC]+1
                thisRead[-1].words.append(thisWord)
            lastWordType = wordType
        if thisTDC==5:
            self.tdcFiveBuffer.append([thisTime,thisRead])
            return
        #Now that we've read in the new events, want to check if we missed any in between TDC reads by taking the mode of maybeMissed
        if len(maybeMissed)>0:
            missedEvts = max(set(maybeMissed), key=maybeMissed.count)
            
        if missedEvts>0:
            #Add to the miss counter if more than one read was made 
            #and the offset occurred at least twice
            if len(maybeMissed)>1 and maybeMissed.count(missedEvts)>1:
                if p:
                    print("Inserting "+str(missedEvts)+" into TDC",thisTDC)
                for missedEvt in range(missedEvts):
                    if missedEvts<100:
                        thisRead.insert(0, tdcEvent((startCount+missedEvt)%(0x3ff+1), thisTime, [], qual=0x8))
                    self.eventCounts[thisTDC] = self.eventCounts[thisTDC]+1
                    self.missItr[thisTDC] = [0,0]
            else:
                #If there was only one TDC read, store it in an iterator and wait for another event to confirm the offset.
                if self.missItr[thisTDC][1] == missedEvts:
                    if p:
                        print("Inserting",missedEvts,"into TDC",thisTDC)
                    for evt in range(missedEvts):
                        if missedEvts<100:
                            self.tdcEventBuffer[thisTDC].insert(self.missItr[thisTDC][0],tdcEvent(0xffff, thisTime, [], qual=0x8))
                        self.missItr[thisTDC] = [0,0]
                        self.eventCounts[thisTDC] = self.eventCounts[thisTDC]+1
                else:
                    self.missItr[thisTDC] = [len(self.tdcEventBuffer[thisTDC]), missedEvts]
        
        #If there were no missed events, tack any crossover words onto the proper event from the previous TDC read
        if missedEvts==0:
            self.missItr[thisTDC]=[0,0]
            if len(splitEvtBuf)>0:
                for word in splitEvtBuf:
                    wordType = self.getWordType(word)
                    if len(self.tdcEventBuffer[thisTDC])>0:
                        if wordType==self.EOB:
                            self.tdcEventBuffer[thisTDC][-1].EOB=word
                        else:
                            self.tdcEventBuffer[thisTDC][-1].words.append(word)
                       
        #Now store the tdc events from the new read into the full buffer, and check whether we can pull out any complete proANUBIS events
        for event in thisRead:
            if len(event.words)==0:
                event.qual = 0xf
            #print(event.qual)
            self.tdcEventBuffer[thisTDC].append(event)
        self.buildFullEvents()
        thisRead = []
    
    def buildFullEvents(self):
        while self.checkBufferForEvents():
            fullEvent = []
            for tdc in range(5):
                if tdc in self.activeTDCs:
                    fullEvent.append(self.tdcEventBuffer[tdc].pop(0))
                else:
                    fullEvent.append(tdcEvent(0,0,[],qual=0xff))
            self.events.append(proAnubEvent(fullEvent))
        return
    
    def checkBufferForEvents(self):
        haveAnEvent = True
        for tdc in self.activeTDCs:
            if len(self.tdcEventBuffer[tdc])==0:
                haveAnEvent = False
            #If the last event is missing an EOB, wait for an extra read to see if it got split
            elif len(self.tdcEventBuffer[tdc])==1 and self.tdcEventBuffer[tdc][-1].EOB==None:
                haveAnEvent = False
            #Also don't write events if we might need to insert missing ones
            elif self.missItr[tdc][1] > 0:
                haveAnEvent = False
        return haveAnEvent
           
    def insertFakeEvent(self, tdc):
        self.tdcEventBuffer[tdc].insert(0,tdcEvent(0xffff, 0, [], qual=0x8))
 
    def getWordType(self, word):
        if word&self.headerMask:
            if word&self.EOBMask:
                return self.CORRUPTHEADER
            else:
                return self.HEADER
        elif word&self.EOBMask and (word>>24)>127:
            return self.EOB
        else:
            return self.DATA
