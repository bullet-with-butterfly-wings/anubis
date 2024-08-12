import math


class rpcHit():
    def __init__(self, channel, time, eta, event_num, rpc):
        self.rpc = rpc
        self.time = time
        self.channel = channel
        self.eta = eta
        self.event_num = event_num
    def __str__(self):
        return f"rpcHit(channel={self.channel}, time={self.time}, eta={self.eta}, event_num={self.event_num}, rpc={self.rpc})"


def tdcChanToRPCHit(word, tdc, event_num):
        tdcChannel = (word >> 24) & 0x7f
        tdcHitTime = word & 0xfffff
        eta = False
        rpcChan = -1
        if tdc == 0:
            if tdcChannel < 32:
                rpcChan = tdcChannel
                eta = True
                rpc = 0
            elif tdcChannel < 96:
                rpcChan = tdcChannel - 32
                eta = False
                rpc = 0
            else:
                rpcChan = tdcChannel - 96
                eta = True
                rpc = 1
        elif tdc == 1:
            if tdcChannel < 64:
                rpcChan = tdcChannel
                eta = False
                rpc = 1
            elif tdcChannel < 96:
                rpcChan = tdcChannel - 64
                eta = True
                rpc = 2
            else:
                rpcChan = tdcChannel - 96
                eta = False
                rpc = 2
        elif tdc == 2:
            if tdcChannel < 32:
                rpcChan = tdcChannel + 32
                eta = False
                rpc = 2
            elif tdcChannel < 64:
                rpcChan = tdcChannel - 32
                eta = True
                rpc = 3
            elif tdcChannel < 128:
                rpcChan = tdcChannel - 64
                eta = False
                rpc = 3
        elif tdc == 3:
            if tdcChannel < 32:
                rpcChan = tdcChannel
                eta = True
                rpc = 4
            elif tdcChannel < 96:
                rpcChan = tdcChannel - 32
                eta = False
                rpc = 4
            else:
                rpcChan = tdcChannel - 96
                eta = True
                rpc = 5
        elif tdc == 4:
            rpcChan = tdcChannel
            eta = False
            rpc = 5
        #was there * 0.8 time idk
        return rpc, rpcHit(rpcChan, (25/32)*tdcHitTime, eta, event_num, rpc)

def find_tdc_alignment_metric(tdc0, tdc1):
    if tdc0 > tdc1:
        tdc0, tdc1 = tdc1, tdc0
    i, j, k, l = None, None, None, None
    if tdc0 == 0:
        if tdc1 == 1:
            i, j, k, l = 1, 2, 0, 1
        if tdc1 == 2:
            i, j, k, l = 3, 0, 3, 0
        if tdc1 == 3:
            i, j, k, l = 4, 0, 4, 0
        if tdc1 == 4:
            i, j, k, l = -1, -1, 5, 0
    if tdc0 == 1:
        if tdc1 == 2:
            i, j, k, l = 3, 2, 3, 1
        if tdc1 == 3:
            i, j, k, l = 5, 2, 4, 1
        if tdc1 == 4:
            i, j, k, l = -1, -1, 5, 1
    if tdc0 == 2:
        if tdc1 == 3:
            i, j, k, l = 4, 3, 4, 3
        if tdc1 == 4:
            i, j, k, l = -1, -1, 5, 3
    if tdc0 == 3:
        if tdc1 == 4:
            i, j, k, l = -1, -1, 5, 4
        
    return i, j, k, l
        
def testAlign(rpc1Hits, rpc2Hits, skipChans = []):
    minTimes = [300,300]
    minChans = [-1,-1]
    if len(rpc1Hits)<1 or len(rpc2Hits)<1:
        return -1
    for hit in rpc1Hits: # might not be the same time for both hits
        if hit.time<minTimes[0] and hit.channel not in skipChans:
            minTimes[0]=hit.time #choice that time hit = min time hit
            minChans[0]=hit.channel
    for hit in rpc2Hits:
        if hit.time<minTimes[1] and hit.channel not in skipChans:
            minTimes[1]=hit.time
            minChans[1]=hit.channel
    if -1 in minChans:
        return -1
    return abs(minChans[1]-minChans[0])


def calcAvgAlign(event_chunk, offSet=0, i = 1, j = 2, k = 0, l = 2, tdc1 =0, tdc0 = 1, processedEvents = 0, skipChans = []):
    mets = []
    for idx, event in enumerate(event_chunk):
        etaHits = [[] for rpc in range(6)]
        phiHits = [[] for rpc in range(6)]
        if (idx+abs(offSet))<len(event_chunk):
            if offSet<=0:
                oneIdx = idx+abs(offSet)
                twoIdx = idx
            else:
                oneIdx = idx
                twoIdx = idx+offSet
            for word in event_chunk[oneIdx].tdcEvents[tdc1].words:
                rpc, thisHit = tdcChanToRPCHit(word,tdc1, processedEvents + idx)
                if thisHit.eta:
                    etaHits[rpc].append(thisHit)

                else:
                    phiHits[rpc].append(thisHit)
            for word in event_chunk[twoIdx].tdcEvents[tdc0].words:
                rpc, thisHit = tdcChanToRPCHit(word,tdc0, processedEvents + idx)
                if thisHit.eta:
                    etaHits[rpc].append(thisHit)

                else:
                    phiHits[rpc].append(thisHit)     
            if i != -1:  
                etOff = testAlign(etaHits[i],etaHits[j], skipChans = skipChans)
                phOff = testAlign(phiHits[k],phiHits[l], skipChans = skipChans)
                if etOff>=0 and phOff>=0:
                    mets.append(math.sqrt(etOff*etOff+phOff*phOff)) # strips do not have different dimensions??
            else:
                phOff = testAlign(phiHits[k],phiHits[l], skipChans = skipChans)
                if phOff>=0:
                    mets.append(math.sqrt(phOff*phOff))
    if len(mets)>0:
        return sum(mets)/len(mets)
    else:
        return -1


def ConstructEventInsertionList(updates, order):
    insertion_list = [0 for _ in range(5)]
    for idx, update in enumerate(updates):
        i, j = order[idx]
        insertion_list[j] += (insertion_list[i] - update)
    min_value = min(insertion_list)
    if min_value < 0:
        addition = abs(min_value)
    else:
        addition = 0
    return [x + addition for x in insertion_list]

