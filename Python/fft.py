
import datetime
import enum
import math
import os
import random
import shutil
import typing
import unittest

from indexSin import SinSource, IndexSin
import histo
import varDbl


class SignalType (enum.StrEnum):
    Linear = 'Linear',
    Sin = 'Sin',
    Cos = 'Cos',
    Aggr = 'Aggr',

class NoiseType (enum.StrEnum):
    Gaussian = 'Gaussian',
    White = 'White',

class TestType (enum.StrEnum):
    Forward = 'Forward',
    Reverse = 'Reverse',
    Roundtrip = 'Roundtrip',



class FFT:
    _bitReversedIndex = {}

    def __init__(self, sinSource:SinSource):
        self.idxSin = IndexSin(sinSource) 


    @staticmethod
    def bitReversedIndices(order:int) -> tuple[int]:
        '''
        Return bit reversed indices
        from NumericalRecipesinC.pdf
        '''
        if sRes := FFT._bitReversedIndex.get(order):
            return sRes
        N = 1 << order
        sRes = [0] * N
        M = N >> 1
        j = 0
        for i in range(N):
            sRes[i] = j
            k = M 
            while (k != 0) and (j >= k):
                j -= k
                k >>= 1
            j += k
        sRes = tuple(sRes)
        FFT._bitReversedIndex[order] = sRes
        return sRes
    
    def transform(self, sInput:tuple[varDbl.VarDbl], forward:bool, 
                  traceSteps=False):
        '''
        {sInput}: a VarDbl array of size (2<<order), with each datum contains (real, image)
        {forward}: true for forware transformation, false for backward transformation.

        When {sSin} and {sCos} are not None, output the result.
        '''        
        order = IndexSin.validateSize(len(sInput) >> 1) 
        size = 1 << order
        sRes = [0] * (size << 1)
        self.ssStep = []
        
        sIndex = FFT.bitReversedIndices(order)
        for i in range(len(sIndex)):
            j = sIndex[i]
            sRes[i << 1], sRes[(i << 1) + 1] = sInput[j << 1], sInput[(j << 1) + 1]
        if traceSteps:
            self.ssStep.append(tuple([varDbl.VarDbl(var) for var in sRes]))

        for i in range(0, len(sIndex), 2):
            j = i << 1
            sRes[j], sRes[j + 2] = sRes[j] + sRes[j + 2], sRes[j] - sRes[j + 2]
            sRes[j + 1], sRes[j + 3] = sRes[j + 1] + sRes[j + 3], sRes[j + 1] - sRes[j + 3]
        if traceSteps:
            self.ssStep.append(tuple([varDbl.VarDbl(var) for var in sRes]))

        for o in range(1, order):
            k = 2 << o
            for j in range(k >> 1):
                cos = self.idxSin.cos(j, o) 
                sin = self.idxSin.sin(j if forward else -j, o)
                for i in range(0, size, k):
                    i0 = (i + j) << 1
                    i1 = i0 + k
                    rd = sRes[i1] * cos - sRes[i1 + 1] * sin
                    id = sRes[i1] * sin + sRes[i1 + 1] * cos
                    sRes[i0], sRes[i1] = sRes[i0] + rd, sRes[i0] - rd
                    sRes[i0 + 1], sRes[i1 + 1] = sRes[i0 + 1] + id, sRes[i0 + 1] - id
            if traceSteps:
                self.ssStep.append(tuple([varDbl.VarDbl(var) for var in sRes]))

        if not forward:
            sz = 1/len(sIndex)
            for i in range(len(sRes)):
                sRes[i] *= sz
        if traceSteps:
            self.ssStep.append(tuple([varDbl.VarDbl(var) for var in sRes]))
        return sRes
    

class Measure:
    def __init__(self, divids=5, devs=3) -> None:
        self.sUncStat = {
            TestType.Forward: histo.Stat(), 
            TestType.Reverse: histo.Stat(), 
            TestType.Roundtrip: histo.Stat()}
        self.sValStat = {
            TestType.Forward: histo.Stat(), 
            TestType.Reverse: histo.Stat(), 
            TestType.Roundtrip: histo.Stat()}
        self.sValStat = {
            TestType.Forward: histo.Stat(), 
            TestType.Reverse: histo.Stat(), 
            TestType.Roundtrip: histo.Stat()}
        self.sHisto = {
            TestType.Forward: histo.Histo(divids, devs), 
            TestType.Reverse: histo.Histo(divids, devs), 
            TestType.Roundtrip: histo.Histo(divids, devs)}


class FFT_Signal (FFT): 
    '''
    A complex signal {sWave} and its spectrum {sFreq} of size (2<<order).
    The signal can be specified by either {signalType}, {order} and {freq}, or {sWave} and {sFreq}.
    The indexed sin and cos values are from either {sinSource} or {sCosSin}.
    '''

    def __init__(self, sinSource:SinSource, signalType:SignalType, order:int, freq:int,
                 sCosSin:tuple[varDbl.VarDbl]=None, sWave:tuple[varDbl.VarDbl]=None, sFreq:tuple[varDbl.VarDbl]=None):
        self.sinSource = sinSource
        self.idxSin = IndexSin(sinSource, sCosSin)
        self.signalType = signalType
        self.order = order
        self.size = 1 << order
        self.freq = freq

        if sWave and len(sWave) == (self.size << 1) and sFreq and len(sFreq) == (self.size << 1):
            self.sWave = sWave
            self.sFreq = sFreq
            return
        self.sWave = []
        self.sFreq = []
        half = self.size >> 1
        match self.signalType:
            case SignalType.Sin:
                if not (1 <= freq <= half):
                    raise RuntimeError(f'Invalid freq={freq} for order={order}')
                for i in range(self.size):
                    self.sWave.append(self.idxSin.sin(freq * i * 2, order))
                    self.sWave.append(varDbl.VarDbl(0))
                    if i == freq:
                        self.sFreq.append(varDbl.VarDbl(0))
                        self.sFreq.append(varDbl.VarDbl(half))
                    elif i == (self.size - freq):
                        self.sFreq.append(varDbl.VarDbl(0))
                        self.sFreq.append(varDbl.VarDbl(-half))
                    else:
                        self.sFreq.append(varDbl.VarDbl(0))
                        self.sFreq.append(varDbl.VarDbl(0))
            case SignalType.Cos:
                if not (1 <= freq <= half):
                    raise RuntimeError(f'Invalid freq={freq} for order={order}')
                for i in range(self.size):
                    self.sWave.append(self.idxSin.cos(freq * i * 2, order))
                    self.sWave.append(varDbl.VarDbl(0))
                    if i == freq:
                        self.sFreq.append(varDbl.VarDbl(half))
                        self.sFreq.append(varDbl.VarDbl(0))
                    elif i == (self.size - freq):
                        self.sFreq.append(varDbl.VarDbl(half))
                        self.sFreq.append(varDbl.VarDbl(0))
                    else:
                        self.sFreq.append(varDbl.VarDbl(0))
                        self.sFreq.append(varDbl.VarDbl(0))
            case SignalType.Linear:
                if freq:
                    raise RuntimeError(f'Invalid freq={freq} for order={order}')
                for i in range(self.size):
                    self.sWave.append(varDbl.VarDbl(i))
                    self.sWave.append(varDbl.VarDbl(0))
                    if i:
                        self.sFreq.append(varDbl.VarDbl(-half))
                        self.sFreq.append(varDbl.VarDbl(-half * self.idxSin.cos(i, order) / self.idxSin.sin(i, order)))
                    else:
                        self.sFreq.append(varDbl.VarDbl(self.size * (self.size - 1) / 2))
                        self.sFreq.append(varDbl.VarDbl(0))
            case _:
                raise RuntimeError(f"Unknown signal={signalType} for sinSource={sinSource} order={order} freq={freq}")
        if len(self.sWave) != (self.size << 1):
            raise RuntimeError(f"Invalid wave length {len(self.sWave)} vs {self.size << 1} for sinSource={sinSource} order={order} signalType={signalType} freq={freq}")
        if len(self.sFreq) != (self.size << 1):
            raise RuntimeError(f"Invalid spec length {len(self.sFreq)} vs {self.size << 1} for sinSource={sinSource} order={order} signalType={signalType} freq={freq}")



class FFT_Order (FFT_Signal):
    '''
    Perform FFT for a FFT_Signal, with noise of {noiseType} and {noise}.

    FFT_Order.dump() dumps the result to file.  
        Due to FFT_Order.read(), FFT_Order.dump() can continue to append the result to existing file.
    FFT_Order.read() reads the result from file.
    FFT_Order.sort() sorts the result in the file in a standard order.
    FFT_Order.compare() compares two result files 
    '''
    MAX_FREQ = 8
    DIVIDS = 5
    DEVS = 3
    MIN_COUNT = 64
    NORMALIZED_ERROR_OUTLIER = 1e15
    ssssAggr = {}

    def __init__(self, signal:FFT_Signal, noiseType:NoiseType, noise:float,
                 sCosSin:tuple[varDbl.VarDbl]=None,
                 sWave:tuple[varDbl.VarDbl]=None, sFreq:tuple[varDbl.VarDbl]=None,
                 sFrwd:tuple[varDbl.VarDbl]=None, sBack:tuple[varDbl.VarDbl]=None,
                 traceSteps=False, minCount=MIN_COUNT):
        super().__init__(SinSource.Limit if sCosSin else signal.sinSource, 
                         signal.signalType, signal.order, signal.freq,
                         sCosSin=sCosSin, sWave=sWave, sFreq=sFreq)

        self.noiseType = noiseType
        self.noise = abs(noise)
        match noiseType:
            case NoiseType.Gaussian:
                pass
            case NoiseType.White:
                pass
            case _:
                raise ValueError(f'Invalid noiseType={noiseType}')

        if sFrwd:
            if len(sFrwd) != (self.size << 1):
                raise RuntimeError(f'Invalid forward input size {len(sFrwd)} for order {self.order}: {sFrwd}')
            else:
                self.sFrwd = sFrwd
        elif self.noise == 0:
            self.sFrwd = [varDbl.VarDbl(self.sWave[i]) for i in range(self.size << 1)]
        else:
            self.sFrwd = [varDbl.VarDbl(self.sWave[i]) + varDbl.VarDbl(self.getNoise(), self.noise) for i in range(self.size << 1)]

        if sBack:
            if len(sBack) != (self.size << 1):
                raise RuntimeError(f'Invalid backward input size {len(sBack)} for order {self.order}: {sBack}')
            else:
                self.sBack = sBack
        elif self.noise == 0:
            self.sBack = [varDbl.VarDbl(self.sFreq[i]) for i in range(self.size << 1)]
        else:
            self.sBack = [varDbl.VarDbl(self.sFreq[i]) + varDbl.VarDbl(self.getNoise(), self.noise) for i in range(self.size << 1)]

        self.measure = Measure(FFT_Order.DIVIDS, FFT_Order.DEVS)
        self.calc(traceSteps)
        while (not traceSteps) and (0 < noise) and (self.measure.sUncStat[TestType.Roundtrip].count() < minCount):
            self.calc(traceSteps)

        if traceSteps:
            self.ssSpecStep.append(self.sFreq)
            self.ssRoundStep.append(self.sFrwd)
            self.ssRevStep.append(self.sWave)
            self.ssSpecStep.append([actual - expected for actual, expected in zip(self.sSpec, self.sFreq)])
            self.ssRoundStep.append([actual - expected for actual, expected in zip(self.sRound, self.sFrwd)])
            self.ssRevStep.append([actual - expected for actual, expected in zip(self.sRev, self.sWave)])

    def calc(self, traceSteps:bool):
        self.sSpec = self.transform(self.sFrwd, True, traceSteps=traceSteps)
        self.ssSpecStep = self.ssStep
        self.sRound = self.transform(self.sSpec, False, traceSteps=traceSteps)
        self.ssRoundStep = self.ssStep
        self.sRev = self.transform(self.sBack, False, traceSteps=traceSteps)
        self.ssRevStep = self.ssStep

        if self.signalType == SignalType.Linear:
            self.aggr = None
        else:
            self.aggr = FFT_Order.ssssAggr.setdefault(self.order, {}).setdefault(self.sinSource, {})\
                            .setdefault(self.noiseType, {}).setdefault(self.noise, Measure(FFT_Order.DIVIDS, FFT_Order.DEVS))     
        for i in range(self.size << 1):
            self.accum(TestType.Forward, i, self.sSpec[i], self.sFreq[i])
            self.accum(TestType.Roundtrip, i, self.sRound[i], self.sFrwd[i])
            self.accum(TestType.Reverse, i, self.sRev[i], self.sWave[i])

    def getNoise(self) -> float:
        match self.noiseType:
            case NoiseType.Gaussian:
                return random.gauss(0, self.noise)
            case NoiseType.White:
                return random.uniform(-self.noise*math.sqrt(3), self.noise*math.sqrt(3))
            case _:
                raise ValueError(f'Invalid noiseType={self.noiseType}')

    def accum(self, test:TestType, index:int, actual:varDbl.VarDbl, expected:varDbl.VarDbl):
        err = actual - expected
        self.measure.sUncStat[test].accum(actual.uncertainty(), index)
        self.measure.sValStat[test].accum(err.value(), index)
        if self.aggr:
            self.aggr.sUncStat[test].accum(actual.uncertainty(), index)
            self.aggr.sValStat[test].accum(err.value(), index)
        if err.uncertainty() > 0:
            try:
                norm = err.value() / err.uncertainty()
            except OverflowError as ex:
                raise ex 
            if abs(norm) < FFT_Order.NORMALIZED_ERROR_OUTLIER:
                self.measure.sHisto[test].accum(norm, index)
                if self.aggr:
                    self.aggr.sHisto[test].accum(norm, index)
            else:
                print(f'For signal={self.signalType} freq={self.freq} noiseType={self.noiseType} noise={self.noise} test={test} index={index}, normaliized error outlier {norm} between {actual} and {expected}')

    @staticmethod
    def title(divids, devs) -> str:
        return "SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest"\
               "\tUncertainty Count\tUncertainty Mean\tUncertainty Deviation\tUncertainty Minimum\tUncertainty Minimum At\tUncertainty Maximum\tUncertainty Maximum At"\
               "\tValue Count\tValue Mean\tValue Deviation\tValue Minimum\tValue Minimum At\tValue Maximum\tValue Maximum At"\
               "\tError Count\tError Mean\tError Deviation\tError Minimum\tError Minimum At\tError Maximum\tError Maximum At"\
               "\tLower Count\tUpper Count\t"\
               + '\t'.join([f'{i/divids}' for i in range(-devs * divids, devs * divids + 1)])\
               + '\n' 

    @staticmethod
    def uncertainty_offset():
        return 7
    
    @staticmethod
    def error_offset():
        return 21
    
    @staticmethod
    def histo_offset():
        return 30

    def dumpMeasure(self, fw, signalType:SignalType, measure:Measure):
        for test in TestType:
            try:
                fw.write(f'{self.sinSource}\t{self.noiseType}\t{self.noise}\t{signalType}\t{self.order}\t{self.freq}\t{test}')
                for stat in (measure.sUncStat[test], measure.sValStat[test], measure.sHisto[test].stat()):
                    fw.write(f'\t{stat.count()}\t{stat.mean()}\t{stat.dev()}\t{stat.min()}\t{stat.minAt()}\t{stat.max()}\t{stat.maxAt()}')
                fw.write(f'\t{measure.sHisto[test].less()}\t{measure.sHisto[test].more()}')
                count = sum([c for c in measure.sHisto[test].histogram()])
                for c in measure.sHisto[test].histogram():
                    if count:
                        fw.write(f'\t{c/count}')
                    else:
                        fw.write(f'\t{c}')
                fw.write('\n')
            except BaseException as ex:
                print(f'{self.sinSource}\t{self.noiseType}\t{self.noise}\t{self.signalType}\t{self.order}\t{self.freq}\t{test}: {ex}')
                raise ex

    @staticmethod
    def dumpPath(sOrder:tuple[int]):
        if os.getcwd().endswith('VarianceArithmetic'):
            return f'./Python/Output/FFT_{min(sOrder)}_{max(sOrder) + 1}.txt'
        elif os.getcwd().endswith('Python'):
            return f'./Output/FFT_{min(sOrder)}_{max(sOrder) + 1}.txt'
        else:
            raise ValueError(f'Invalid cwd {os.getcwd()}')
    
    @staticmethod
    def readLine(line:str, ln:int, dumpPath:str) -> tuple[SinSource, NoiseType, float, SignalType, int, int, TestType]:
        sWord = line.split('\t')[:FFT_Order.uncertainty_offset()]
        sinSource = SinSource(sWord[0])
        noiseType = NoiseType(sWord[1])
        noise = float(sWord[2])
        if noise < 0:
            raise RuntimeError(f'Invalid noise={noise} in #{ln} line of {dumpPath}: {line}')
        signalType = SignalType(sWord[3])
        order = int(sWord[4])
        if not (IndexSin.MIN_ORDER <= order <= IndexSin.MAX_ORDER):
            raise RuntimeError(f'Invalid order={order} in #{ln} line of {dumpPath}: {line}')
        freq = int(sWord[5])
        if signalType in (SignalType.Linear, SignalType.Aggr): 
            if freq:
                raise RuntimeError(f'Invalid freq={freq} in #{ln} line of {dumpPath}: {line}')
        elif signalType in (SignalType.Sin, SignalType.Cos): 
            if not (0 <= freq < min(FFT_Order.MAX_FREQ, 1 << (order - 1))):
                raise RuntimeError(f'Invalid freq={freq} in #{ln} line of {dumpPath}: {line}')
        test = TestType(sWord[6])
        return sinSource,noiseType,noise,signalType,order,freq,test

    @staticmethod
    def read(dumpPath:str, filterSignalType=SignalType.Aggr, filterFreq=0, 
             checkEnding=True, minOrder=3, minNoise=1e-15) \
        -> dict[int, dict[SinSource, dict[NoiseType, dict[float, dict[TestType, tuple[float, float, int]]]]]]:
        '''
        Read the error deviation, the uncertainty mean, and the uncertainty count 
            of a particular {filterSignalType} and {filterFreq}
            from file {dumpPath} as a dictionary.
        If {checkEnding} is true, the last line must be of {filterSignalType}.
            This to make sure that the ending is SignalType.Aggr, as a complete calculation.
        '''
        if not os.path.isfile(dumpPath):
            return
        sssssAggr = {}
        with open(dumpPath) as f:
            hdr = next(f).split('\t')
            header = FFT_Order.title(FFT_Order.DIVIDS, FFT_Order.DEVS).split('\t')
            if hdr != header:
                if hdr[:FFT_Order.histo_offset()] != header[:FFT_Order.histo_offset()]:
                    print(f'Warning: Invalid title of {dumpPath}: {hdr[:FFT_Order.histo_offset()]} vs {header[:FFT_Order.histo_offset()]}')  
                if [float(hist) for hist in hdr[FFT_Order.histo_offset():]] != [float(hist) for hist in header[FFT_Order.histo_offset():]]:
                    print(f'Warning: Invalid histogram of {dumpPath}: {hdr[FFT_Order.histo_offset():]} vs {header[FFT_Order.histo_offset()]:}')  
            n = 0
            for ln, line in enumerate(f):
                try:
                    sinSource,noiseType,noise,signalType,order,freq,test = FFT_Order.readLine(line, ln, dumpPath)
                    if order <= 1:
                        continue
                except BaseException as ex:
                   raise ex
                if signalType == filterSignalType and freq == filterFreq:
                    sWord = line.split('\t')
                    cntUnc = int(sWord[FFT_Order.uncertainty_offset()])
                    if cntUnc < (2 << order):
                        raise AssertionError(f'Invalid Uncertainty Count {cntUnc} < {(2 << order)} for order {order} at #{ln}: {line}')
                    cntErr = int(sWord[FFT_Order.error_offset()])
                    if (minNoise < noise) and (minOrder < order):
                        if cntErr != cntUnc:
                            raise AssertionError(f'Invalid Error Count {cntErr} vs Uncertainty Count {cntUnc} for order {order} at #{ln} with offset {FFT_Order.error_offset()}: {line}')
                    else:
                        if cntErr > cntUnc:
                            raise AssertionError(f'Invalid Error Count {cntErr} vs Uncertainty Count {cntUnc} for order {order} at #{ln}: {line}')
                    uncMean = float(sWord[FFT_Order.uncertainty_offset() + 1])
                    errDev = float(sWord[FFT_Order.error_offset() + 2])
                    sssssAggr.setdefault(order, {}).setdefault(sinSource, {})\
                             .setdefault(noiseType, {}).setdefault(noise, {})\
                             .setdefault(test, (errDev, uncMean, cntErr))
                    n += 1     
            if checkEnding and signalType != SignalType.Aggr:
                raise RuntimeError(f'Invalid end line of {dumpPath}: {line}')
        print(f'Finish reading {ln} lines for {n} {filterSignalType} records from {dumpPath}')
        return sssssAggr

    @staticmethod
    def dump(sOrder=range(2, IndexSin.MAX_ORDER + 1), 
             sSinSource=(SinSource.Prec, SinSource.Quart, SinSource.Full, SinSource.Lib, SinSource.Fixed),
             sFreq = range(1, MAX_FREQ),
             sNoise=[0] + [math.pow(10, n) for n in range(-17, 1)],
             sNoiseType=(NoiseType.Gaussian, NoiseType.White)):
        '''
        Dump the FFT for order between [ndexSin.MIN_ORDER, IndexSin.MAX_ORDER), freq between [1, MAX_FREQ)
            noise of type (NoiseType.Gaussian, NoiseType.White) and level between [0, 1e-17, ... 1]
        For statistical significancy, the minimal count is FFT_Order.MIN_COUNT
        '''
        dumpPath = FFT_Order.dumpPath(sOrder)
        sssssAggr = FFT_Order.read(dumpPath)

        with open(dumpPath, 'a' if sssssAggr else 'w') as fw:
            if not sssssAggr:
                fw.write(FFT_Order.title(FFT_Order.DIVIDS, FFT_Order.DEVS))
            for noiseType in sNoiseType:
                for order in sOrder:
                    half = 1 << (order - 1)
                    for sinSource in sSinSource:
                        sSignal = None
                        for noise in sNoise:
                            if sssssAggr and (ssssAggr := sssssAggr.get(order)) and (sssAggr := ssssAggr.get(sinSource)) \
                                    and (ssAggr := sssAggr.get(noiseType)) and (sAggr := ssAggr.get(noise)) and (len(sAggr) == 3):
                                continue
                            print(f'{datetime.datetime.now()}: Start calulation order={order}, sinSource={sinSource}, noiseType={noiseType}, noise={noise}')
                            if not sSignal:
                                sSignal = [FFT_Signal(sinSource, SignalType.Sin, order, freq) for freq in sFreq if freq < half] +\
                                          [FFT_Signal(sinSource, SignalType.Cos, order, freq) for freq in sFreq if freq < half] +\
                                          [FFT_Signal(sinSource, SignalType.Linear, order, 0)]
                                print(f'{datetime.datetime.now()}: Finish create signal for order={order}, sinSource={sinSource}')
                            for signal in sSignal:
                                calc = FFT_Order(signal, noiseType, noise)
                                calc.dumpMeasure(fw, calc.signalType, calc.measure)
                            # the last one is linear   
                            calc.dumpMeasure(fw, SignalType.Aggr, FFT_Order.ssssAggr[order][sinSource][noiseType][noise])
                            fw.flush()

    @staticmethod
    def sort(dumpPath:str=None):
        '''
        Obselete: Order the FFT dump to be comparable with those from C++ and Java
        '''
        if not dumpPath:
            dumpPath = FFT_Order.dumpPath((IndexSin.MIN_ORDER, IndexSin.MAX_ORDER - 1))
        if not os.path.isfile(dumpPath):
           return
        shutil.copy(dumpPath, dumpPath + '.bak')
        
        sssssssSort = {}
        with open(dumpPath) as f:
            title = next(f)
            if FFT_Order.title(FFT_Order.DIVIDS, FFT_Order.DEVS) != title:
                raise RuntimeError(f'Invalid title line in {dumpPath}: {title}')
            for ln, line in enumerate(f):
                try:
                    sinSource,noiseType,noise,signalType,order,freq,test = FFT_Order.readLine(line, ln, dumpPath)
                except BaseException as ex:
                   raise ex
                sssssssSort.setdefault(order, {}).setdefault(sinSource, {})\
                           .setdefault(noiseType, {}).setdefault(noise, {})\
                           .setdefault(signalType, {}).setdefault(freq, {})\
                           .setdefault(test, line)
            if signalType != SignalType.Aggr:
                raise RuntimeError(f'Invalid end line of {dumpPath}: {line}')
                
        with open(dumpPath, 'w') as f:
            f.write(FFT_Order.title(FFT_Order.DIVIDS, FFT_Order.DEVS))
            for order in sorted(sssssssSort):
                ssssssSort = sssssssSort[order]
                for sinSource in sorted(ssssssSort):
                    sssssSort = ssssssSort[sinSource]
                    for noiseType in sorted(sssssSort):
                        ssssSort = sssssSort[noiseType]
                        for noise in sorted(ssssSort):
                            sssSort = ssssSort[noise]
                            for signalType in (SignalType.Sin, SignalType.Cos, SignalType.Linear, SignalType.Aggr):
                                ssSort = sssSort[signalType]
                                for freq in sorted(ssSort):
                                    sSort = ssSort[freq]
                                    for test in (TestType.Forward, TestType.Roundtrip, TestType.Reverse):
                                        f.write(sSort[test])

    @staticmethod
    def compare(testCase:unittest.TestCase, thisPath:str, thatPath:str, 
                filterSignalType=SignalType.Aggr, filterFreq=0,
                precErr=5e-1, precUnc=1e-6, matchCount=False, 
                minNoise=1e-15, minCount=MIN_COUNT, minOrder=2):
        '''
        Compare two FFT dump files {thisPath} and {thatPath} for uncertainty mean and error deviation.
        When {matchCount} is true, also compare the uncertainty count.
        Only compare the error deviation for noise >= {minNoise}.
        Only compare the error deviation for uncertainty count >= {minCount}
        ''' 

        sKey = []
        sssssThis = FFT_Order.read(thisPath, filterSignalType, filterFreq, checkEnding=False)
        sssssThat = FFT_Order.read(thatPath, filterSignalType, filterFreq, checkEnding=False)
        sThis = set(sssssThis.keys())
        sThat = set(sssssThat.keys())
        if sThis != sThat:
            print(f'Diff orders:'
                  f' {sThis - sThat} only in {thisPath}, {sThat - sThis} only in {thatPath}')
        for order in sorted(sThis & sThat):
            if order <= minOrder:
                continue
            ssssThis = sssssThis[order]
            ssssThat = sssssThat[order]
            sThis = set(ssssThis.keys())
            sThat = set(ssssThat.keys())
            if sThis != sThat:
                sKey.append(f'Diff sinSource for order={order}: '
                            f'{sThis - sThat} only in {thisPath}, {sThat - sThis} only in {thatPath}')
            for sinSource in sorted(sThis & sThat):
                sssThis = ssssThis[sinSource]
                sssThat = ssssThat[sinSource]
                sThis = set(sssThis.keys())
                sThat = set(sssThat.keys())
                if sThis != sThat: 
                    sKey.append(f'Diff noiseType for order={order} sinSource={sinSource}: '
                                f'{sThis - sThat} only in {thisPath}, {sThat - sThis} only in {thatPath}')
                for noiseType in sorted(sThis & sThat):
                    ssThis = sssThis[noiseType]
                    ssThat = sssThat[noiseType]
                    sThis = set(ssThis.keys())
                    sThat = set(ssThat.keys())
                    if sThis != sThat: 
                        sKey.append(f'Diff noise for order={order} sinSource={sinSource} noiseType={noiseType}: '
                                    f'{sThis - sThat} only in {thisPath}, {sThat - sThis} only in {thatPath}')
                    for noise in sorted(sThis & sThat):
                        sThis = set(ssThis[noise].keys())
                        sThat = set(ssThat[noise].keys())
                        if sThis != sThat:
                            sKey.append(f'Diff test for order={order} sinSource={sinSource} noiseType={noiseType} noise={noise}: '
                                        f'{sThis - sThat} only in {thisPath}, {sThat - sThis} only in {thatPath}')
                        sCommon = sThis & sThat
                        sThis = ssThis[noise]
                        sThat = ssThat[noise]
                        for test in sCommon:
                            this = sThis[test]
                            that = sThat[test]
                            if matchCount:
                                try:
                                    testCase.assertEqual(this[2], that[2])
                                except AssertionError as ex:
                                    testCase.fail(f'Diff count {this[2]} vs {that[2]} for order={order} sinSource={sinSource} noiseType={noiseType} noise={noise} test={test}: {ex}')                           
                            try:
                                testCase.assertAlmostEqual(this[1] / that[1], 1, delta = precUnc)
                            except AssertionError as ex:
                                testCase.fail(f'Diff uncertainty {this[1]} vs {that[1]} for order={order} sinSource={sinSource} noiseType={noiseType} noise={noise} test={test}: {ex}')
                            if (test != TestType.Roundtrip) and (minNoise <= noise) and (minCount <= this[2]) and (minCount <=that[2]):
                                try:
                                    testCase.assertAlmostEqual(this[0] / that[0], 1, delta = precErr)
                                except AssertionError as ex:
                                    testCase.fail(f'Diff error {this[0]} vs {that[0]} for order={order} sinSource={sinSource} noiseType={noiseType} noise={noise} test={test}: {ex}')

class FFT_Signal_Param: 
    def __init__(self, sinSource:SinSource, signalType:SignalType, order:int, freq:int):
        self.sinSource = sinSource
        self.signalType = signalType
        self.order = order
        self.size = 1 << order
        self.freq = freq


class FFT_Order_Param:
    def __init__(self, signal:FFT_Signal_Param, noiseType:NoiseType, noise:float):
        self.signal = signal
        self.noiseType = noiseType
        self.noise = noise


class FFT_Step (FFT_Order):
    '''
    Test the FFT step by step.

    FFT_Step.dump() dumps the result to file.
    FFT_Step.test() tests the result from the dump file by recalculating.
    '''
    EXTRA_STEPS = ('Output', 'Expected', 'Error')

    def __init__(self, signal:FFT_Signal, noiseType:NoiseType, noise:float,
                 sCosSin:tuple[varDbl.VarDbl]=None,
                 sWave:tuple[varDbl.VarDbl]=None, sFreq:tuple[varDbl.VarDbl]=None,
                 sFrwd:tuple[varDbl.VarDbl]=None, sBack:tuple[varDbl.VarDbl]=None):
        super().__init__(signal, noiseType, noise,
                         sCosSin, sWave, sFreq, sFrwd, sBack, True)

    @staticmethod
    def header(order:int):
        return 'SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest\tStep\tImag\tValue\t' +\
            '\t'.join([f'{i}' for i in range(1 << order)])
    
    @staticmethod
    def dataOffset():
        return 10

    @staticmethod
    def dumpPath(order, sinSource):
        if os.getcwd().endswith('VarianceArithmetic'):
            return f'./Python/Output/FFT_Step_{order}_{sinSource}.txt'
        elif os.getcwd().endswith('Python'):
            return f'./Output/FFT_Step_{order}_{sinSource}.txt'
        else:
            raise ValueError(f'Invalid cwd {os.getcwd()}')

    @staticmethod
    def dump(order:int, sinSource: SinSource,
             sNoiseType = (NoiseType.Gaussian, NoiseType.White),
             sNoise:tuple[float]=(0, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1), 
             sFreq=range(1, FFT_Order.MAX_FREQ),
             dumpStepPath:str=None,
             dumpOrderPath:str=None):
        if not IndexSin.MIN_ORDER <= order < IndexSin.MAX_ORDER:
            raise RuntimeError(f'Invalid order {order}')
        size = 1 << order
        fft = FFT(sinSource)
        if dumpStepPath:
            dumpPath = dumpStepPath
        else:
            dumpPath = FFT_Step.dumpPath(order, sinSource) 
        if os.path.isfile(dumpPath):
            os.remove(dumpPath)

        def writeData(fw, sData:tuple[varDbl.VarDbl], 
                      noiseType:NoiseType, noise:float, 
                      signal:SignalType, freq:int, test:TestType, step):
            for imag in (0, 1):
                for val in (1, 0):
                    fw.write(f'\n{sinSource}\t{noiseType}\t{noise}\t{signal}\t{order}\t{freq}\t{test}\t{step}\t{imag}\t{val}')
                    for i in range(imag, size << 1, 2):
                        fw.write(f'\t{sData[i].value() if val else sData[i].uncertainty()}')
                fw.flush()

        fo = open(dumpOrderPath, 'w') if dumpOrderPath else None
        if fo:
            fo.write(FFT_Order.title(FFT_Order.DIVIDS, FFT_Order.DEVS))
        print(f'{datetime.datetime.now()}: Start dump to {dumpPath}')
        with open(dumpPath, 'w') as fw:
            fw.write(FFT_Step.header(order))
            sCosSin = []
            for i in range(size):
                sCosSin.append(fft.idxSin.cos(i, order))
                sCosSin.append(fft.idxSin.sin(i, order))
            writeData(fw, sCosSin, '', 0, '', '', 'CosSin', '')
            
            print(f'{datetime.datetime.now()}: Start calulation order={order}, sinSource={sinSource}')
            half = size >> 1
            sSignal = [FFT_Signal(sinSource, SignalType.Sin, order, freq) for freq in sFreq if freq < half] +\
                      [FFT_Signal(sinSource, SignalType.Cos, order, freq) for freq in sFreq if freq < half] +\
                      [FFT_Signal(sinSource, SignalType.Linear, order, 0)]

            for noiseType in sNoiseType:
                for noise in sNoise:
                    for signal in sSignal:
                        calc = FFT_Step(signal, noiseType, noise)
                        for test, sInput, ssStep in ((TestType.Forward, calc.sFrwd, calc.ssSpecStep),
                                                     (TestType.Roundtrip, calc.sSpec, calc.ssRoundStep),
                                                     (TestType.Reverse, calc.sBack, calc.ssRevStep)):
                            writeData(fw, sInput, calc.noiseType, calc.noise,
                                      calc.signalType, calc.freq, test, 'Input')
                            for step, sStep in enumerate(ssStep):
                                writeData(fw, sStep, calc.noiseType, calc.noise,
                                          calc.signalType, calc.freq, test, step)
                        if fo:
                            calc.dumpMeasure(fo, calc.signalType, calc.measure) 
                    if fo:
                        calc.dumpMeasure(fo, SignalType.Aggr, FFT_Order.ssssAggr[order][sinSource][noiseType][noise])
        if fo:
            print(f'{datetime.datetime.now()}: Finish dump to {dumpOrderPath}')
            fo.close()

    @staticmethod
    def readCosSin(testCase:unittest.TestCase, order:int, sinSource:SinSource, fr) -> tuple[varDbl.VarDbl]:
        size = 1 << order
        sStep = [0] * (size << 1)
        for imag in (0, 1):
            for val in (1, 0):
                line = next(fr)
                sWord = line.split('\t')
                testCase.assertEqual(FFT_Step.dataOffset() + size, len(sWord))
                testCase.assertEqual(sinSource, SinSource(sWord[0]))
                testCase.assertEqual('0', sWord[2])
                testCase.assertEqual(order, int(sWord[4]))
                testCase.assertEqual('CosSin', sWord[6])
                for i in (1,3,5,7):
                    testCase.assertEqual('', sWord[i])
                testCase.assertEqual(imag, int(sWord[8]))
                testCase.assertEqual(val, int(sWord[9]))
                for i in range(size):
                    try:
                        if val:
                            sStep[(i << 1) + imag] = float(sWord[i + FFT_Step.dataOffset()])
                        else:
                            sStep[(i << 1) + imag] = varDbl.VarDbl(sStep[(i << 1) + imag], float(sWord[i + FFT_Step.dataOffset()]))
                    except ValueError:
                        testCase.fail(f'Invalid imag={imag} val={val} [{i}]={sWord[i + FFT_Step.dataOffset()]} in CosSin line: {sWord}')
        return sStep
    
    @staticmethod
    def readFFTOrderParam(testCase:unittest.TestCase, sWord:tuple[str], dumpStepPath:str) -> tuple[FFT_Order_Param]:
        try:
            sinSource = SinSource(sWord[0])
            noiseType = NoiseType(sWord[1])
            noise = float(sWord[2])
            signalType = SignalType(sWord[3])
            order = int(sWord[4])
            freq = int(sWord[5])
            signal = FFT_Signal_Param(sinSource, signalType, order, freq)
            return FFT_Order_Param(signal, noiseType, noise)
        except ValueError as ex:
            testCase.fail(f'Invalid parameter {ex} for {dumpStepPath}: {sWord}')

    @staticmethod
    def readStep(testCase:unittest.TestCase, order:int, sinSource:SinSource, 
                 dumpStepPath:str, fr, param:FFT_Order_Param, test:TestType, step:typing.Union[str, int]) \
            -> tuple[tuple[varDbl.VarDbl], FFT_Order_Param]:
        '''
        Read one step from {fr}, which is the file reader of {dumpStepPath}.
        If {param} is None, read the parameters from the first line.  Otherwise, verify the parameters.
        Return the step data and the parameters.
        If the first line is missing, and {test} is Forward, {step} is 'Input', return None, None.
        '''
        size = 1 << order
        sStep = [0] * (size << 1)
        for imag in (0, 1):
            for val in (1, 0):
                try:
                    line = next(fr)
                except StopIteration:
                    if test == TestType.Forward and step == 'Input' and imag == 0 and val == 1:
                        return None, None
                    else:
                        testCase.fail(f'Invalid empty line for {dumpStepPath} at test={test} step={step}, imag={imag}, val={val}')
                sWord = line.split('\t')
                if not param:
                    testCase.assertEqual(step, 'Input')
                    param = FFT_Step.readFFTOrderParam(testCase, sWord, dumpStepPath)
                    testCase.assertEqual(param.signal.order, order)
                    testCase.assertEqual(param.signal.sinSource, sinSource)
                else:
                    try:
                        testCase.assertEqual(param.signal.sinSource, SinSource(sWord[0]))
                        testCase.assertEqual(param.noiseType,  NoiseType(sWord[1]))
                        testCase.assertEqual(param.noise, float(sWord[2]))
                        testCase.assertEqual(param.signal.signalType, SignalType(sWord[3]))
                        testCase.assertEqual(param.signal.order, int(sWord[4]))
                        testCase.assertEqual(param.signal.freq, int(sWord[5]))
                        testCase.assertEqual(test, TestType(sWord[6]))
                        testCase.assertEqual(step, sWord[7])
                        testCase.assertEqual(imag, int(sWord[8]))
                        testCase.assertEqual(val, int(sWord[9]))
                    except (AssertionError, ValueError) as ex:
                        testCase.fail(f'Invalid parameter {ex} for {dumpStepPath}: {line}')
                for i in range(size):
                    try:
                        if val:
                            sStep[(i << 1) + imag] = float(sWord[i + FFT_Step.dataOffset()])
                        else:
                            sStep[(i << 1) + imag] = varDbl.VarDbl(sStep[(i << 1) + imag], float(sWord[i + FFT_Step.dataOffset()]))
                    except ValueError:
                        testCase.fail(f'Invalid imag={imag} val={val} [{i}]={sWord[i + FFT_Step.dataOffset()]} in {step} line for {dumpStepPath}: {sWord}')
        return tuple(sStep), param
                
    @staticmethod
    def readSteps(testCase:unittest.TestCase, order:int, sinSource:SinSource, 
                  dumpStepPath:str, fr, test:TestType) -> tuple[tuple[varDbl.VarDbl], tuple[tuple[varDbl.VarDbl]], FFT_Signal, NoiseType, float]:
        sInput, param = FFT_Step.readStep(testCase, order, sinSource, dumpStepPath, fr, None, test, 'Input')
        if not sInput:
            return None, None, None
        ssStep = []
        for i in range(order + len(FFT_Step.EXTRA_STEPS) + 1):
            sStep, _ = FFT_Step.readStep(testCase, order, sinSource, dumpStepPath, fr, param, test, str(i))
            ssStep.append(sStep)
        return param, sInput, ssStep

    @staticmethod
    def assertStep(testCase:unittest.TestCase, context:str, order:int, step:int,
                   sData1:tuple[varDbl.VarDbl], sData2:tuple[varDbl.VarDbl],
                   precDiff=0):
        '''
        compare {sData1} and {sData2} item by item.
        If {precDiff} > 0, the difference must be within {precDiff} times the uncertainty.
        If {precDiff} == 0, the value and uncertainty must be exactly the same.
        If {precDiff} < 0, the value must be exactly the same.
        '''
        testCase.assertEqual(len(sData1), len(sData2))
        for i in range(len(sData1)):
            if precDiff > 0:
                diff = sData1[i] - sData2[i]
                try:
                    testCase.assertLessEqual(abs(diff.value()), precDiff * diff.uncertainty())
                except AssertionError as ex:
                    print(f'Diff at order={order} step={step} {context} [{i}] {sData1[i]} vs {sData2[i]}')
                    raise ex
            elif precDiff == 0:
                try:
                    varDbl.assertVarDblEqual(testCase, sData1[i], sData2[i])
                except AssertionError as ex:
                    diff =  sData1[i] - sData2[i]
                    try:
                        testCase.assertLessEqual(abs(diff.value()), math.ulp(1))
                    except AssertionError as ex:
                        print(f'Diff at order={order} step={step} {context} [{i}] {sData1[i]} vs {sData2[i]}')
                        raise ex
            else:
                try:
                    testCase.assertAlmostEqual(sData1[i].value(), sData2[i].value(), 
                                                delta=math.ulp(1) * (sData1[i].value() + sData2[i].value()))
                except AssertionError as ex:
                    print(f'Diff value at order={order} step={step} {context} [{i}] {sData1[i]} vs {sData2[i]}')
                    raise ex
                try:
                    testCase.assertAlmostEqual(sData1[i].uncertainty(), sData2[i].uncertainty(),
                                                delta=math.ulp(1) * (sData1[i].uncertainty() + sData2[i].uncertainty()))
                except AssertionError as ex:
                    print(f'Diff uncertainty at order={order} step={step} {context} [{i}] {sData1[i]} vs {sData2[i]}')
                    raise ex


    @staticmethod
    def assertSteps(testCase:unittest.TestCase, context:str, order:int,
                    ssData1:tuple[tuple[varDbl.VarDbl]], ssData2:tuple[tuple[varDbl.VarDbl]],
                    precDiff=0):
        for i in range(order + len(FFT_Step.EXTRA_STEPS) + 1):
            FFT_Step.assertStep(testCase, context, order, i, ssData1[i], ssData2[i], precDiff=precDiff)

    @staticmethod                
    def assertFFTOrderParam(testCase:unittest.TestCase, fft1:FFT_Order_Param, fft2:FFT_Order_Param):
        testCase.assertEqual(fft1.signal.sinSource, fft2.signal.sinSource)
        testCase.assertEqual(fft1.noiseType, fft2.noiseType)
        testCase.assertEqual(fft1.noise, fft2.noise)
        testCase.assertEqual(fft1.signal.signalType, fft2.signal.signalType)
        testCase.assertEqual(fft1.signal.order, fft2.signal.order)
        testCase.assertEqual(fft1.signal.freq, fft2.signal.freq)


    @staticmethod
    def recalc(testCase:unittest.TestCase, order:int, sinSource: SinSource, 
                dumpStepPath:str=None, 
                dumpOrderPath=f'./Python/Output/FFT_{IndexSin.MIN_ORDER}_{IndexSin.MAX_ORDER}.txt',
                uncPrec=1e-2, errPrec=1, precDiff=0,
                minCount=FFT_Order.MIN_COUNT, minNoise=1e-14,
                verbose=False):
        '''
        Test the FFT step from the dump file {dumpStepPath} by recalculation FFT_Step using data from the dump file, 
            to verify the floating arithmetic from different languages.

        Compare the individual step calculation with the aggregated result {dumpOrderPath}.
            The error deviation must be within {errPrec} and the uncertainty mean within {uncPrec}.
            Only test the error deviation for noise >= {minNoise} and the count > minCount.
        '''
        IndexSin.validateOrder(order)
        size = 1 << order

        if not dumpStepPath:
            dumpStepPath = FFT_Step.dumpPath(order, sinSource)
        if not os.path.isfile(dumpStepPath):
            testCase.fail(f'Invalid dumpPath {dumpStepPath} for reading FFT data')

        if dumpOrderPath:
            sssssAggr = FFT_Order.read(dumpOrderPath, SignalType.Aggr, 0, checkEnding=False)
            sssssLine = FFT_Order.read(dumpOrderPath, SignalType.Linear, 0, checkEnding=False)
        else:
            sssssAggr = {}
            sssssLine = {}

        def assertIndexSin(testCase, idxSin:IndexSin):
            testCase.assertEqual(idxSin.sinSource, SinSource.Limit)
            half = size >> 1
            for i in range(size):
                cos = idxSin.cos(i, order)
                sin = idxSin.sin(i + half, order)
                testCase.assertEqual(cos.value(), sin.value())
                testCase.assertEqual(cos.uncertainty(), sin.uncertainty())
        
        with open(dumpStepPath) as fr:
            hdr = next(fr)      
            size = len(hdr.split('\t')) - FFT_Step.dataOffset()
            testCase.assertEqual(size, 1 << order)
            testCase.assertEqual(hdr.strip(), FFT_Step.header(order))
            sCosSin = FFT_Step.readCosSin(testCase, order, sinSource, fr)

            while True:
                forward, sFrwd, ssSpecStep = FFT_Step.readSteps(testCase, order, sinSource, dumpStepPath, fr, TestType.Forward)
                if not forward:
                    break
                print(f'{datetime.datetime.now()}: Start assert order={order} sinSource={sinSource} signal={forward.signal.signalType} freq={forward.signal.freq} noiseType={forward.noiseType} noise={forward.noise}')
                roundtrip, sRound, ssRoundStep = FFT_Step.readSteps(testCase, order, sinSource, dumpStepPath, fr, TestType.Roundtrip)
                FFT_Step.assertFFTOrderParam(testCase, forward, roundtrip)
                FFT_Step.assertStep(testCase, f'signal={forward.signal.signalType} freq={forward.signal.freq} noise={forward.noise} roundtrip', 
                                    order, order + 1, sRound, ssSpecStep[order + 1])
                reverse, sBack, ssRevStep = FFT_Step.readSteps(testCase, order, sinSource, dumpStepPath, fr, TestType.Reverse)
                FFT_Step.assertFFTOrderParam(testCase, forward, reverse)
                sFreq = ssSpecStep[order + 2]
                sWave = ssRevStep[order + 2]

                fft = FFT_Step(forward.signal, forward.noiseType, forward.noise,
                               sCosSin=sCosSin, sWave=sWave, sFreq=sFreq, sFrwd=sFrwd, sBack=sBack)
                testCase.assertEqual(len(fft.ssSpecStep), order + len(FFT_Step.EXTRA_STEPS) + 1)
                testCase.assertEqual(len(fft.ssRoundStep), order + len(FFT_Step.EXTRA_STEPS) + 1)
                testCase.assertEqual(len(fft.ssRevStep), order + len(FFT_Step.EXTRA_STEPS) + 1)
                assertIndexSin(testCase, fft.idxSin)
                if verbose:
                    print(f'Start assert for order={order} sinSource={sinSource} signal={forward.signal.signalType} freq={forward.signal.freq} noiseType={forward.noiseType} noise={forward.noise}')
                FFT_Step.assertSteps(testCase, f'signal={forward.signal.signalType} freq={forward.signal.freq} noise={forward.noise} forward', 
                                     order, ssSpecStep, fft.ssSpecStep, precDiff=precDiff)
                FFT_Step.assertSteps(testCase, f'signal={forward.signal.signalType} freq={forward.signal.freq} noise={forward.noise} roundtrip', 
                                     order, ssRoundStep, fft.ssRoundStep, precDiff=precDiff)
                FFT_Step.assertSteps(testCase, f'signal={forward.signal.signalType} freq={forward.signal.freq} noise={forward.noise} reverse', 
                                     order, ssRevStep, fft.ssRevStep, precDiff=precDiff)


                sssssOrder = sssssLine if forward.signal.signalType == SignalType.Linear else sssssAggr
                if sssssOrder and (ssssOrder := sssssOrder.get(order)) and (sssOrder := ssssOrder.get(sinSource)) \
                        and (ssOrder := sssOrder.get(forward.noiseType)) and (sOrder := ssOrder.get(forward.noise)) and (len(sOrder) == 3):
                    for test in TestType:
                        errDev, uncMean, cntErr = sOrder[test]
                        try:
                            testCase.assertLessEqual(fft.measure.sHisto[test].stat().count(), cntErr)
                        except AssertionError as ex:
                            raise ex
                        if forward.noise < minNoise:
                            continue
                        try:
                            unc = fft.measure.sUncStat[test].mean()
                            testCase.assertAlmostEqual(unc /uncMean, 1, delta = uncPrec)
                        except AssertionError as ex:
                            raise ex
                        if minCount <= cntErr and test != TestType.Roundtrip and 0 < forward.noise:
                            try:
                                err = fft.measure.sHisto[test].stat().dev()
                                testCase.assertAlmostEqual(err / errDev, 1, delta = errPrec)
                            except AssertionError as ex:
                                raise ex

                                     

    @staticmethod
    def compare(testCase:unittest.TestCase, sinSource:SinSource, dumpStepPath1:str, dumpStepPath2:str, precDiff= 1):
        '''
        Direct compare the FFT step result from different folder {dumpDir}
        '''
        print(f'{datetime.datetime.now()}: Start compare {dumpStepPath1} vs {dumpStepPath2}')
        with open(dumpStepPath1) as f1, open(dumpStepPath2) as f2:
            hdr = next(f1)      
            size = len(hdr.split('\t')) - FFT_Step.dataOffset()
            order = IndexSin.validateSize(size)
            testCase.assertEqual(hdr, next(f2))
            sCosSin1 = FFT_Step.readCosSin(testCase, order, sinSource, f1)
            sCosSin2 = FFT_Step.readCosSin(testCase, order, sinSource, f2)
            FFT_Step.assertStep(testCase, 'CosSin', order, 'CosSin', sCosSin1, sCosSin2, precDiff=0)

            while True:
                forward1, sFrwd1, ssSpecStep1 = FFT_Step.readSteps(testCase, order, sinSource, dumpStepPath1, f1, TestType.Forward)
                if not forward1:
                    break
                print(f'{datetime.datetime.now()}: Start compare order={order} sinSource={sinSource} signal={forward1.signal.signalType} freq={forward1.signal.freq} noiseType={forward1.noiseType} noise={forward1.noise}')    
                testCase.assertEqual(forward1.signal.order, order)
                context = f'signal={forward1.signal.signalType} freq={forward1.signal.freq} noise={forward1.noise} forward'
                forward2, sFrwd2, ssSpecStep2 = FFT_Step.readSteps(testCase, order, sinSource, dumpStepPath2, f2, TestType.Forward)
                FFT_Step.assertFFTOrderParam(testCase, forward1, forward2)
                FFT_Step.assertStep(testCase, context, order, 'sFrwd', sFrwd1, sFrwd2, precDiff=precDiff)
                FFT_Step.assertStep(testCase, context, order, 'sFreq', ssSpecStep1[-2], ssSpecStep2[-2], precDiff=-1)
                FFT_Step.assertSteps(testCase, context, order, ssSpecStep1, ssSpecStep2, precDiff=precDiff)

                context = f'signal={forward1.signal.signalType} freq={forward1.signal.freq} noise={forward1.noise} roundtrip'
                roundtrip1, sRound1, ssRoundStep1 = FFT_Step.readSteps(testCase, order, sinSource, dumpStepPath1, f1, TestType.Roundtrip)
                roundtrip2, sRound2, ssRoundStep2 = FFT_Step.readSteps(testCase, order, sinSource, dumpStepPath2, f2, TestType.Roundtrip)
                FFT_Step.assertFFTOrderParam(testCase, forward1, roundtrip1)
                FFT_Step.assertFFTOrderParam(testCase, forward2, roundtrip2)
                FFT_Step.assertStep(testCase, context, order, 'sRound', sRound1, sRound2, precDiff=precDiff)
                FFT_Step.assertStep(testCase, context, order, 'sFrwd', ssRoundStep1[-2], ssRoundStep2[-2], precDiff=precDiff)
                FFT_Step.assertSteps(testCase, context, order, ssRoundStep1, ssRoundStep2, precDiff=precDiff)
                
                context = f'signal={forward1.signal.signalType} freq={forward1.signal.freq} noise={forward1.noise} reverse'
                reverse1, sBack1, ssRevStep1 = FFT_Step.readSteps(testCase, order, sinSource, dumpStepPath1, f1, TestType.Reverse)
                reverse2, sBack2, ssRevStep2 = FFT_Step.readSteps(testCase, order, sinSource, dumpStepPath2, f2, TestType.Reverse)
                FFT_Step.assertFFTOrderParam(testCase, forward1, reverse1)
                FFT_Step.assertFFTOrderParam(testCase, forward2, reverse2)
                FFT_Step.assertStep(testCase, context, order, 'sBack', sBack1, sBack2, precDiff=precDiff)
                FFT_Step.assertStep(testCase, context, order, 'sWave', ssRevStep1[-2], ssRevStep2[-2], precDiff=-1)
                FFT_Step.assertSteps(testCase, context, order, ssRevStep1, ssRevStep2, precDiff=precDiff)




