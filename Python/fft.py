
import abc
import enum
import math
import os
import random
import traceback

import indexSin
import histo
import varDbl


class FFTSinSource (enum.StrEnum):
    IndexSin = 'IndexSin',
    LibSin = 'LibSin',


class FFTBase (abc.ABC):
    MAX_ORDER = 19
    _bitReversedIndex = {}

    @abc.abstractmethod
    def source(self) ->FFTSinSource:
        '''
        sin to use for FFT phase
        '''

    @abc.abstractmethod
    def sin(self, freq:int, order:int) -> float:
        '''
        sin(Math.pi * freq / (1 << order))
        '''

    @abc.abstractmethod
    def cos(self, freq:int, order:int) -> float:
        '''
        cos(Math.pi * freq / (1 << order))
        '''

    @staticmethod
    def bitReversedIndices(order:int) -> tuple[int]:
        '''
        Return bit reversed indices
        from NumericalRecipesinC.pdf
        '''
        if sRes := FFTBase._bitReversedIndex.get(order):
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
        FFTBase._bitReversedIndex[order] = sRes
        return sRes
    
    def transform(self, sData:list, forward:bool, useOriginalArray:bool=False):
        '''
        "sData": an array of size (2<<order), with each datum contains (real, image)
        "forward": true for forware transformation, false for backward transformation.
        "useOriginalArray": change sData without allocate a result array
        '''
        for order in range(2, FFTBase.MAX_ORDER):
            if (2 << order) == len(sData):
                break
        if order > FFTBase.MAX_ORDER:
            raise RuntimeError(f'Invalid input array size {len(sData)} which is not 2^{order}')

        sRes = sData if useOriginalArray else [0] * (2 << order)
        
        sIndex = FFTBase.bitReversedIndices(order)
        for i in range(len(sIndex)):
            j = sIndex[i]
            sRes[i << 1], sRes[(i << 1) + 1] = sData[j << 1], sData[(j << 1) + 1]

        for i in range(0, len(sIndex), 2):
            j = i << 1
            sRes[j], sRes[j + 2] = sRes[j] + sRes[j + 2], sRes[j] - sRes[j + 2]
            sRes[j + 1], sRes[j + 3] = sRes[j + 1] + sRes[j + 3], sRes[j + 1] - sRes[j + 3]

        for o in range(2, order + 1):
            k = 1 << o
            for j in range(k >> 1):
                cos = self.cos(j, o)
                sin = self.sin(j, o) if forward else - self.sin(j, o)
                for i in range(0, len(sIndex), k):
                    i0 = (i + j) << 1
                    i1 = i0 + k
                    rd = sRes[i1] * cos - sRes[i1 + 1] * sin
                    id = sRes[i1] * sin + sRes[i1 + 1] * cos
                    sRes[i0], sRes[i1] = sRes[i0] + rd, sRes[i0] - rd
                    sRes[i0 + 1], sRes[i1 + 1] = sRes[i0 + 1] + id, sRes[i0 + 1] - id

        if not forward:
            sz = 1/len(sIndex)
            for i in range(len(sRes)):
                sRes[i] *= sz

        return sRes
    

class FFTIndexSin (FFTBase):
    _indexSin = indexSin.IndexSin(FFTBase.MAX_ORDER)

    def source(self):
        return FFTSinSource.IndexSin
    
    def sin(self, freq:int, order:int):
        return FFTIndexSin._indexSin.sin(freq *(1 <<(FFTBase.MAX_ORDER - order + 1)))
    
    def cos(self, freq:int, order:int):
        return FFTIndexSin._indexSin.cos(freq *(1 <<(FFTBase.MAX_ORDER - order + 1)))


class FFTLibSin (FFTBase):
    def source(self):
        return FFTSinSource.LibSin
    
    def sin(self, freq:int, order:int):
        return math.sin(math.pi *freq /(1 << (order - 1)))
    
    def cos(self, freq:int, order:int):
        return math.cos(math.pi *freq /(1 << (order - 1)))
    

class FFTLimitedSin (FFTBase):
    _indexSin = indexSin.IndexSin(FFTBase.MAX_ORDER - 1)

    def source(self):
        return FFTSinSource.LimitedSin
    
    def sin(self, freq:int, order:int):
        idx = FFTLimitedSin._indexSin.get_index(freq)
        return math.sin(math.pi *idx /(1 << (order - 1))) if idx >= 0 else\
              -math.sin(math.pi *-idx /(1 << (order - 1)))
    
    def cos(self, freq:int, order:int):
        idx = FFTLimitedSin._indexSin.get_index(freq)
        return math.cos(math.pi *idx /(1 << (order - 1))) if idx >= 0 else\
              -math.sin(math.pi *idx /(1 << (order - 1)))


class SignalType (enum.StrEnum):
    Sin = 'Sin',
    Cos = 'Cos',
    Linear = 'Linear',
    Aggr = 'Aggr',

class NoiseType (enum.StrEnum):
    Gaussian = 'Gaussian',
    White = 'White',

class TestType (enum.StrEnum):
    Forward = 'Forward',
    Reverse = 'Reverse',
    Roundtrip = 'Roundtrip',

class Measure:
    def __init__(self, divids=5, devs=3) -> None:
        self.sStat = {
            TestType.Forward: histo.Stat(), 
            TestType.Reverse: histo.Stat(), 
            TestType.Roundtrip: histo.Stat()}
        self.sHisto = {
            TestType.Forward: histo.Histo(divids, devs), 
            TestType.Reverse: histo.Histo(divids, devs), 
            TestType.Roundtrip: histo.Histo(divids, devs)}

   


class FFTTest:
    DIVIDS = 5
    DEVS = 3

    ssssAggr = {}

    def __init__(self, sinSource:FFTSinSource, noiseType:NoiseType, noise:float,
                 signal:SignalType, order:int, freq:int,
                 divids=DIVIDS, devs=DEVS) -> None:
        self.sinSource = sinSource
        self.noiseType = noiseType
        self.noise = abs(noise)
        self.signal = signal
        self.order = order
        self.freq = freq
        self.size = 1 << order
        if (signal != SignalType.Linear) and (freq * 2 >= self.size):
            raise ValueError(f'Invalid freq {freq} for order {order}')
        self.sWave = [0.0] * (self.size << 1)
        self.sFreq = [0.0] * (self.size << 1)

        match sinSource:
            case FFTSinSource.IndexSin:
                self.fft = FFTIndexSin()
            case FFTSinSource.LibSin:
                self.fft = FFTLibSin()
            case _:
                raise ValueError(f'Invalid sinSource={sinSource}')

        match noiseType:
            case NoiseType.Gaussian:
                pass
            case NoiseType.White:
                pass
            case _:
                raise ValueError(f'Invalid noiseType={noiseType}')

        peak = self.size >> 1
        match signal:
            case SignalType.Sin:
                for i in range(self.size):
                    if freq == (self.size // 4):
                        match (i % 4):
                            case 1:
                                self.sWave[i << 1] = 1.0
                            case 3:
                                self.sWave[i << 1] = -1.0
                    else:
                        self.sWave[i << 1] = self.fft.sin(freq * i, order)
                self.sFreq[(freq << 1) + 1] = peak
                self.sFreq[((self.size - freq) << 1) + 1] = - peak
            case SignalType.Cos:
                for i in range(self.size):
                    if freq == (self.size // 4):
                        match (i % 4):
                            case 0:
                                self.sWave[i << 1] = 1.0
                            case 2:
                                self.sWave[i << 1] = -1.0
                    else:
                        self.sWave[i << 1] = self.fft.cos(freq * i, order)
                self.sFreq[freq << 1] = peak
                self.sFreq[(self.size - freq) << 1] = peak
            case SignalType.Linear:
                self.sFreq[0] = self.size * (self.size - 1) / 2
                for i in range(1, self.size):
                    self.sWave[i << 1] = i
                    self.sFreq[i << 1] = -peak
                    if order < FFTBase.MAX_ORDER:
                        self.sFreq[(i << 1) + 1] = -peak / self.fft.sin(i, order + 1) * self.fft.cos(i, order + 1)
                    else:
                        self.sFreq[(i << 1) + 1] = -peak / math.sin(math.pi * i / self.size) * math.cos(math.pi * i / self.size)
            case _:
                raise ValueError(f"Unknown signal={signal}")

        if noise == 0:
            self.sData = [varDbl.VarDbl(self.sWave[i]) for i in range(self.size << 1)]
            self.sBack = [varDbl.VarDbl(self.sFreq[i]) for i in range(self.size << 1)]
        else:
            self.sData = [varDbl.VarDbl(self.sWave[i] + self.getNoise(), self.noise) for i in range(self.size << 1)]
            self.sBack = [varDbl.VarDbl(self.sFreq[i] + self.getNoise(), self.noise) for i in range(self.size << 1)]

        self.sSpec = self.fft.transform(self.sData, True)
        self.sRound = self.fft.transform(self.sSpec, False)
        self.sRev = self.fft.transform(self.sBack, False)

        self.measure = Measure(divids, devs)
        if signal == SignalType.Linear:
            self.aggr = None
        else:
            FFTTest.ssssAggr.setdefault(sinSource, {}).setdefault(noiseType, {}).setdefault(noise, {}).setdefault(order, Measure())
            self.aggr = FFTTest.ssssAggr[sinSource][noiseType][noise][order]     

        for i in range(self.size << 1):
            try:
                unc1 = self.sSpec[i].uncertainty()
                self.measure.sStat[TestType.Forward].accum(unc1)
                if self.aggr:
                    self.aggr.sStat[TestType.Forward].accum(unc1)
                if unc1 > 0:
                    if type(self.sFreq[i]) == varDbl.VarDbl:
                        self.measure.sHisto[TestType.Forward].accum((self.sSpec[i].value() - self.sFreq[i].value())/unc1)
                        if self.aggr:
                            self.aggr.sHisto[TestType.Forward].accum((self.sSpec[i].value() - self.sFreq[i].value())/unc1)
                    else:                    
                        self.measure.sHisto[TestType.Forward].accum((self.sSpec[i].value() - self.sFreq[i])/unc1)
                        if self.aggr:
                            self.aggr.sHisto[TestType.Forward].accum((self.sSpec[i].value() - self.sFreq[i])/unc1)

                unc2 = self.sRound[i].uncertainty()
                self.measure.sStat[TestType.Roundtrip].accum(unc2)
                if self.aggr:
                    self.aggr.sStat[TestType.Roundtrip].accum(unc2)
                if unc2 > 0:
                    self.measure.sHisto[TestType.Roundtrip].accum((self.sRound[i].value() - self.sData[i].value())/unc2)
                    if self.aggr:
                        self.aggr.sHisto[TestType.Roundtrip].accum((self.sRound[i].value() - self.sData[i].value())/unc2)

                unc3 = self.sRev[i].uncertainty()
                self.measure.sStat.get(TestType.Reverse).accum(unc3)
                if self.aggr:
                    self.aggr.sStat.get(TestType.Reverse).accum(unc3)
                if unc3 > 0:
                    if type(self.sWave[i]) == varDbl.VarDbl:
                        self.measure.sHisto.get(TestType.Reverse).accum((self.sRev[i].value() - self.sWave[i].value())/unc3)
                        if self.aggr:
                            self.aggr.sHisto.get(TestType.Reverse).accum((self.sRev[i].value() - self.sWave[i].value())/unc3)
                    else:
                        self.measure.sHisto.get(TestType.Reverse).accum((self.sRev[i].value() - self.sWave[i])/unc3)
                        if self.aggr:
                            self.aggr.sHisto.get(TestType.Reverse).accum((self.sRev[i].value() - self.sWave[i])/unc3)
            except BaseException as ex:
                print(traceback.format_exc())
                raise ex

    def getNoise(self) -> float:
        match self.noiseType:
            case NoiseType.Gaussian:
                return random.gauss(0, self.noise)
            case NoiseType.White:
                return random.uniform(-self.noise*math.sqrt(3), self.noise*math.sqrt(3))
            case _:
                raise ValueError(f'Invalid noiseType={self.noiseType}')
    
    @staticmethod
    def dumpSpectrumHeader(fw):
        fw.write('SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq'
            '\tIndex\tImag\tWave\tSpec'
            '\tForward Error\tForward Uncertainty'
            '\tReverse Error\tReverse Uncertainty'
            '\tRoundtrip Error\tRoundtrip Uncertainty\n')

    def dumpSpectrum(self, fw):
        for i in range(self.size << 1):
            fw.write(f'{self.sinSource}\t{self.noiseType}\t{self.noise}\t{self.signal}\t{self.order}\t{self.freq}')
            fw.write(f'\t{i >> 1}\t{i % 2}\t{self.sWave[i]}\t{self.sFreq[i]}')
            if type(self.sFreq[i]) == varDbl.VarDbl:
                fw.write(f'\t{self.sSpec[i].value() - self.sFreq[i].value()}\t{self.sSpec[i].uncertainty()}')
            else:                    
                fw.write(f'\t{self.sSpec[i].value() - self.sFreq[i]}\t{self.sSpec[i].uncertainty()}')
            if type(self.sWave[i]) == varDbl.VarDbl:
                fw.write(f'\t{self.sRev[i].value() - self.sWave[i].value()}\t{self.sRev[i].uncertainty()}')
            else:                    
                fw.write(f'\t{self.sRev[i].value() - self.sWave[i]}\t{self.sRev[i].uncertainty()}')
            fw.write(f'\t{self.sRound[i].value() - self.sData[i].value()}\t{self.sRound[i].uncertainty()}\n')
                
    @staticmethod
    def dumpSpectra(fw, sOrder:tuple[int], sNoise:tuple[int]=[i*1e-16 for i in range(10)], sFreq=range(1, 8)):
        for sinSource in FFTSinSource:
            for noiseType in NoiseType:
                for noise in sNoise:
                    for order in sOrder:
                        for signal in SignalType:
                            if signal == SignalType.Aggr:
                                continue
                            if signal == SignalType.Linear:
                                fftTest = FFTTest(sinSource, noiseType, noise, signal, order, 0)
                                fftTest.dumpSpectrum(fw)
                                continue
                            freq = 1 << (order - 2)
                            fftTest = FFTTest(sinSource, noiseType, noise, signal, order, freq)
                            fftTest.dumpSpectrum(fw)
                            for freq in sFreq:
                                if freq == 1 << (order - 2):
                                    continue
                                fftTest = FFTTest(sinSource, noiseType, noise, signal, order, freq)
                                fftTest.dumpSpectrum(fw)
                        fw.flush()

    @staticmethod
    def title(divids, devs):
        return "SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest"\
               "\tUncertainty Mean\tUncertainty Deviation\tUncertainty Minimum\tUncertainty Maximum"\
               "\tError Mean\tError Deviation\tError Minimum\tError Maximum\t"\
               + '\t'.join([f'{i/divids}' for i in range(-devs * divids, devs * divids + 1)])\
               + '\n' 

    @staticmethod
    def dumpMeasure(fw, sinSource:FFTSinSource, noiseType:NoiseType, noise:float,
                    signal:SignalType, order:int, freq:int, measure:Measure):
        for test in TestType:
            try:
                fw.write(f'{sinSource}\t{noiseType}\t{noise}\t{signal}\t{order}\t{freq}\t{test}')
                stat = measure.sStat[test]
                fw.write(f'\t{stat.mean()}\t{stat.dev()}\t{stat.min()}\t{stat.max()}')
                stat = measure.sHisto[test].stat()
                fw.write(f'\t{stat.mean()}\t{stat.dev()}\t{stat.min()}\t{stat.max()}')
                count = sum([c for c in measure.sHisto[test].histogram()])
                for c in measure.sHisto[test].histogram():
                    if count:
                        fw.write(f'\t{c/count}')
                    else:
                        fw.write(f'\t{c}')
                fw.write('\n')
            except BaseException as ex:
                print(f'{sinSource}\t{noiseType}\t{noise}\t{signal}\t{order}\t{freq}\t{test}: {ex}')
                raise ex
            
    @staticmethod
    def dumpOrders(sOrder=[o for o in range(5, FFTBase.MAX_ORDER)], 
                   sNoise=[math.pow(10,-n-1) for n in range(16)] + [0]):
        path = f'./Python/Output/FFT_{min(sOrder)}_{max(sOrder)}.txt'
        exist = os.path.isfile(path)
        sssssAggr = {}

        def readExist():
            nonlocal exist
            with open(path) as f:
                title = next(f)
                if FFTTest.title(FFTTest.DIVIDS, FFTTest.DEVS) != title:
                    exist = False
                    return
                for line in f:
                    sinSource,noiseType,noise,signal,order,freq,test = line.split('\t')[:7]
                    if signal != "Aggr":
                        continue
                    try:
                        sinSource = FFTSinSource(sinSource)
                    except BaseException as ex:
                        exist = False
                        return
                    try:
                        noiseType = NoiseType(noiseType)
                    except BaseException as ex:
                        exist = False
                        return
                    noise = float(noise)
                    order = int(order)
                    if not (4 <= order < FFTBase.MAX_ORDER):
                        continue
                    freq = int(freq)
                    if freq != 0:
                        exist = False
                        return
                    try:
                        test = TestType(test)
                    except BaseException as ex:
                        exist = False
                        return
                    sssssAggr.setdefault(sinSource, {}).setdefault(noiseType, {}).setdefault(noise, {})\
                              .setdefault(order, set()).add(test)

        if exist:
            readExist()

        with open(path, 'a' if exist else 'w') as fw:
            if not exist:
                fw.write(FFTTest.title(FFTTest.DIVIDS, FFTTest.DEVS))
            for sinSource in FFTSinSource:
                for noiseType in NoiseType:
                    for noise in sNoise:
                        for order in sOrder:
                            if (ssssAggr := sssssAggr.get(sinSource)) and (sssAggr := ssssAggr.get(noiseType)) \
                                    and (ssAggr := sssAggr.get(noise)) and (sAggr := ssAggr.get(order)) and (len(sAggr) == 3):
                                continue
                            for signal in SignalType:
                                if signal == SignalType.Linear:
                                    fftTest = FFTTest(sinSource, noiseType, noise, signal, order, 0)
                                    FFTTest.dumpMeasure(fw, sinSource, noiseType, noise, signal, order, 0, fftTest.measure)
                                    continue
                                for freq in range(1, 8):
                                    if signal == SignalType.Aggr:
                                        continue
                                    fftTest = FFTTest(sinSource, noiseType, noise, signal, order, freq)
                                    FFTTest.dumpMeasure(fw, sinSource, noiseType, noise, signal, order, freq, fftTest.measure)
                            aggr = FFTTest.ssssAggr[sinSource][noiseType][noise][order]
                            FFTTest.dumpMeasure(fw, sinSource, noiseType, noise, SignalType.Aggr, order, 0, aggr)
                            fw.flush()



