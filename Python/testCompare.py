'''
Compare Taylor expansion from different implementation
'''
import datetime
import math
import os
import unittest
import re
import sys
import traceback

from histo import Histo, Stat
from taylor import Taylor
from varDbl import VarDbl, assertVarDblEqual
from testManual import TestConvergence
from fft import SinSource, SignalType, TestType, FFT_Order, FFT_Step


class TestIndexSine (unittest.TestCase):
    '''
    Compare the indexed sine between Python, Java and C++
    It shows that Java and C++ are almost identical while Python differs slightly due to different sine implementation
    '''
    def diff(self, dumpPathBase:str, dumpPathCode:str, dumpPathDiff:str, context:str, dev:float,
             compareUncertainty:bool=True, denom=2**18):
        '''
        compare the indexed sine dump file between {dumpPathBase} and {dumpPathCode}, 
        and dump the diff histogram to {dumpPathDiff}
        '''
        histo = Histo(5, 3)
        sDiff = []
        with open(dumpPathBase) as fb, open(dumpPathCode) as fc:
            self.assertEqual(next(fb), next(fc))
            ln = 0
            for lineBase, lineCode in zip(fb, fc):
                ln += 1
                sWordBase = tuple(map(float, lineBase.split('\t')))
                sWordCode = tuple(map(float, lineCode.split('\t')))
                self.assertTupleEqual(sWordBase[:2], sWordCode[:2])
                sin1 = VarDbl(sWordCode[2], sWordCode[3])
                sin2 = VarDbl(sWordBase[2], sWordBase[3])
                if compareUncertainty:
                    self.assertEqual(sin1.uncertainty(), sin2.uncertainty())
                err = sin1 - sin2
                if err.uncertainty() > 0:
                    norm = err.value() /err.uncertainty()
                    if compareUncertainty:
                        norm *= math.sqrt(2/3)  # to LSB
                    histo.accum(norm, ln)
                    if norm:
                        sDiff.append((ln - 1, norm))
                else:
                    self.assertEqual(err.value(), 0)
                   
        with open(dumpPathDiff, 'w') as fw:
            fw.write(f'Context\tCount\tMean\tDev\tMin\tMin At\tMax\tMax At\tLower\tUpper')
            for b in histo.buckets():
                fw.write(f'\t{b}')
            fw.write('\n')
            stat = histo.stat()
            fw.write(f'{context}\t{stat.count()}\t{stat.mean()}\t{stat.dev()}'
                     f'\t{stat.min()}\t{stat.minAt()}\t{stat.max()}\t{stat.maxAt()}'
                     f'\t{histo.less()}\t{histo.more()}')
            for c in histo.histogram():
                fw.write(f'\t{c / stat.count()}')
            fw.write('\nIndex\tNorm\n')
            for ln, norm in sDiff:
                fw.write(f'{ln/denom}\t{norm}\n')
        self.assertAlmostEqual(histo.stat().dev(), dev)
            
    def test_cpp_vs_python(self):
        dumpPathBase = './Python/Output/IndexSin_Quart_18.txt'
        dumpPathCode = './Cpp/Output/IndexSin_Quart_18.txt'
        dumpPathDiff = './Python/Output/IndexSin_Quart_18_Cpp_Python.txt'
        self.diff(dumpPathBase, dumpPathCode, dumpPathDiff, 'Quart C++ vs Python', 0.1815195,)

    def test_java_vs_python(self):
        dumpPathBase = './Python/Output/IndexSin_Quart_18.txt'
        dumpPathCode = './Java/Output/IndexSin_Quart_18.txt'
        dumpPathDiff = './Python/Output/IndexSin_Quart_18_Java_Python.txt'
        self.diff(dumpPathBase, dumpPathCode, dumpPathDiff, 'Quart Java vs Python', 0.1821751)

    def test_cpp_vs_java(self):
        dumpPathBase = './Cpp/Output/IndexSin_Quart_18.txt'
        dumpPathCode = './Java/Output/IndexSin_Quart_18.txt'
        dumpPathDiff = './Python/Output/IndexSin_Quart_18_Java_Cpp.txt'
        self.diff(dumpPathBase, dumpPathCode, dumpPathDiff, 'Quart Java vs C++', 0.0354802)

    def test_quart_vs_prec_python(self):
        dumpPathBase = './Python/Output/IndexSin_Quart_18.txt'
        dumpPathCode = './Cpp/Output/IndexSin_Prec_18.txt'
        dumpPathDiff = './Python/Output/IndexSin_Quart_Prec_18_Python.txt'
        self.diff(dumpPathBase, dumpPathCode, dumpPathDiff, 'Quart vs Prec Python', 0.6431504,
                  compareUncertainty=False)

    def test_quart_vs_prec_cpp(self):
        dumpPathBase = './Cpp/Output/IndexSin_Quart_18.txt'
        dumpPathCode = './Cpp/Output/IndexSin_Prec_18.txt'
        dumpPathDiff = './Python/Output/IndexSin_Quart_Prec_18_Cpp.txt'
        self.diff(dumpPathBase, dumpPathCode, dumpPathDiff, 'Quart vs Prec C++', 0.6207008,
                  compareUncertainty=False)

    def test_quart_vs_prec_java(self):
        dumpPathBase = './Java/Output/IndexSin_Quart_18.txt'
        dumpPathCode = './Cpp/Output/IndexSin_Prec_18.txt'
        dumpPathDiff = './Python/Output/IndexSin_Quart_Prec_18_Java.txt'
        self.diff(dumpPathBase, dumpPathCode, dumpPathDiff, 'Quart vs Prec Java', 0.6210673,
                  compareUncertainty=False)


class TestDumpPath (unittest.TestCase):

    def assertExtensionEqual(self, l, r, valPrec=1e-6, uncPrec=1e-6):
        self.assertEqual(l.monotonics, r.monotonics)
        if r.exp == 0:
            self.assertEqual(l.exp, 0)
        else:
            self.assertAlmostEqual(l.exp / r.exp, 1)
        if r.momentum == 0:
            self.assertEqual(l.momentum, 0) 
        else:
            self.assertAlmostEqual(l.momentum /r.momentum, 1)
        assertVarDblEqual(self, l.taylor, r.taylor, valPrec, uncPrec)
        assertVarDblEqual(self, l.val, r.val, valPrec, uncPrec)
        assertVarDblEqual(self, l.var, r.var, valPrec, uncPrec)
        assertVarDblEqual(self, l.newVal, r.newVal, valPrec, uncPrec)
        assertVarDblEqual(self, l.newVar, r.newVar, valPrec, uncPrec)

    def dumpPath(self, pyPath, javaPath, valPrec=1e-6, uncPrec=1e-6):
        self.assertTrue(os.path.isfile(pyPath))
        self.assertTrue(os.path.isfile(javaPath))
        psIn, psExt, pRes = Taylor.verifyDumpFile(self, pyPath)
        jsIn, jsExt, jRes = Taylor.verifyDumpFile(self, javaPath)

        self.assertEqual(set(psIn), set(jsIn))
        for k in psIn:
            if type(psIn[k]) == str:
                pass
            elif type(psIn[k]) == VarDbl:
                assertVarDblEqual(self,psIn[k], jsIn[k], valPrec, uncPrec)
            elif type(psIn[k]) == float:
                self.assertAlmostEqual(psIn[k], jsIn[k])
            else:
                self.assertEqual(psIn[k], jsIn[k])

        self.assertEqual(len(psExt), len(jsExt))
        for n in range(1, len(psExt)):
            pExt = psExt[n]
            jExt = jsExt[n]
            try:
                self.assertExtensionEqual(pExt, jExt, valPrec, uncPrec)
            except AssertionError as ex:
                raise ex

        if type(pRes) == VarDbl:
            assertVarDblEqual(self, pRes, jRes, valPrec, uncPrec)
        elif type(pRes) == str:
            self.assertEqual(pRes, jRes)
        else:
            self.fail(f'Unknown result type {type(pRes)}: {pRes}')
    

    def test_inversion(self):
        self.dumpPath('./Python/Output/Pow_1_0.2_-1.txt', 
                      './Java/Output/pow_1.0_0.2_-1.0.txt')
        self.dumpPath('./Python/Output/Pow_1_0.2_-1.txt', 
                      './Cpp/Output/pow_1_0.2_-1.txt')
        
    def test_inversion_NotMonotonic(self):
        self.dumpPath('./Python/Output/Pow_1_0.2_-2.txt', 
                      './Java/Output/pow_1.0_0.2_-2.0.txt')
        self.dumpPath('./Python/Output/Pow_1_0.2_-2.txt', 
                      './Cpp/Output/pow_1_0.2_-2.txt')
        
    def test_exp(self):
        self.dumpPath('./Python/Output/exp_1_0.1.txt',
                      './Java/Output/exp_1.0_0.1.txt')
        self.dumpPath('./Python/Output/exp_1_0.1.txt',
                      './Cpp/Output/exp_1_0.1.txt')

    def test_log(self):
        self.dumpPath('./Python/Output/log_1_0.1.txt',
                      './Java/Output/log_1.0_0.1.txt')
        self.dumpPath('./Python/Output/log_1_0.1.txt',
                      './Cpp/Output/log_1_0.1.txt')

    def test_sin_0(self):
        self.dumpPath('./Python/Output/sin_0_0.1.txt',
                      './Java/Output/sin_0_0.1.txt')
        self.dumpPath('./Python/Output/sin_0_0.1.txt',
                      './Cpp/Output/sin_0_0.1.txt')
        
        self.dumpPath('./Python/Output/sin_0_1.0.txt',
                      './Java/Output/sin_0_1.txt')
        self.dumpPath('./Python/Output/sin_0_1.0.txt',
                      './Cpp/Output/sin_0_1.txt')

    def test_sin_half_pi(self):
        '''
        math.sin(pi/2) and math.cos(pi/2) seems to have larger rounding error
        '''
        self.dumpPath('./Python/Output/sin_0.5_0.1.txt',
                      './Java/Output/sin_0.5_0.1.txt')
        self.dumpPath('./Python/Output/sin_0.5_0.1.txt',
                      './Cpp/Output/sin_1.5708_0.1.txt',
                      valPrec=1e-4, uncPrec=1e-4)
        
        self.dumpPath('./Python/Output/sin_0.5_1.txt',
                      './Java/Output/sin_0.5_1.txt')
        self.dumpPath('./Python/Output/sin_0.5_1.txt',
                      './Cpp/Output/sin_1.5708_1.txt',
                      valPrec=1e-4, uncPrec=1e-4)

    def test_sin_quad_pi(self):
        self.dumpPath('./Python/Output/sin_0.25_0.1.txt',
                      './Java/Output/sin_0.25_0.1.txt')
        self.dumpPath('./Python/Output/sin_0.25_0.1.txt',
                      './Cpp/Output/sin_0.785398_0.1.txt')
        
        self.dumpPath('./Python/Output/sin_0.25_1.txt',
                      './Java/Output/sin_0.25_1.txt')
        self.dumpPath('./Python/Output/sin_0.25_1.txt',
                      './Cpp/Output/sin_0.785398_1.txt')


class TestConvergeEdge (unittest.TestCase):

    def convergeEdge(self, pyPath, javaPath, prec):
        self.assertTrue(os.path.isfile(pyPath))
        self.assertTrue(os.path.isfile(javaPath))
        with open(pyPath) as pf, open(javaPath) as jf:
            self.assertEqual(next(pf), TestConvergence.HEADER)
            self.assertEqual(next(jf), TestConvergence.HEADER)
            for pLn, jLn, in zip(pf.readlines(), jf.readlines()):
                psWd = pLn.split('\t')
                jsWd = jLn.split('\t')
                self.assertEqual(len(psWd), len(jsWd))
                psVal = map(float, psWd[:-1])
                jsVal = map(float, jsWd[:-1])
                for i, pv, jv in zip(range(len(psWd)), psVal, jsVal):
                    try:
                        if pv == 0:
                            self.assertLessEqual(abs(jv), sys.float_info.epsilon)
                        elif jv == 0:
                            self.assertLessEqual(abs(pv), sys.float_info.epsilon)
                        else:
                            self.assertAlmostEqual(pv / jv, 1, delta=prec)
                    except AssertionError as ex:
                        raise
        
    def test_PowEdge(self):
        self.convergeEdge('./Python/Output/PowEdge.txt', 
                          './Java/Output/PowEdge.txt',
                          prec=1e-9)
        self.convergeEdge('./Python/Output/PowEdge.txt', 
                          './Cpp/Output/PowEdge.txt',
                          prec=1e-9)

    def test_SinEdge(self):
        '''
        cos(PI/2)=6.12323399573676E-17 in Python differs from cos(PI/2)=6.12303176911189E-17 from Java or C++
        '''
        self.convergeEdge('./Python/Output/SinEdge.txt', 
                          './Java/Output/SinEdge.txt', 
                          prec=5e-3)
        self.convergeEdge('./Java/Output/SinEdge.txt', 
                          './Cpp/Output/SinEdge.txt', 
                          prec=5e-5)
        

class Test_FFT_Step_Prec (unittest.TestCase):

    def test_java(self):
        for order in range(2, 6):
            dumpPathPy = f'./Python/Output/FFT_Step_{order}_Prec.txt'
            dumpPathCode = f'./Java/Output/FFT_Step_{order}_Prec.txt'
            FFT_Step.compare(self, SinSource.Prec, order, dumpPathPy, dumpPathCode, precDiff=-1)

    def test_cpp(self):
        for order in range(2, 6):
            dumpPathPy = f'./Python/Output/FFT_Step_{order}_Prec.txt'
            dumpPathCode = f'./Cpp/Output/FFT_Step_{order}_Prec.txt'
            FFT_Step.compare(self, SinSource.Prec, order, dumpPathPy, dumpPathCode, precDiff=-1)


class Test_FFT_Step_Quart (unittest.TestCase):

    def test_java(self):
        sOrder_precDiff = {}
        for order in range(2, 6):
            dumpPathPy = f'./Python/Output/FFT_Step_{order}_Quart.txt'
            dumpPathCode = f'./Java/Output/FFT_Step_{order}_Quart.txt'
            for precDiff in range(-1, 10):
                try:
                    FFT_Step.compare(self, SinSource.Quart, order, dumpPathPy, dumpPathCode, precDiff=precDiff)
                    sOrder_precDiff[order] = precDiff
                    break
                except AssertionError as ex:
                    sOrder_precDiff[order] = ex
        self.assertDictEqual(sOrder_precDiff, {2: 0, 3: 4, 4: 5, 5: 5})

    def test_cpp(self):
        sOrder_precDiff = {}
        for order in range(2, 6):
            dumpPathPy = f'./Python/Output/FFT_Step_{order}_Quart.txt'
            dumpPathCode = f'./Cpp/Output/FFT_Step_{order}_Quart.txt'
            for precDiff in range(-1, 10):
                try:
                    FFT_Step.compare(self, SinSource.Quart, order, dumpPathPy, dumpPathCode, precDiff=precDiff)
                    sOrder_precDiff[order] = precDiff
                    break
                except AssertionError as ex:
                    sOrder_precDiff[order] = ex
        self.assertDictEqual(sOrder_precDiff, {2: 0, 3: 4, 4: 5, 5: 5})


class Test_FFT_Step_Lib (unittest.TestCase):

    def test_cos(self):
        py25 = VarDbl(   7.07106781186547572737e-01, 6.40987562127854728723e-17)
        java25 = VarDbl( 7.07106781186547600e-01,    6.40987562127854700e-17)
        cpp25 = VarDbl(  7.07106781186547572737e-01, 6.40987562127854728723e-17)

        py50 = VarDbl(   6.12323399573676603587e-17, 7.11639149972692318350e-33)
        java50 = VarDbl( 1.22464679914735320e-16,    1.42327829994538460e-32)
        cpp50 = VarDbl(  1.22460635382237725821e-16, 0.0)

        py75 = VarDbl(  -7.07106781186547461715e-01, 6.40987562127854728723e-17)
        java75 = VarDbl(-7.07106781186547500e-01,    6.40987562127854700e-17)
        cpp75 = VarDbl( -7.07106781186547461715e-01, 6.40987562127854728723e-17)


class Test_FFT_Order (unittest.TestCase):

    @staticmethod
    def printDiff(sssDiff):
        for sinSource in sssDiff:
            for test in sssDiff[sinSource]:
                for i, context in enumerate(['ErrDev', 'UncMean', 'Count']):
                    stat = sssDiff[sinSource][test][i]
                    print(f'{sinSource}\t{test}\t{context}\tmin={stat.min()}, {stat.minAt()}\tmax={stat.max()}, {stat.maxAt()}\tmean={stat.mean()}\tdev{stat.dev()}')

    def test_Cpp_Aggr(self):
        sssDiff = FFT_Order.compare(self, './Java/Output/FFT_2_19.txt', './Cpp/Output/FFT_2_19.txt')
        try:
            for sinSource in (SinSource.Prec, SinSource.Quart, SinSource.Lib):
                for test in TestType:
                    self.assertLess( sssDiff[sinSource][test][0].max(), 0.3 if test == TestType.Roundtrip else 0.1)
                    self.assertLess(-sssDiff[sinSource][test][0].min(), 0.3 if test == TestType.Roundtrip else 0.1)
                    self.assertLess(abs(sssDiff[sinSource][test][0].mean()), 2e-3 if test == TestType.Roundtrip else 1e-3)
                    self.assertLess(sssDiff[sinSource][test][1].max(), 1e-15)
                    self.assertLess(sssDiff[sinSource][test][2].max(), 1e-6)
        except AssertionError as ex:
            Test_FFT_Order.printDiff(sssDiff)
            raise ex

    def test_Cpp_Linear(self):
        sssDiff = FFT_Order.compare(self, './Java/Output/FFT_2_19.txt', './Cpp/Output/FFT_2_19.txt', SignalType.Linear)
        try:
            for sinSource in (SinSource.Prec, SinSource.Quart, SinSource.Lib):
                for test in TestType:
                    self.assertLess( sssDiff[sinSource][test][0].max(), 0.5 if test == TestType.Roundtrip else 0.4)
                    self.assertLess(-sssDiff[sinSource][test][0].max(), 0.5 if test == TestType.Roundtrip else 0.4)
                    self.assertLess(abs(sssDiff[sinSource][test][0].mean()), 6e-3 if test == TestType.Roundtrip else 3e-3)
                    self.assertLess(sssDiff[sinSource][test][1].max(), 1e-15)
                    self.assertLess(sssDiff[sinSource][test][2].max(), 1e-6)
        except AssertionError as ex:
            Test_FFT_Order.printDiff(sssDiff)
            raise ex


class Test_Exe_Time (unittest.TestCase):

    def readExeTime(self, logFile:str, regEx:str, 
                    idxOrder:int, idxNoise:int, idxSinSource:int=None):
        '''
        read execution time in the format of <time>: <regEx for order>
        '''
        pattern = re.compile(regEx)
        sExeTime:dict[int, Stat] = {}
        prevTime = None
        prevOrder = 0
        with open(logFile) as f:
            for line in f:
                if ': ' not in line:
                    prevTime = None
                    continue
                try:
                    sWord = line.split(': ')
                    time = datetime.datetime.fromisoformat(sWord[0])
                    match = pattern.match(sWord[1])
                    if not match:
                        prevTime = None
                        continue
                    order = int(match.group(idxOrder))
                    noise = float(match.group(idxNoise))
                    sinSource = None if idxSinSource is None else match.group(idxSinSource)
                except Exception:
                    prevTime = None
                    continue
                if prevTime:
                    if order not in sExeTime:
                        sExeTime[order] = Stat()
                    if sinSource:
                        sExeTime[prevOrder].accum((time - prevTime).total_seconds(), (noise, sinSource))
                    else:
                        sExeTime[prevOrder].accum((time - prevTime).total_seconds(), noise)
                prevTime = time
                prevOrder = order
        
        with open(logFile + '.txt', 'w') as f:
            f.write(f'Order\tCount\tMean\tDev\tMin\tMin At\tMax\tMax At\n')
            for order in sorted(sExeTime.keys()):
                stat = sExeTime[order]
                self.assertGreaterEqual(stat.min(), 0)
                self.assertLessEqual(stat.min(), stat.mean())
                self.assertGreaterEqual(stat.max(), stat.mean())
                f.write(f'{order}\t{stat.count()}\t{stat.mean()}\t{stat.dev()}\t{stat.min()}\t{stat.minAt()}\t{stat.max()}\t{stat.maxAt()}\n')

    def test_Adjugate(self):
        self.readExeTime('./Python/Output/testAdjugate.log', 
                          '^Start size=(\d+), noise=(\d+e-\d+|\d+.\d+|\d+)$', 1, 2)

    def test_FFT_Order_Cpp(self):
        self.readExeTime('./Cpp/Output/FFT_2_19.txt.log', 
                          '^Start calulation order=(\d+), sinSource=(Prec|Quart|Lib), noiseType=(Gaussian|White), noise=(\d+e-\d+|\d+.\d+|\d+)$',
                          1, 4, 2)
    
    def test_FFT_Order_Java(self):
        self.readExeTime('./Java/Output/FFT_2_19.txt.log', 
                          '^Starting noiseType=(Gaussian|White) noise=(\d+e-\d+|\d+.\d+|\d+) order=(\d+) sinSource=(Prec|Quart|Lib)$',
                          3, 2, 4)
    
    def test_FFT_Order_Python(self):
        self.readExeTime('./Python/Output/FFT_2_19.txt.log', 
                          '^Start calulation order=(\d+), sinSource=(Prec|Quart|Lib), noiseType=(Gaussian|White), noise=(\d+e-\d+|\d+.\d+|\d+)$',
                          1, 4, 2)


if __name__ == '__main__':
    unittest.main()