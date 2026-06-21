"""Bounded moments for distributions used by Taylor expansion:
abstract Moment base class with Normal, Uniform, and Exponential
implementations returning the moment table truncated by the bounding factor.
Normal and Uniform are symmetric (odd central moments vanish); Exponential
is asymmetric (all odd central moments are nonzero).
"""
import abc
import datetime
import math
import os
import scipy.special
import scipy.stats
import sympy

from indexSin import OUTDIR

class Moment (abc.ABC):

    @property
    @abc.abstractmethod
    def bounding(self):
        pass

    @property
    @abc.abstractmethod
    def leakage(self):
        pass

    @property
    @abc.abstractmethod
    def maxOrder(self):
        pass

    @abc.abstractmethod
    def __getitem__(self, n:int) -> float:
        pass


class Normal (Moment):
    '''
    Calculate variance moment for the given {bounding}.
    Detect the {_maxOrder} for the bounding, which is 448 when {bounding}=5.
    '''
    __slots__ = ('_sMoment', '_maxOrder', '_bounding')

    @staticmethod
    def readPreciseNorm(filePath:str):
        sMoment = []
        with open(filePath) as f:
            hdr = next(f)
            sHdr = hdr.split('\t')
            if (len(sHdr) != 4) or (sHdr[0] != 'n') or (sHdr[1] != 'Moment') or (sHdr[2] != 'Bounding:'):
                raise ValueError(f'Invalid header "{hdr}" in {filePath}')
            bounding = float(sHdr[3])
            if not (0 < bounding <= 8):
                raise ValueError(f'Invalid bounding {bounding} in header "{hdr}" in {filePath}')
            for i, line in enumerate(f):
                n, mmt = map(float, line.split('\t'))
                if n != i * 2:
                    raise ValueError(f'Invalid index {n} at line {i+2} in {filePath}, expected {i*2}: {line}')
                sMoment.append(mmt)
        return bounding, sMoment

    @staticmethod
    def calcPreciseNorm(bounding:float=5.0, maxOrder:int=1000000):
        filePath = f'{OUTDIR}/Python/NormalMoment_{bounding}.txt'
        HEADER = f'n\tMoment\tBounding={bounding}'
        b, sMoment = Normal.readPreciseNorm(filePath)
        if b != bounding:
            sMoment = []
        x = sympy.symbols("x", is_real=True)
        # Normalized per Formula (2.2): ζ(n, κ) = ∫z^n ρ dz / ∫ρ dz.
        density_const = 1/sympy.sqrt(2*sympy.pi) * sympy.exp(-x**2/2)
        norm = sympy.integrate(density_const, (x, -bounding, bounding)).evalf()
        while ((n := len(sMoment)*2) < maxOrder):
            print(f'At {datetime.datetime.now()}, Calculate {n}/{maxOrder} moment')
            mmt = sympy.integrate(density_const * x**n, (x, - bounding, bounding))
            if not math.isfinite(mmt):
                break
            sMoment.append((mmt / norm).evalf())
            with open(filePath, 'w') as f:
                f.write(HEADER)
                for i, mmt in enumerate(sMoment):
                    f.write(f'{i*2}\t{mmt}\n')

    def __init__(self, bounding:float=5, maxOrder:int=1000000):
        self._bounding = bounding
        filePath = f'{OUTDIR}/Python/Output/NormalMoment_{bounding}.txt'
        HEADER = 'Order\tMoment\n'
        if os.path.isfile(filePath):
            try:
                with open(filePath) as f:
                    hdr = next(f)
                    if hdr != HEADER:
                        raise NotImplementedError(f'Invalid header {hdr} vs {HEADER}')
                    n = -1
                    prevVal = 0
                    sMoment = []
                    for line in f.readlines():
                        n += 1
                        nn, val = map(float, line.strip().split('\t'))
                        if nn != n:
                            raise NotImplementedError(f'Invalid index {nn} vs {n}')
                        if (n & 1) == 0:
                            if val <= 0:
                                raise NotImplementedError(f'Invalid value {val} for index={n}')
                            if val <= prevVal:
                                raise NotImplementedError(f'Invalid {val} vs {prevVal} for index={n}')
                            sMoment.append(val)
                            prevVal = val
                        else:
                            if val != 0:
                                raise NotImplementedError(f'Invalid value {val} for index={n}')
                    self._maxOrder = len(sMoment) * 2
                    self._sMoment = sMoment
                    return
            except BaseException as ex:
                os.remove(filePath)

        term = 2 * scipy.stats.norm.pdf(bounding) * self._bounding
        bounding2 = self._bounding**2
        sTerm = []
        for n in range(maxOrder):
            try:
                sTerm.append(term * (1/(2*n + 1)))
                if not math.isfinite(sTerm[-1]):
                    del sTerm[-1]
                    break
                term *= bounding2
            except:
                break
        self._maxOrder = n * 2
        self._sMoment = sTerm[:]
        for j in range(2, maxOrder):
            for i in range(n):
                sTerm[i] *= 1/(2*i - 1 + 2*j) * bounding2
                prev = self._sMoment[i]
                self._sMoment[i] += sTerm[i]
                if prev == self._sMoment[i]:
                    n = i
                    break
            if n == 0:
                break
        # Normalize per Formula (2.2): divide unnormalized moments by ∫ρ dz = 1 - leakage.
        # Snapshot the 0th moment first; otherwise the first loop iteration sets sMoment[0]
        # to 1 and the remaining moments would be divided by 1 (no-op).
        norm_factor = self._sMoment[0]
        for i in range(len(self._sMoment)):
            self._sMoment[i] = self._sMoment[i] / norm_factor
        with open(filePath, 'w') as f:
            f.write(HEADER)
            for n in range(self.maxOrder):
                mmt = self[n]
                f.write(f'{n}\t{mmt}\n')

    @property
    def bounding(self):
        return self._bounding

    @property
    def leakage(self):
        return 1 - scipy.special.erf(self.bounding/math.sqrt(2))

    @property
    def maxOrder(self):
        return self._maxOrder

    def __getitem__(self, n:int) -> float:
        if n < 0:
            return IndexError()
        if (n % 2) == 1:
            return 0
        n //= 2
        if n >= len(self._sMoment):
            return IndexError()
        return self._sMoment[n]


class Uniform:
    '''
    Pre-calculated variance moment for uniform distribution [-1, 1].
    '''
    __slots__ = ('_bounding', '_sMoment', '_maxOrder')

    def __init__(self, bounding=math.sqrt(3)):
        self._bounding = bounding
        self._sMoment = []
        fac = 1
        for n in range(10000):
            try:
                mmt = fac/(2*n + 1)
            except OverflowError:
                break
            if not math.isfinite(mmt):
                break
            self._sMoment.append(mmt)
            fac *= bounding**2
        self._maxOrder = len(self._sMoment)

    @property
    def bounding(self):
        return self._bounding

    @property
    def leakage(self):
        # Two-sided leakage for the Uniform on [-√3, √3] truncated to [-κ, κ]:
        # 1 - (κ/√3). Default κ=√3 gives leakage 0.
        return 1 - self.bounding / math.sqrt(3)

    @property
    def maxOrder(self):
        return self._maxOrder

    def __getitem__(self, n:int) -> float:
        if n < 0 or n >= self._maxOrder or (n % 2) == 1:
            return 0
        return self._sMoment[n >> 1]


class Exponential(Moment):
    '''
    Central moments of the standardized Exponential distribution truncated to
    [ρ, κ] with the mean-reverting bounding condition ζ(1, κ) = 0.

    Density on the standardized variable z = λx - 1 (where x ~ Exp(λ)):
        ρ(z) = exp(-(1 + z))   for z ∈ [-1, ∞), 0 otherwise.
    Mean-reverting bounding has the closed-form solution
        a ≡ ρ + 1 = -W_0(-b · e^(-b))     where b = κ + 1
    (W_0 is the principal Lambert W branch). For κ > 0, -b·e^(-b) ∈ (-1/e, 0),
    where W_0 is real-valued in (-1, 0); the trivial root W_{-1}(...) = -b is
    rejected.

    With J_n(ρ, κ) = ∫_ρ^κ z^n exp(-(1+z)) dz, integration by parts unrolls into
    a closed sum:
        J_n = n! · [ e^(-a) · Σ_{j=0}^{n} ρ^j / j!  -  e^(-b) · Σ_{j=0}^{n} κ^j / j! ]
        J_0 = e^(-a) - e^(-b)
        μ_n(κ) = J_n / J_0
    By construction μ_1(κ) = 0 (mean-reverting).

    Unlike Normal/Uniform, the Exponential is *asymmetric*: ζ(odd, κ) ≠ 0 for
    n ≥ 3. Callers that assume symmetric moments must be revised accordingly.
    '''
    __slots__ = ('_bounding', '_lower', '_sMoment', '_maxOrder')

    @staticmethod
    def _solveLower(b):
        """Mean-reverting lower bound a = -W_0(-b·e^(-b)) ∈ (0, 1). Requires b > 1."""
        if b <= 1.0:
            raise ValueError(f'Exponential bounding must give b=κ+1 > 1; got b={b}')
        w = scipy.special.lambertw(-b * math.exp(-b), k=0)
        # For arg ∈ (-1/e, 0), W_0 is real; the imaginary part is numerical noise.
        if abs(w.imag) > 1e-12 * max(abs(w.real), 1.0):
            raise ValueError(f'Lambert W returned non-real value for b={b}: {w}')
        return float(-w.real)

    def __init__(self, bounding=15.0, maxOrder=256):
        if bounding <= 0:
            raise ValueError(f'bounding must be positive, got {bounding}')
        self._bounding = bounding
        b = bounding + 1.0
        a = Exponential._solveLower(b)
        self._lower = a - 1.0
        rho = a - 1.0
        kappa = bounding
        e_neg_a = math.exp(-a)
        e_neg_b = math.exp(-b)
        norm = e_neg_a - e_neg_b
        # Build μ_n from the closed-form J_n = n! · (e^(-a)·rho_sum - e^(-b)·kap_sum),
        # where rho_sum and kap_sum are partial exponential series accumulated in-place.
        sMoment = []
        rho_sum = 0.0
        kap_sum = 0.0
        rho_j = 1.0     # ρ^j
        kap_j = 1.0     # κ^j
        fact_j = 1.0    # j!
        n_fact = 1.0    # n!
        for n in range(maxOrder + 1):
            rho_sum += rho_j / fact_j
            kap_sum += kap_j / fact_j
            J_n = n_fact * (e_neg_a * rho_sum - e_neg_b * kap_sum)
            mu_n = J_n / norm
            if not math.isfinite(mu_n):
                break
            sMoment.append(mu_n)
            rho_j *= rho
            kap_j *= kappa
            fact_j *= (n + 1)
            n_fact *= (n + 1)
        self._sMoment = sMoment
        self._maxOrder = len(sMoment)

    @property
    def bounding(self):
        return self._bounding

    @property
    def lower(self):
        """Derived lower bound ρ in the standardized z-space, chosen by the
        mean-reverting condition ζ(1, κ) = 0. Always in (-1, 0)."""
        return self._lower

    @property
    def leakage(self):
        # Two-sided: lower tail (1 - e^(-a)) plus upper tail e^(-b).
        a = self._lower + 1.0
        b = self._bounding + 1.0
        return (1.0 - math.exp(-a)) + math.exp(-b)

    @property
    def maxOrder(self):
        return self._maxOrder

    def __getitem__(self, n:int) -> float:
        if n < 0 or n >= self._maxOrder:
            return IndexError()
        return self._sMoment[n]


NORMAL = Normal(bounding=5.0)

UNIFORM = Uniform()

EXPONENTIAL = Exponential()
