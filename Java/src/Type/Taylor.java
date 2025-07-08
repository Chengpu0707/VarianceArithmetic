package Type;

import java.io.FileWriter;
import java.io.IOException;


public class Taylor {
    static final NormalMomentum IDEAL_MOMENTUM = new NormalMomentum(5.0);

    static final int MIN_MONOTONIC_COUNT = 20;

    static String INPUT_HEADER =
            "result\tvalue\tuncertainty\tinPrec\toutPrec\tbounding\tmaxOrder\tMinMonotonic" +
            "\tcheckMonotonic\tcheckStability\tcheckReliablity\tcheckPositive\tName";
    static String EXTENSION_HEADER = 
            "Order\tTaylor Value\tTaylor Uncertainty\tExponent\tMomentum\tMonotonics" +
            "\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty" +
            "\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty";
    static String OUTPUT_HEADER =
            "Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\tException";

    static void writeResult(final FileWriter fw, final String exception,
            final VarDbl value, final VarDbl variance) throws IOException {
        if (fw != null) {
            fw.write(OUTPUT_HEADER + "\n");
            fw.write(String.format("%.15e\t%.15e\t%.15e\t%.15e\t%s\n", 
                     value.value(), value.uncertainty(), variance.value(), variance.uncertainty(), 
                     exception));
            fw.flush();
        }

    }
    /*
    1d Taylor expansion for differential series {s1dTaylor} at {input}, with {name} for logging.
    When {inPrec} is true, calculate Taylor expnasion against the precision of {input}.
    When {outPrec} is true, the result of the Taylor expnasion is the precision.
    s1dTaylor[n] should already normalized by /n!. 
    The max order of expansion is {maxOrder}, which should not exceed momentum.Normal.MAX_ORDER

    When {checkMonotonic} is true, raise {NotMonotonicException} if 
        after full expansion, the monotonic count is still less than {MIN_MONOTONIC_COUNT}.
    It should always be True.

    When {checkStability} is true, raise {NotStableException} if 
        after full expansion, the value for the last expansion term is more than 
        LEAKAGE_5-fold of the expansion uncertainty.
    It should always be True.

    When {checkReliablity} is true, raise {NotReliableException} if
        the precision of the result variance is more than 1/5, in which 5 is the binding factor.
    It should always be True.

    When {checkPositive} is true, raise {NotPosive} if
        the expansion variance at any order becomes negative
    It should always be True.

    Both the result value and variance are guaranteed to be finite, otherwise
        raise {NotFiniteException} 

    Dump the expansion to {dumpPath} when it is provided.
    {dumpPath} can be read back and tested using verifyDumpFile()
     */
    static private VarDbl taylor1d(
            final VarDbl in, final String name, final UnionArray s1dTaylor, boolean inPrec, boolean outPrec,
            final String dumpPath, final Momentum momentum,
            boolean checkMonotonic, boolean checkStability, boolean checkPositive, boolean checkReliablity) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                InitException, IOException {
        final int length;
        final double val;
        if (s1dTaylor.sDbl != null) {
            length = s1dTaylor.sDbl.length;
            for (int i = 0; i < length; ++i) {
                if (!Double.isFinite(s1dTaylor.sDbl[i]))
                    throw new IllegalArgumentException(String.format("Invalid Taylor coefficient [%d/%d]=%s", 
                        i, length, s1dTaylor.sDbl[i]));
            }
            val = s1dTaylor.sDbl[0];
        }
        else if (s1dTaylor.sVar != null) {
            length = s1dTaylor.sVar.length;
            for (int i = 0; i < length; ++i) {
                if (s1dTaylor.sVar == null)
                    throw new IllegalArgumentException(String.format("Null Taylor coefficient [%d/%d]", 
                        i, length));
                if (!Double.isFinite(s1dTaylor.sVar[i].value()) || !Double.isFinite(s1dTaylor.sVar[i].variance())) {
                    throw new IllegalArgumentException(String.format("Invalid Taylor coefficient [%d/%d]=%s", 
                        i, length, s1dTaylor.sDbl[i]));
                }
            }
            val = s1dTaylor.sVar[0].value();
        }
        else
            throw new IllegalArgumentException("Invalid Taylor coefficient");
        if (in.uncertainty() == 0) {
            if (s1dTaylor.sDbl != null)
                return new VarDbl(val);
            else 
                return new VarDbl(s1dTaylor.sVar[0]);
        }

        FileWriter fw = (dumpPath == null)? null : new FileWriter(dumpPath);
        if (fw != null) {
            fw.write(INPUT_HEADER + "\n");
            fw.write(String.format("%.15e\t%.15e\t%.15e\t%b\t%b\t%.3f\t%d\t%d\t%b\t%b\t%b\t%b\t%s\n", 
                     (s1dTaylor.sDbl != null)? s1dTaylor.sDbl[0] : s1dTaylor.sVar[0].value(), 
                     in.value(), in.uncertainty(), inPrec, outPrec, 
                     momentum.bounding, momentum.maxOrder, MIN_MONOTONIC_COUNT,
                     checkMonotonic, checkStability, checkReliablity, checkPositive, name));
            fw.write(EXTENSION_HEADER + "\n");
        }

        int monotonics = 0;
        boolean monotonicPrev = true;

        VarDbl value = outPrec? new VarDbl(1) : new VarDbl((s1dTaylor.sDbl != null)? s1dTaylor.sDbl[0] : s1dTaylor.sVar[0].value());
        VarDbl variance = new VarDbl();
        final double unc = inPrec? in.uncertainty() /in.value() : in.uncertainty();
        double uncN = unc;
        VarDbl prevVariance = new VarDbl();
        VarDbl newValue = null, newVariance = null;
        int n = 1;
        for (; (n < momentum.maxOrder) && (n < length) && Double.isFinite(uncN) && (uncN > 0); ++n, uncN *= unc) {
            newValue = new VarDbl();
            newVariance = new VarDbl();
            String infinite = null;
            int j = 0;
            try {
                if (s1dTaylor.sDbl != null) {
                    if (uncN < 1)
                        newValue.addInPlace(s1dTaylor.sDbl[n] * (uncN * momentum.get(n)));
                    else
                        newValue.addInPlace(s1dTaylor.sDbl[n] * uncN * momentum.get(n));
                    double newVar = 0;
                    for (j = 1; j < n; ++j) {
                        if (uncN < 1)
                            newVar += s1dTaylor.sDbl[j] * s1dTaylor.sDbl[n - j] * (uncN *
                                (momentum.get(n) - momentum.get(j) * momentum.get(n - j)));
                        else
                            newVar += s1dTaylor.sDbl[j] * s1dTaylor.sDbl[n - j] * uncN *
                                (momentum.get(n) - momentum.get(j) * momentum.get(n - j));
                    }
                    newVariance.addInPlace(newVar);
                } else {
                    if (uncN < 1)
                        newValue.addInPlace(s1dTaylor.sVar[n].multiply(uncN * momentum.get(n)));
                    else
                        newValue.addInPlace(s1dTaylor.sVar[n].multiply(uncN).multiply(momentum.get(n)));
                    for (j = 1; j < n; ++j) {
                        if (uncN < 1)
                            newVariance.addInPlace(s1dTaylor.sVar[j].multiply(s1dTaylor.sVar[n - j]).multiply(uncN *
                                (momentum.get(n) - momentum.get(j) * momentum.get(n - j))));
                        else
                            newVariance.addInPlace(s1dTaylor.sVar[j].multiply(s1dTaylor.sVar[n - j]).multiply(uncN)
                                .multiply(momentum.get(n) - momentum.get(j) * momentum.get(n - j)));
                    }
                }
                value.addInPlace(newValue);
                variance.addInPlace(newVariance);
            } catch(InitException e) {
                infinite = "NotFinteException\t" + e.getMessage();
            } catch(Throwable e) {
                infinite = "NotFinteException\t" + e.getMessage();
            }
            if (! Double.isFinite(value.value()))
                infinite = "NotFinteException\tvalue";
            else if (! Double.isFinite(value.variance() + variance.value()))
                infinite = "NotFinteException\tvariance";
            else if (! Double.isFinite(variance.variance()))
                infinite = "NotFinteException\tVariance uncertainty";
            else if ((n & 1) == 0) {
                if (Math.abs(newVariance.value()) <= Math.abs(prevVariance.value()))
                    monotonics += 1;
                else if ((monotonics > MIN_MONOTONIC_COUNT) && monotonicPrev)
                    monotonicPrev = false;
                else
                    monotonics = 0;
                prevVariance = newVariance;
            }
            if (fw != null) {
                fw.write(String.format("%d\t%.15e\t%.15e\t%.15e\t%.15e\t%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", 
                        n, 
                        (s1dTaylor.sDbl != null)? s1dTaylor.sDbl[n] : s1dTaylor.sVar[n].value(),
                        (s1dTaylor.sDbl != null)? 0 : s1dTaylor.sVar[n].uncertainty(),
                        uncN, momentum.get(n), monotonics, 
                        value.value(), value.uncertainty(), variance.value(), variance.uncertainty(),
                        newValue.value(), newValue.variance(), newVariance.value(), newVariance.uncertainty()));
                fw.flush();
            }
            if (infinite != null) {
                writeResult(fw, infinite, value, variance);
                throw new NotFiniteException("UnionArray: " + infinite,
                        name, s1dTaylor, inPrec, outPrec,
                        in, value, variance, n, newValue, newVariance, monotonics);
            }
            if (checkPositive && (variance.value() < 0)) {
                writeResult(fw, "NotPositiveException", value, variance);
                throw new NotPositiveException("negative variance",
                        name, s1dTaylor, inPrec, outPrec,
                        in, value, variance, n, newValue, newVariance, monotonics);
            }
        }

        if (checkMonotonic && (uncN > 0) && (monotonics < MIN_MONOTONIC_COUNT)) {
            writeResult(fw, "NotMonotonicException", value, variance);
            throw new NotMonotonicException(String.format("Taylor1d: MIN_MONOTONIC_COUNT=%d", MIN_MONOTONIC_COUNT),
                    name, s1dTaylor, inPrec, outPrec,
                    in, value, variance, n, newValue, newVariance, monotonics);
        }
        if (checkStability) {
            if (newValue.value() > Math.sqrt(variance.value()) * (1 - momentum.get(0))) {
                writeResult(fw, "NotStableException\tvariance", value, variance);
                throw new NotStableException(String.format("Taylor1d: limit=%e", unc),
                        name, s1dTaylor, inPrec, outPrec,
                        in, value, variance, n, newValue, newVariance, monotonics);
            }
            if (Math.abs(newValue.value()) > Math.ulp(value.value())) {
                writeResult(fw, "NotStableException\tvalue", value, variance);
                throw new NotStableException(String.format("Taylor1d: limit=%e", unc),
                        name, s1dTaylor, inPrec, outPrec,
                        in, value, variance, n, newValue, newVariance, monotonics);
            }
        } 

        if (outPrec) {
            value.multiplyInPlace(val);
            variance.multiplyInPlace(val * val);
        }
        if (! Double.isFinite(value.value())) {
            writeResult(fw, "NotFinteException\tvalue", value, variance);
            throw new NotFiniteException(String.format("Taylor1d: value overflow %e", val),
                    name, s1dTaylor, inPrec, outPrec,
                    in, value, variance, n, newValue, newVariance, monotonics);
        }
        if (! Double.isFinite(variance.value() + value.variance())) {
            writeResult(fw, "NotFinteException\tvariance", value, variance);
            throw new NotFiniteException(String.format("Taylor1d: variance overflow %e",val),
                    name, s1dTaylor, inPrec, outPrec,
                    in, value, variance, n, newValue, newVariance, monotonics);
        }
        writeResult(fw, "", value, variance);
        return new VarDbl(value.value(), Math.sqrt(variance.value() + value.variance()));
    } 

    static public VarDbl taylor1d(final VarDbl in, final String name, VarDbl[] s1dTaylor, boolean inPrec, boolean outPrec, 
                                  final String dumpPath, final Momentum momment) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                InitException, IOException {
        return taylor1d(in, name, new UnionArray(s1dTaylor), inPrec, outPrec, 
            dumpPath, momment, 
            true, true, true, true);
    }
    static public VarDbl taylor1d(final VarDbl in, final String name, VarDbl[] s1dTaylor, boolean inPrec, boolean outPrec, 
                                  final Momentum momment) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                InitException, IOException {
        return taylor1d(in, name, new UnionArray(s1dTaylor), inPrec, outPrec, 
            null, momment, 
            true, true, true, true);
    }
    static public VarDbl taylor1d(final VarDbl in, final String name, VarDbl[] s1dTaylor, boolean inPrec, boolean outPrec, 
                                  final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                InitException, IOException {
        return taylor1d(in, name, new UnionArray(s1dTaylor), inPrec, outPrec, 
            dumpPath, IDEAL_MOMENTUM, 
            true, true, true, true);
    }
    static public VarDbl taylor1d(final VarDbl in, final String name, VarDbl[] s1dTaylor, boolean inPrec, boolean outPrec) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException {
        try {
            return taylor1d(in, name, s1dTaylor, inPrec, outPrec, IDEAL_MOMENTUM);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    static public VarDbl taylor1d(final VarDbl in, final String name, double[] s1dTaylor, boolean inPrec, boolean outPrec, 
                                  final String dumpPath, final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return taylor1d(in, name, new UnionArray(s1dTaylor), inPrec, outPrec, 
            dumpPath, momentum, 
            true, true, true, true);
    }
    static public VarDbl taylor1d(final VarDbl in, final String name, double[] s1dTaylor, boolean inPrec, boolean outPrec, 
                                  final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return taylor1d(in, name, new UnionArray(s1dTaylor), inPrec, outPrec, 
            dumpPath, IDEAL_MOMENTUM, 
            true, true, true, true);
    }
    static public VarDbl taylor1d(final VarDbl in, final String name, double[] s1dTaylor, boolean inPrec, boolean outPrec, 
                                  final Momentum momment) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return taylor1d(in, name, new UnionArray(s1dTaylor), inPrec, outPrec, 
            null, momment, 
            true, true, true, true);
    }
    static public VarDbl taylor1d(final VarDbl in, final String name, double[] s1dTaylor, boolean inPrec, boolean outPrec) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException  {
        return taylor1d(in, name, s1dTaylor, inPrec, outPrec, IDEAL_MOMENTUM);
    }

    /*
    1d Taylor expansion for poly1d at "input" with "sCoeff".
    Allow input.value() +- input.uncertainty() to include 0
    */
    static private VarDbl poly1d(final VarDbl in, final UnionArray s1dCoeff, 
                                 final String dumpPath, final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        final int length;
        final UnionArray s1dTaylor;
        if (s1dCoeff.sDbl != null) {
            length = s1dCoeff.sDbl.length;
            s1dTaylor = new UnionArray(new double[length * 2 - 1]);
        }
        else if (s1dCoeff.sVar != null) {
            length = s1dCoeff.sVar.length;
            s1dTaylor = new UnionArray(new VarDbl[length * 2 - 1]);
        }
        else
            throw new IllegalArgumentException("Invalid Taylor coefficient");
        if (length * 2 >= momentum.maxOrder) {
            throw new IllegalArgumentException(String.format("coefficient length %d >= %d", 
                    length, momentum.maxOrder / 2));
        }
        final int exp = length - 1;
        if (s1dCoeff.sVar != null) {
            s1dTaylor.sVar[0] = s1dCoeff.sVar[0];
            for (int k = 1; k < s1dTaylor.sVar.length; ++k)
                s1dTaylor.sVar[k] = new VarDbl();
        } else {
            s1dTaylor.sDbl[0] = s1dCoeff.sDbl[0];
        }
        
        final double[] sPow = new double[exp + 1];
        sPow[0] = 1;
        if (exp >= 1)
            sPow[1] = in.value();
        final int[] sTaylor = new int[2*exp + 1];
        sTaylor[0] = 1;

        for (int j = 1; j < length; ++j) {
            VarDbl coeff = null;
            if (s1dCoeff.sVar != null) {
                coeff = s1dCoeff.sVar[j];
                if ((coeff == null)) {
                    throw new IllegalArgumentException(String.format("null coefficient at index %d/%d is too large", 
                            j, length));
                }
                if (!Double.isFinite(coeff.value()) || !Double.isFinite(coeff.variance())) {
                    throw new InitException(String.format("poly1d coefficent at %d", j), 
                            coeff.value(), coeff.variance());
                }
            }
            sTaylor[1] = j;
            for (int k = 2; k <= j; ++k) {
                sTaylor[k] = sTaylor[k - 1] * (j + 1 - k) / k;
                if (sPow[k] == 0)
                    sPow[k] = sPow[k - 1] * in.value();
            }
            for (int k = 0; k <= j; ++k) {
                if (coeff != null)
                    s1dTaylor.sVar[k].addInPlace(coeff.multiply(sTaylor[k] * sPow[j - k]));
                else
                    s1dTaylor.sDbl[k] += s1dCoeff.sDbl[j] * sTaylor[k] * sPow[j - k];
            }
        }
        return taylor1d(in, String.format("Poly[%d]", exp), 
                    s1dTaylor, false, false, dumpPath, momentum,
                    false, false, true, true);
    }

    static public VarDbl poly1d(final VarDbl in, final VarDbl[] sCoeff, 
                                final String dumpPath, final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return poly1d(in, new UnionArray(sCoeff), dumpPath, momentum);
    }
    static public VarDbl poly1d(final VarDbl in, final VarDbl[] sCoeff, 
                                final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return poly1d(in, new UnionArray(sCoeff), dumpPath, IDEAL_MOMENTUM);
    }
    static public VarDbl poly1d(final VarDbl in, final VarDbl[] sCoeff, 
                                final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return poly1d(in, new UnionArray(sCoeff), null, momentum);
    }
    static public VarDbl poly1d(final VarDbl in, final VarDbl[] sCoeff) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return poly1d(in, sCoeff, IDEAL_MOMENTUM);
    }

    static public VarDbl poly1d(final VarDbl in, final double[] sCoeff, 
                                final String dumpPath, final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return poly1d(in, new UnionArray(sCoeff), dumpPath, momentum);
    }
    static public VarDbl poly1d(final VarDbl in, final double[] sCoeff, 
                                final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return poly1d(in, new UnionArray(sCoeff), dumpPath, IDEAL_MOMENTUM);
    }
    static public VarDbl poly1d(final VarDbl in, final double[] sCoeff, 
                                final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return poly1d(in, new UnionArray(sCoeff), null, momentum);
    }
    static public VarDbl poly1d(final VarDbl in, final double[] sCoeff) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return poly1d(in, sCoeff, IDEAL_MOMENTUM);
    }


    /*
     * @return in^exponenet
     */
    static public VarDbl pow(final VarDbl in, final double exponent, 
                             final String dumpPath, final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        if (exponent == 0) {
            return new VarDbl(1, 0);
        }
        if (exponent == 1) {
            return new VarDbl(in);
        }
        if ((0 < exponent) && (Math.floor(exponent) == Math.ceil(exponent))) {
            final int exp = (int) exponent;
            final double[] sCoeff = new double[exp + 1];
            for (int i = 0; i < exp; ++i) {
                sCoeff[i] = 0;
            }
            sCoeff[exp] = 1;
            return poly1d(in, sCoeff, dumpPath, momentum);
        }
        if (in.value() < 0) {
            throw new IllegalArgumentException(String.format("pow(%s, %s): negative base", in, exponent));
        }
        final VarDbl[] sTaylor = new VarDbl[momentum.maxOrder];
        sTaylor[0] = new VarDbl(Math.pow(in.value(), exponent));
        sTaylor[1] = new VarDbl(exponent);
        double exp = exponent - 1;
        for (int i = 2; i < momentum.maxOrder; ++i, --exp) {
            sTaylor[i] = sTaylor[i - 1].multiply(exp / i);
        }
        return taylor1d(in, String.format("(%s)^%s", in, exponent), 
                        sTaylor, true, true, dumpPath, momentum);
    }
    static public VarDbl pow(final VarDbl in, double exponent, 
                             final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return pow(in, exponent, dumpPath, IDEAL_MOMENTUM);
    }
    static public VarDbl pow(final VarDbl in, double exponent, 
                             final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return pow(in, exponent, null, momentum);
    }
    static public VarDbl pow(final VarDbl in, double exponent) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return pow(in, exponent, IDEAL_MOMENTUM);
    }

    /*
     * @return sin(x)
     */
    static public VarDbl sin(final VarDbl in, 
                             final String dumpPath, final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        final double[] sTaylor = new double[momentum.maxOrder];
        final double x = in.value();
        sTaylor[0] = Math.sin(x);
        double fac = 1;
        for (int i = 1; i < momentum.maxOrder; ++i, fac /= i) {
            switch (i%4) {
            case 0:
                sTaylor[i] = Math.sin(x) * fac; 
                break;
            case 1:
                sTaylor[i] = Math.cos(x) * fac;   
                break;
            case 2:
                sTaylor[i] = -Math.sin(x) * fac; 
                break;
            case 3:
                sTaylor[i] = -Math.cos(x) * fac; 
                break;
            }
        }
        return taylor1d(in, String.format("sin(%s)", in), sTaylor, false, false, dumpPath, momentum);
    }
    static public VarDbl sin(final VarDbl in, 
                             final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return sin(in, dumpPath, IDEAL_MOMENTUM);
    }
    static public VarDbl sin(final VarDbl in, 
                             final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return sin(in, null, momentum);
    }
    static public VarDbl sin(final VarDbl in) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return sin(in, null, IDEAL_MOMENTUM);
    }

    /*
     * @return exp(x)
     */
    static public VarDbl exp(final VarDbl in, 
                             final String dumpPath, final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        final double[] sTaylor = new double[momentum.maxOrder];
        sTaylor[0] = Math.exp(in.value());
        double fac = 1;
        for (int i = 1; i < momentum.maxOrder; ++i, fac /= i) {
            sTaylor[i] = 1.0 *fac;
        }
        return taylor1d(in, String.format("exp(%s)", in), sTaylor, false, true, 
                        dumpPath, momentum);
    }
    static public VarDbl exp(final VarDbl in, 
                             final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return exp(in, dumpPath, IDEAL_MOMENTUM);
    }
    static public VarDbl exp(final VarDbl in, 
                             final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return exp(in, null, momentum);
    }
    static public VarDbl exp(final VarDbl in) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return exp(in, null, IDEAL_MOMENTUM);
    }

    /*
     * @return log(x)
     */
    static public VarDbl log(final VarDbl in,
                             final String dumpPath, final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        final double[] sTaylor = new double[momentum.maxOrder];
        sTaylor[0] = Math.log(in.value());
        for (int i = 1; i < momentum.maxOrder; ++i) {
            sTaylor[i] =((i%2) == 1)? +1.0/i : -1.0/i;
        }
        return taylor1d(in, String.format("log(%s)", in), sTaylor, true, false, dumpPath, momentum);
    }
    static public VarDbl log(final VarDbl in,
                             final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
            return log(in, dumpPath, IDEAL_MOMENTUM);
    }
    static public VarDbl log(final VarDbl in,
                             final Momentum momentum) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
            return log(in, null, momentum);
    }
    static public VarDbl log(final VarDbl in) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return log(in, null, IDEAL_MOMENTUM);
    }


}
