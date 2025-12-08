package Func;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import Type.InitException;
import Type.VarDbl;


enum SinSource {
    Prec,
    Quart,
    Lib
}


/*
 * A cached sine function with resolution down to PI /2^order().
 * It has integer input to sin(), cos() and tan() to avoid rounding errors.
 */
public class IndexSin {

    /*
     * The maximual order of FFT calculation.
     * The size of FFT is limited to MAX_QUARD_ORDER
     */

    static public final int MIN_ORDER = 1;
    static public final int MAX_ORDER = 18;

    private static final VarDbl[] _sPrecSin = new VarDbl[(1 << (MAX_ORDER - 1)) + 1];

    private static final VarDbl[] _sQuartSin = new VarDbl[(1 << (MAX_ORDER - 1)) + 1];
    static {
        final int size = 1 << MAX_ORDER;
        final int half = size >> 1;
        final int quart = half >> 1;
        _sQuartSin[0] = new VarDbl(0);
        _sQuartSin[half] = new VarDbl(1);
        try {
            _sQuartSin[quart] = new VarDbl(Math.sin(Math.PI/4));
        } catch(InitException ex) {
            throw new RuntimeException(
                String.format("The index %d/%d for IndexSin has exception %s", size, MAX_ORDER, ex));
        }
        final double f = Math.PI/size;
        for (int i = 1; i < quart; ++i) {
            try {
                _sQuartSin[i] = new VarDbl(Math.sin(f*i));
                _sQuartSin[half - i] = new VarDbl(Math.cos(f*i));
            } catch(InitException ex) {
                throw new RuntimeException(
                    String.format("The index %d/%d for IndexSin has exception %s", i, MAX_ORDER, ex));
            }
        }
    }


    public final SinSource sinSource;
 
    public IndexSin(SinSource sinSource) {
        if (sinSource == SinSource.Prec) {
            if (_sPrecSin[0] == null)
                readPrec();
        } else if (sinSource != SinSource.Quart && sinSource != SinSource.Lib) {
            throw new IllegalArgumentException(
                String.format("Unknown SinSource %s for IndexSin", sinSource));
        }
        this.sinSource = sinSource;
    }

    void readPrec() {
        final String dumpPath = "./Cpp/Output/IndexSin_Prec_18.txt";
        try (final BufferedReader br = new BufferedReader(new FileReader(dumpPath))) {
            String line = br.readLine();
            if (! line.equals("Index\tX\tValue\tUncertainty"))
                throw new IllegalArgumentException(String.format("Invalid header in file %s: %s", dumpPath, line));
            for (int i = 0; i < _sPrecSin.length; ++i) {
                line = br.readLine();
                if (line == null)
                    throw new IllegalArgumentException(String.format("Invalid %d null line in file %s", i, dumpPath));
                final String[] sWord = line.trim().split("\t");
                if (sWord.length != 4)
                    throw new IllegalArgumentException(String.format("Invalid %d line with %d colums in file %s: %s", i, sWord.length, dumpPath, line));
                try {
                    final int n = Integer.parseInt(sWord[0]);
                    if (n != i)
                        throw new IllegalArgumentException(String.format("Invalid %d line index %d in file %s: %s", i, n, dumpPath, line));
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException(String.format("Invalid %d line at index count in file %s: %s", i, dumpPath, line));
                }
                try {
                    final double x = Double.parseDouble(sWord[1]);
                    final double rad = ((double) i) / (1 << MAX_ORDER);
                    if (Math.abs(x - rad) > Math.ulp(1))
                        throw new IllegalArgumentException(String.format("Invalid %d line x %e vs %e in file %s: %s", i, x, rad, dumpPath, line));
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException(String.format("Invalid %d line at index count in file %s: %s", i, dumpPath, line));
                }
                final double val, unc;
                try {
                    val = Double.parseDouble(sWord[2]);
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException(String.format("Invalid %d line at value '%s' in file %s: %s", i, sWord[2], dumpPath, line));
                }
                try {
                    unc = Double.parseDouble(sWord[3]);
                    if ((unc < 0) || (unc > 2.*Math.ulp(val)))
                       throw new IllegalArgumentException(String.format("Invalid %d line uncertainty %e in file %s: %s", i, unc, dumpPath, line));
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException(String.format("Invalid %d line at value '%s' in file %s: %s", i, sWord[2], dumpPath, line));
                }
                try {
                    _sPrecSin[i] = new VarDbl(val, unc);
                } catch (InitException e) {
                    throw new IllegalArgumentException(String.format("Invalid %d line of VarDbl(%e, %e) in file %s: %s", i, val, unc, dumpPath, line));
                }
            }
        } catch (IOException e) {   
            throw new IllegalArgumentException(String.format("Invalid Prec file '%s': %s", dumpPath, e.getMessage()));
        } 
    }

    int get_index(long freq, int order) {
        if (order < MIN_ORDER) {
             throw new IllegalArgumentException(
                String.format("The order %d < %d for IndexSin", order, MIN_ORDER));
        }
        if (MAX_ORDER < order) {
            throw new IllegalArgumentException(
                String.format("The order %d > %d for IndexSin", order, MAX_ORDER));
        }
        final int half = 1 << (order - 1);
        boolean neg = freq < 0;
        long div = Math.abs(freq) / half;
        long rem = Math.abs(freq) % half;
        if ((div & 1) == 1) {
            rem = half - rem;
        }
        if ((div & 2) == 2) {
            neg = !neg;
        }
        return (int) (neg? -rem : rem);
    }
    
    public VarDbl sin(long freq, int order) {
        if (sinSource == SinSource.Quart) {
            final int idx = this.get_index(freq, order);
            final VarDbl ret = new VarDbl(_sQuartSin[Math.abs(idx) << (MAX_ORDER - order)]);
            return (idx >= 0)? ret : ret.negate();
        } else if (sinSource == SinSource.Prec) {
            final int idx = this.get_index(freq, order);
            final VarDbl ret = new VarDbl(_sPrecSin[Math.abs(idx) << (MAX_ORDER - order)]);
            return (idx >= 0)? ret : ret.negate();
        } else if (sinSource == SinSource.Lib) {
            try {
                double ret = Math.sin(Math.PI / (1 << order) * freq);
                return new VarDbl(ret);
            } catch (InitException e) {
                throw new RuntimeException(String.format(
                    "The freq %d/%d for IndexSin.sin()", freq, order));
            }
        } else {
            throw new IllegalArgumentException(
                String.format("Unknown SinSource %s for IndexSin.sin()", sinSource));
        }
    }

    public VarDbl cos(long freq, int order) {
        return sin(freq + (1 << (order - 1)), order);
    }

    public boolean dump(final String dumpPath, int order) {
        try (FileWriter writer = new FileWriter(dumpPath)) {
            writer.write("Index\tX\tValue\tUncertainty");
            final double f = 1.0 / (1 << order);
            for (int i = 0; i < _sQuartSin.length; ++i) {
                final VarDbl sin = sin(i, order);
                writer.write(String.format("\n%d\t%.20e\t%.20e\t%.20e", 
                             i, i *f, sin.value(), sin.uncertainty()));
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }            
        return true;
    }
}
