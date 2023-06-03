package Stats;

import static org.junit.Assert.fail;

import java.io.FileWriter;
import java.io.IOException;

import org.junit.Test;

public class TestPower {
    /*
     * To see if the x^p has a minimal at ether p=11/7 or p=5/3 using Mont Calo method
     */

    static final double[] sPower = new double[] {
        -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
        2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
        3.0
    };
    static final double[] sPrec = new double[] {
        1E-1, 1E-2, 1E-3, 1E-4, 1E-5, 1E-6
    };
    static final int count = 1000000;

    static final Stat[] sStat = new Stat[sPower.length];
    static {
        for (int i = 0; i < sPower.length; ++i) {
            sStat[i] = new Stat();
        }    
    }

    static final String pathToBias0 = "C:\\Users\\Cheng\\Documents\\Proj\\VarianceArithemtic\\Java\\Output\\bias.0.csv";
    static final String pathToCost1 = "C:\\Users\\Cheng\\Documents\\Proj\\VarianceArithemtic\\Java\\Output\\var.csv";
    static final String pathToCost0 = "C:\\Users\\Cheng\\Documents\\Proj\\VarianceArithemtic\\Java\\Output\\var.0.csv";
    static final String pathToBias1 = "C:\\Users\\Cheng\\Documents\\Proj\\VarianceArithemtic\\Java\\Output\\var.1.csv";


    private void writeHeader(final String path) {
        final int start = path.lastIndexOf('\\');
        final int end = path.lastIndexOf(".csv");
        if ((start < 0) || (end < 0)) {
            fail();
        }
        final String title = path.substring(start + 1, end);
        try {
            FileWriter writer = new FileWriter(path);
            writer.write(title);
            for (int i = 0; i < sPower.length; ++i) {
                writer.write(String.format(",%.3e", sPower[i]));
            }
            writer.write("\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }

    private void writeRow(final String path, double prec) {
        final double prec2 = prec * prec;
        try {
            FileWriter writer = new FileWriter(path, true);
            writer.write(String.format("%.3e", prec));
            for (int i = 0; i < sPower.length; ++i) {
                final double pw = sPower[i];
                final double value;
                switch(path) {
                case pathToBias0:
                    value = (sStat[i].avg() - 1) / prec2;
                    break;
                case pathToBias1:
                    value = (sStat[i].avg() - 1 - pw * (pw - 1) * prec2 / 2.0) / (prec2 * prec2);
                    break;
                case pathToCost0:
                    value = sStat[i].var() / prec2;
                    break;
                case pathToCost1:
                    value = (sStat[i].var() - pw * pw *  prec2) / prec2;
                    break;
                default:
                    value = 0;
                    fail();
                }
                writer.write(String.format(",%.3e", value));
            }
            writer.write("\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }

    }

    @Test
    public void test() {
        writeHeader(pathToBias0);
        writeHeader(pathToBias1);
        writeHeader(pathToCost0);
        writeHeader(pathToCost1);

        final Noise noise = new Noise();

        for (double prec : sPrec) {
            for (int i = 0; i < sPower.length; ++i) {
                sStat[i].clear();
            }
            for (int j = 0; j < count; ++j) {
                final double op = 1 + noise.gaussian(prec);
                for (int i = 0; i < sPower.length; ++i) {
                    final double res = Math.pow(op, sPower[i]);
                    if (!sStat[i].accum(res)) {
                        fail();
                    }
                }
            }
            writeRow(pathToBias0, prec);
            writeRow(pathToBias1, prec);
            writeRow(pathToCost0, prec);
            writeRow(pathToCost1, prec);
        }
    }
}
