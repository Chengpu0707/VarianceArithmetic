package Type;

import static org.junit.Assert.fail;
import org.junit.Test;

import java.io.FileWriter;
import java.io.IOException;


public class TestConvergence {
    static final String EDGE_HEADER = "Edge Value\tEdge Uncertainty\tValue\tUncertainty\tException\n";

    @Test
    public void testPower() {
        try (FileWriter f = new FileWriter("./Java/Output/PowerAtOneEdge.txt")) {
            f.write(TestConvergence.EDGE_HEADER);
            for (int i = -60; i <= 60; i += 1) {
                final double exp = i/10.;
                String except = "";
                for (int j = 300; j >= 100; --j) {
                    final double dx = j/1000.;
                    try {
                        final VarDbl res = new VarDbl(1, dx).power(exp);
                        f.write(String.format("%e\t%e\t%e\t%e\t%s\n", 
                                exp, dx, res.value(), res.uncertainty(), except));
                        break;
                    } catch (InitException | DivergentException | NotReliableException | NotMonotonicException
                            | NotStableException e) {
                        except = e.getClass().getName();
                    }
                }
            }
        } catch (IOException ex) {
            fail(ex.getMessage());
        }
    }

    @Test
    public void testSin() {
        try (FileWriter f = new FileWriter("./Java/Output/SinEdge.txt")) {
            f.write(TestConvergence.EDGE_HEADER);
            for (int i = -16; i <= 16; i += 1) {
                final double x = i/16. * Math.PI;
                String except = "";
                for (int j = 150; j >= 50; --j) {
                    final double dx = j/100.;
                    try {
                        final VarDbl res = new VarDbl(x, dx).sin();
                        f.write(String.format("%e\t%e\t%e\t%e\t%s\n", 
                                x, dx, res.value(), res.uncertainty(), except));
                        break;
                    } catch (InitException | DivergentException | NotReliableException | NotMonotonicException
                            | NotStableException e) {
                        except = e.getClass().getName();
                    }
                }
            }
        } catch (IOException ex) {
            fail(ex.getMessage());
        }
    }

    @Test
    public void testLog() {
        try (FileWriter f = new FileWriter("./Java/Output/LogEdge.txt")) {
            f.write(TestConvergence.EDGE_HEADER);
            for (int i = -15; i <= 16; i += 1) {
                final double x = i/16. + 1;
                String except = "";
                for (int j = 250; j >= 0; --j) {
                    final double dx = j/100.;
                    try {
                        final VarDbl res = new VarDbl(x, dx).log();
                        f.write(String.format("%e\t%e\t%e\t%e\t%s\n", 
                                x, dx, res.value(), res.uncertainty(), except));
                        break;
                    } catch (InitException | DivergentException | NotReliableException | NotMonotonicException
                            | NotStableException e) {
                        except = e.getClass().getName();
                    }
                }
            }
        } catch (IOException ex) {
            fail(ex.getMessage());
        }
    }
}
