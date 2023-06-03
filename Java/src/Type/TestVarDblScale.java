package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

public class TestVarDblScale extends VarDbl {
    
    @Test
    public void testNegate() {
        try {
            init(1.01, 0.01, false);
            final VarDbl other = (VarDbl) this.negate();
            assertEquals(exp(), other.exp());
            assertEquals(!neg(), other.neg());
            assertEquals(val(), other.val());
            assertEquals(var(), other.var());
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testShift() {
        VarDbl other;
        try {
            init(1, 1.0 / 1024 / 1024, false);
            assertEquals(-36, exp());
            assertEquals(1L << 36, val());
            assertEquals(Dbl.DOUBLE_VAL_EXTRA, var());
            
            other = (VarDbl) shift(0);
            assertEquals(-36, other.exp());
            assertEquals(1L << 36, other.val());
            assertEquals(Dbl.DOUBLE_VAL_EXTRA, other.var());

            other = (VarDbl) shift(4);
            assertEquals(-32, other.exp());
            assertEquals(1L << 36, other.val());
            assertEquals(Dbl.DOUBLE_VAL_EXTRA, other.var());

            other = (VarDbl) shift(-10);
            assertEquals(-46, other.exp());
            assertEquals(1L << 36, other.val());
            assertEquals(Dbl.DOUBLE_VAL_EXTRA, other.var());

        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testShifAtMax() {
        VarDbl other;
        try {
            init(1, 1.0 / 1024 / 1024, false);
            assertEquals(-36, exp());
            assertEquals(1L << 36, val());
            assertEquals(Dbl.DOUBLE_VAL_EXTRA, var());

            other = (VarDbl) shift(Dbl.DOUBLE_EXP_MAX + 36);
            assertEquals(Dbl.DOUBLE_EXP_MAX, other.exp());
            assertEquals(1L << 36, other.val());
            assertEquals(Dbl.DOUBLE_VAL_EXTRA, other.var());
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            other = (VarDbl) shift(Dbl.DOUBLE_EXP_MAX + 37);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail();
        }
     }

    @Test
    public void testShifAtMin() {
        VarDbl other;
        try {
            init(1, 1.0 / 1024 / 1024, false);
            assertEquals(-36, exp());
            assertEquals(1L << 36, val());
            assertEquals(Dbl.DOUBLE_VAL_EXTRA, var());

            other = (VarDbl) shift(Dbl.DOUBLE_EXP_MIN + 36);
            assertEquals(Dbl.DOUBLE_EXP_MIN, other.exp());
            assertEquals(1L << 36, other.val());
            assertEquals(Dbl.DOUBLE_VAL_EXTRA, other.var());

            other = (VarDbl) shift(Dbl.DOUBLE_EXP_MIN + 35);
            assertEquals(Dbl.DOUBLE_EXP_MIN, other.exp());
            assertEquals(1L << 35, other.val());
            assertEquals(Dbl.DOUBLE_VAL_EXTRA >>> 2, other.var());

            other = (VarDbl) shift(Dbl.DOUBLE_EXP_MIN + 11);
            assertEquals(Dbl.DOUBLE_EXP_MIN, other.exp());
            assertEquals(1 << 11, other.val());
            assertEquals(4, other.var());

            other = (VarDbl) shift(Dbl.DOUBLE_EXP_MIN + 10);
            assertEquals(Dbl.DOUBLE_EXP_MIN, other.exp());
            assertEquals(1 << 10, other.val());
            assertEquals(1, other.var());

            // the uncertainty can not be getten rid of by shift
            other = (VarDbl) shift(Dbl.DOUBLE_EXP_MIN + 9);
            assertEquals(Dbl.DOUBLE_EXP_MIN, other.exp());
            assertEquals(1 << 9, other.val());
            assertEquals(0, other.var());

            // when the value becomes precicse
            other = (VarDbl) shift(-2000);
            assertEquals(Dbl.DOUBLE_EXP_MIN, other.exp());
            assertEquals(0, other.val());
            assertEquals(0, other.var());
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testScale() {
        VarDbl other;
        try {
            init(Math.sqrt(2), 2.0 / 1024 / 1024, false);  // 10 + 15 bits, 

            other = (VarDbl) scale(-1);
            assertEquals(-1.0, other.value() / Math.sqrt(2), 2E-8);
            assertEquals( 1.0, other.uncertainty() / (Math.sqrt(2) / 1024), 1E-16);

            other = (VarDbl) scale(Math.sqrt(2));
            assertEquals(1.0, other.value() / 2.0, 3E-8);
            assertEquals(1.0, other.uncertainty() / (2.0 / 1024), 1E-16);

            other = (VarDbl) scale(1.0 / Math.sqrt(2));
            assertEquals(1.0, other.value() / 1.0, 3E-8);
            assertEquals(1.0, other.uncertainty() / (1.0 / 1024), 2E-16);

            other = (VarDbl) scale(0);
            assertEquals(0.0, other.value(), 1E-16);
            assertEquals(0.0, other.uncertainty(), 1E-16);
            
            other = (VarDbl) scale(1.0 / (1 << 15));
            assertEquals(1.0, other.value() / (Math.sqrt(2) / (1 << 15)), 2E-8);
            assertEquals(1.0, other.uncertainty() /(Math.sqrt(2) / (1 << 25)), 1E-16);
            
            other = (VarDbl) scale(Double.MIN_VALUE);
            assertEquals(1.0, other.value() / (Math.sqrt(2) * Double.MIN_VALUE), 2E-8);
            assertEquals(0.0, other.uncertainty(), 1E-16);
            
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            scale(Double.MAX_VALUE / 2);
        } catch (ValueException e) {
            fail();;
        } catch (UncertaintyException e) {
        }
    }
}
