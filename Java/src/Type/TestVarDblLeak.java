package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

public class TestVarDblLeak extends VarDbl {
    @Test
    public void test0() {
        try {
            assertEquals( 0L, VarDbl.leak(0) );
            assertEquals( 0L, VarDbl.leak(VarDbl.LEAK_FACTOR * 0.5) );
            pack(0, 0, false, false, 0);
            assertEquals( 0, leakage(), 0 );
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testMin() {
        try {
            assertEquals( 1L, VarDbl.leak(VarDbl.LEAK_FACTOR * 0.500001) );
            assertEquals( 1L, VarDbl.leak(VarDbl.LEAK_FACTOR) );
            pack(0, 0, false, false, 1);
            assertEquals( VarDbl.LEAK_FACTOR, leakage(), 0 );
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testIntermedia() {
        try {
            final long leak = (6 << LEAK_EXP_SHIFT) | 8;
            assertEquals( leak, VarDbl.leak(1.0/1024) );
            assertEquals( leak, VarDbl.leak(0.001) );
            pack(0, 0, false, false, leak);
            assertEquals( 1.0/1024, leakage(), 0 );
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testMax() {
        try {
            final long maxLeak = (VarDbl.LEAK_EXP_MAX << LEAK_EXP_SHIFT) | VarDbl.LEAK_VAL_MAX;
            assertEquals( maxLeak, VarDbl.leak(VarDbl.LEAK_MAX) );
            assertEquals( maxLeak, VarDbl.leak( 1.0 ) );
            pack(0, 0, false, false, maxLeak);
            assertEquals( VarDbl.LEAK_MAX, leakage(), 0 );
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testAdd() {
        try {
            final long leak = (6 << LEAK_EXP_SHIFT) | 8;
            pack(0, 0, false, false, leak);
            assertEquals( 1.0/1024, leakage(), 0 );
            VarDbl op = new VarDbl();
            op.pack(0, 0, false, false, leak);
            assertEquals( 1.0/1024, op.leakage(), 0 );
            this.add(op);
            assertEquals( 2.0/1024, leakage(), 0 );
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testMultiple() {
        try {
            final long leak = (6 << LEAK_EXP_SHIFT) | 8;
            pack(0, 0, false, false, leak);
            assertEquals( 1.0/1024, leakage(), 0 );
            VarDbl op = new VarDbl();
            op.pack(0, 0, false, false, leak);
            assertEquals( 1.0/1024, op.leakage(), 0 );
            this.multiply(op);
            assertEquals( 2.0/1024, leakage(), 0 );
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

}
