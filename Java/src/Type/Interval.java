package Type;

/**
 * Best-effort double-precision interval arithmetic with ULP outward padding
 * after each operation. Not IEEE-1788 rigorous (no directed-rounding mode),
 * but tight enough for diagnostic comparison against variance arithmetic.
 *
 * Conversion from VarDbl: [value - max(unc, ulp(value)), value + max(unc, ulp(value))],
 * so degenerate (unc=0) sources like IndexSin Quart still get a proper ULP-wide interval.
 */
public final class Interval {

    private final double lo;
    private final double hi;

    public Interval()                            { this.lo = 0.0; this.hi = 0.0; }
    public Interval(double point)                { this.lo = point; this.hi = point; }
    public Interval(double lo, double hi) {
        this.lo = Math.min(lo, hi);
        this.hi = Math.max(lo, hi);
    }
    public Interval(VarDbl v) {
        final double w = Math.max(v.uncertainty(), safeUlp(v.value()));
        this.lo = v.value() - w;
        this.hi = v.value() + w;
    }
    public static Interval centered(double value, double unc) {
        if (unc < 0)
            throw new IllegalArgumentException("Interval.centered: unc must be >= 0");
        final double w = Math.max(unc, safeUlp(value));
        return new Interval(value - w, value + w);
    }

    public double lo()  { return lo; }
    public double hi()  { return hi; }
    public double mid() { return 0.5 * (lo + hi); }
    public double rad() { return 0.5 * (hi - lo); }
    public boolean contains(double x)     { return lo <= x && x <= hi; }
    public boolean encloses(Interval o)   { return lo <= o.lo && o.hi <= hi; }

    public Interval negate() { return new Interval(-hi, -lo); }

    public Interval add(Interval r) {
        return padded(lo + r.lo, hi + r.hi);
    }
    public Interval sub(Interval r) {
        return padded(lo - r.hi, hi - r.lo);
    }
    public Interval mul(Interval r) {
        final double a = lo * r.lo, b = lo * r.hi;
        final double c = hi * r.lo, d = hi * r.hi;
        return padded(Math.min(Math.min(a, b), Math.min(c, d)),
                      Math.max(Math.max(a, b), Math.max(c, d)));
    }
    public Interval mul(double s) {
        final double a = lo * s, b = hi * s;
        return padded(Math.min(a, b), Math.max(a, b));
    }

    private static double safeUlp(double x) {
        final double u = Math.ulp(x);
        return (u > 0.0) ? u : Double.MIN_NORMAL;
    }
    private static Interval padded(double lo, double hi) {
        return new Interval(lo - safeUlp(lo), hi + safeUlp(hi));
    }

    @Override
    public String toString() {
        return "[" + lo + ", " + hi + "]";
    }
}
