package Stats;

/**
 * Random noise generators used by tests and analysis: gaussian (Normal) and
 * white (Uniform-equivalent variance) noise scaled to a given standard
 * deviation.
 */
import java.util.Random;

public class Noise {
    Random rand = new Random();

    public double gaussian(final double dev) {
        return rand.nextGaussian() * Math.abs(dev);
    }

    public double white(final double dev) {
        return (rand.nextDouble() - 0.5) * Math.sqrt(12) * Math.abs(dev);
    }

}
