package Stats;

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
