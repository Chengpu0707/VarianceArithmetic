package Stats;

import java.util.Random;

public class Noise {
    Random rand = new Random();

    double gaussian(final double dev) {
        return rand.nextGaussian() * Math.abs(dev);
    }
}
