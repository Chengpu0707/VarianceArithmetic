package Func;

import java.util.HashMap;
import java.util.Map;

import Stats.Histogram;
import Stats.Stat;

enum RealType {
    Var,
    Intv
}

enum NoiseType {
    Gaussian,
    White
}

enum TestType {
    Forward,
    Reverse,
    Roundtrip,
}

class Measure {
    static final int BINDING = 3, DIVIDS = 5;

    final Map<TestType, Stat> sStat = new HashMap<>();
    final Map<TestType, Histogram> sHisto = new HashMap<>();

    Measure() {
        for (TestType test: TestType.values()) {
            sStat.put(test, new Stat());
            sHisto.put(test, new Histogram(BINDING, DIVIDS));
        }
    }
}


