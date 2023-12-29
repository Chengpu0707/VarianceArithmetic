package fitnesse;

import Type.VarDbl;

public class SlimVarDbl extends VarDbl {
    private double initValue;
    private double initDev;

    private double value;
    private double dev;
    static private double valueTolerance;
    static private double devTolerance;

    public void setValueTolerance(double tolerance) {
        valueTolerance = tolerance;
    }
    
    public void setDevTolerance(double tolerance) {
        devTolerance = tolerance;
    }
    
    public void setComment(String comment) {
    }

    public void setInitValue(double value) {
        initValue = value;
    }

    public void setInitDev(double dev) {
        initDev = dev;
    }
    
    private String compare() {
        try {
            if ((valueTolerance <= 0) ||  (value == 0)) {
                if (0 != value())
                   return String.format("The value {0} differs from 0", value());
            } else {
                if (Math.abs(value() / value - 1) > valueTolerance)
                   return String.format("The value {0} differs from the init {1} by more than {2}", 
                        value(), value, valueTolerance);
            }

            if ((devTolerance <= 0) || (dev == 0)) {
                if (0 != uncertainty())
                   return String.format("The dev {0} differs from 0", uncertainty());
            } else {
                if (Math.abs(uncertainty() / dev - 1) > devTolerance)
                   return String.format("The dev {0} differs from the init {1} by more than {2}", 
                        uncertainty(), dev, devTolerance);
            }

             if ((0 < valueTolerance) && (0 < devTolerance) && (0 != value) && (0 < dev)) {
                if(Math.abs(precSq() / (dev / value) - 1) > Math.max(valueTolerance, devTolerance))
                    return null;
            }
        } catch (ValueException e) {
            return e.toString();
        } catch (UncertaintyException e) {
            return e.toString();
        }
        return null;
    }

    public String init() {
        try {
            pack(initValue, initDev * initDev, false, false, 255);
        } catch (ValueException e) {
            return e.toString();
        } catch (UncertaintyException e) {
            return e.toString();
        }
        value = initValue;
        dev = initDev;
        return compare();
    }


    public static void main(String[] args) {
        final SlimVarDbl fix = new SlimVarDbl();
        fix.setInitValue(0);
    }

}
