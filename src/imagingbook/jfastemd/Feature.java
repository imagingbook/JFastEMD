package imagingbook.jfastemd;

import java.util.Arrays;

/**
 * This class represents a n-dimensional feature (position in feature space).
 * @author WB
 *
 */
public class Feature {
	
    private final double[] x;

    public Feature(double ... x) {
        this.x = x;
    }
    
    public double[] getX() {
    	return this.x;
    }
    
    public double getElement(int k) {
    	return this.x[k];
    }
    
    public double groundDist(Feature other) {
    	if (this.x.length != other.x.length) {
    		throw new IllegalArgumentException("groundDist(): non-matching feature lengths!");
    	}
    	double sum = 0;
    	for (int i = 0; i < x.length; i++) {
	        double d = this.x[i] - other.x[i];
	        sum = sum + d * d;
    	}
    	return Math.sqrt(sum);
    }
    
    public String toString() {
    	return Arrays.toString(x);
    }
    
}