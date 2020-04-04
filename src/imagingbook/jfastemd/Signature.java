package imagingbook.jfastemd;

/**
 * @author Telmo Menezes (telmo@telmomenezes.com)
 * @author WILBUR (modified)
 * 
 *
 */
public class Signature {
	
    private final Feature[] features;	//
    private final double[] weights;
    
    public Signature(Feature[] features, double[] weights) {
    	if (features.length != weights.length)
    		throw new IllegalArgumentException("Must be of same length: features, weights");
    	this.features = features;
    	this.weights = weights;
    }
    
    public int getNumberOfFeatures() {
    	return features.length;
    }

    public Feature[] getFeatures() {
        return features;
    }

    public double[] getWeights() {
        return weights;
    }
    
    public void print(String name) {
    	System.out.println("Signature: " + name);
    	for (int i = 0; i < features.length; i++) {
    		System.out.format("   %3d: %s w = %.3f\n", i, features[i].toString(), weights[i]);
    	}
    }
    
    // ----------------------------------------------------------------------
    
    /**
     * Creates a {@link Signature} from a 1-dimensional map.
     * @param map a 1D sequence of weight values
     * @return a new {@link Signature} instance
     */
    public static Signature create(double[] map) {
    	int n = 0;
    	for (int x = 0; x < map.length; x++) {
    		if (map[x] > 0) n++;
    	}

    	Feature[] features = new Feature[n];
    	double[] weights = new double[n];

    	int i = 0;
    	for (int x = 0; x < map.length; x++) {
    		double val = map[x];
    		if (val > 0) {
    			features[i] = new Feature(x);
    			weights[i] = val;
    			i++;
    		}
    	}
    	return new Signature(features, weights);
    }
    
    /**
     * Creates a {@link Signature} from a 2-dimensional map.
     * @param map a 2D array of weight values (assumed to be rectangular)
     * @return a new {@link Signature} instance
     */
    public static Signature create(double[][] map) {
    	final int rows = map.length;
    	final int cols = map[0].length;
    	int n = 0;
    	for (int y = 0; y < rows; y++) {
	    	for (int x = 0; x < cols; x++) {
	    		if (map[x][y] > 0) n++;
	    	}
    	}

    	Feature[] features = new Feature[n];
    	double[] weights = new double[n];

    	int i = 0;
    	for (int y = 0; y < rows; y++) {
	    	for (int x = 0; x < cols; x++) {
	    		double val = map[x][y];
	    		if (val > 0) {
	    			features[i] = new Feature(x, y);
	    			weights[i] = val;
	    			i++;
	    		}
	    	}
    	}
    	return new Signature(features, weights);
    }

}