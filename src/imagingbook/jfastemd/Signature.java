package imagingbook.jfastemd;


/**
 * @author Telmo Menezes (telmo@telmomenezes.com)
 * @author WILBUR (modified)
 * 
 *
 */
public class Signature {
	
    private final Feature[] features;
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

}