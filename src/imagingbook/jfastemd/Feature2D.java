package imagingbook.jfastemd;

/**
 * @author Telmo Menezes (telmo@telmomenezes.com)
 *
 */
public class Feature2D implements Feature {
    private final double x;
    private final double y;

    public Feature2D(double x, double y) {
        this.x = x;
        this.y = y;
    }
    
    @Override
    public double groundDist(Feature f) {
    	if (f instanceof Feature2D) {
	        Feature2D f2d = (Feature2D)f;
	        double deltaX = x - f2d.x;
	        double deltaY = y - f2d.y;
	        return Math.sqrt((deltaX * deltaX) + (deltaY * deltaY));
    	}
    	else {
    		throw new IllegalArgumentException("groundDist(): non-matching feature types!");
    	}
    }
}