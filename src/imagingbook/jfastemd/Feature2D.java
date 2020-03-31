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
    
    public double groundDist(Feature f) {
        Feature2D f2d = (Feature2D)f;	// TODO: this sucks
        double deltaX = x - f2d.x;
        double deltaY = y - f2d.y;
        return Math.sqrt((deltaX * deltaX) + (deltaY * deltaY));
    }
}