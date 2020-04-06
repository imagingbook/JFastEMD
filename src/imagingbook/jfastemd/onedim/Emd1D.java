package imagingbook.jfastemd.onedim;

/**
 * https://commons.apache.org/proper/commons-math/javadocs/api-3.6/src-html/org/apache/commons/math3/ml/distance/EarthMoversDistance.html
 * 
 * @author WB
 *
 */
public class Emd1D {

	// Apache version
	static double compute0(double[] P, double[] Q) {
		double EMDi = 0;
		double totalDistance = 0;
		for (int i = 0; i < P.length; i++) {
			double di = (P[i] + EMDi) - Q[i];
			totalDistance += Math.abs(di);
			EMDi = di;
		}
		return totalDistance;
	}
	
	// Apache version renamed variables
	static double compute1(double[] P, double[] Q) {
		double d = 0;
		double emd = 0;
		for (int i = 0; i < P.length; i++) {
			d = (P[i] + d) - Q[i];
			emd += Math.abs(d);
		}
		return emd;
	}

	// Wilbur's version A
	static double compute2(double[] P, double[] Q) {
		double dk = 0;
		double emd = 0;
		for (int k = 1; k <= P.length; k++) {
			dk = dk + P[k-1] - Q[k-1];
			emd += Math.abs(dk);
		}
		return emd;
	}
	
	// Wilbur's version B
	static double compute3(double[] P, double[] Q) {
		double dk = P[0] - Q[0];
		double emd = Math.abs(dk);
		for (int k = 1; k < P.length; k++) {
			dk = dk + P[k] - Q[k];
			emd += Math.abs(dk);
		}
		return emd;
	}
	
	// Wilbur's version C
	static double compute4(double[] P, double[] Q) {
		double PP = 0;	// cumulative P
		double QQ = 0;	// cumulative Q
		double emd = 0;
		for (int k = 0; k < P.length; k++) {
			PP = PP + P[k];
			QQ = QQ + Q[k];
			emd += Math.abs(PP - QQ);
		}
		return emd;
	}
	
	
	static double[] P = {0, 1, 3, 4, 0, 9};
	static double[] Q = {3, 1, 4, 0, 4, 5};
	
	public static void main(String[] args) {
		
		System.out.println("compute1 = " + compute1(P, P));
		System.out.println("compute1 = " + compute1(Q, Q));
		System.out.println("compute1 = " + compute1(P, Q));
		System.out.println("compute1 = " + compute1(Q, P));
		
		System.out.println();
		System.out.println("compute1 = " + compute2(P, P));
		System.out.println("compute1 = " + compute2(Q, Q));
		System.out.println("compute2 = " + compute2(P, Q));
		System.out.println("compute2 = " + compute2(Q, P));
		
		System.out.println();
		System.out.println("compute1 = " + compute3(P, P));
		System.out.println("compute1 = " + compute3(Q, Q));
		System.out.println("compute2 = " + compute3(P, Q));
		System.out.println("compute2 = " + compute3(Q, P));
		
		System.out.println();
		System.out.println("compute1 = " + compute4(P, P));
		System.out.println("compute1 = " + compute4(Q, Q));
		System.out.println("compute2 = " + compute4(P, Q));
		System.out.println("compute2 = " + compute4(Q, P));
		
	}

}
