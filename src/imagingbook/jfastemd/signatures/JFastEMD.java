package imagingbook.jfastemd.signatures;

import imagingbook.jfastemd.JFastEMD2;

public abstract class JFastEMD {
	
	public static JFastEMD2 create(Signature signature1, Signature signature2, double extraMassPenalty) {
		int n1 = signature1.getNumberOfFeatures();
		int n2 = signature2.getNumberOfFeatures();
		double[] P = new double[n1 + n2];
		double[] Q = new double[n1 + n2];
		double[][] C = new double[P.length][P.length];
	
		for (int i = 0; i < n1; i++) {
			P[i] = signature1.getWeights()[i];
		}

		for (int j = 0; j < n2; j++) {
			Q[j + n1] = signature2.getWeights()[j];
		}

		Feature[] f1 = signature1.getFeatures();
		Feature[] f2 = signature2.getFeatures();

		for (int i = 0; i < n1; i++) {
			for (int j = 0; j < n2; j++) {
				double dist = f1[i].groundDist(f2[j]);
				C[i][j + n1] = dist;
				C[j + n1][i] = dist;
			}
		}
		
		return new JFastEMD2(P, Q, C, extraMassPenalty);
	}
	
	public static JFastEMD2 create(Signature signature1, Signature signature2) {
//		return create(signature1, signature2, -1);
		return create(signature1, signature2, 0);
	}

}
