/**
 * This class computes the Earth Mover's Distance, using the EMD-HAT algorithm
 * created by Ofir Pele and Michael Werman.
 * 
 * This implementation is strongly based on the C++ code by the same authors,
 * that can be found here:
 * http://www.cs.huji.ac.il/~ofirpele/FastEMD/code/
 * 
 * Some of the author's comments on the original were kept or edited for 
 * this context.
 */

package imagingbook.jfastemd;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import imagingbook.jfastemd.Edges.Edge;
import imagingbook.jfastemd.Edges.Edge0;

/**
 * @author Telmo Menezes (telmo@telmomenezes.com)
 * @author Ofir Pele
 *
 */
public class JFastEMD {
	/**
	 * This interface is similar to Rubner's interface. See:
	 * http://www.cs.duke.edu/~tomasi/software/emd.htm
	 *
	 * To get the same results as Rubner's code you should set extra_mass_penalty to 0,
	 * and divide by the minimum of the sum of the two signature's weights. However, I
	 * suggest not to do this as you lose the metric property and more importantly, in my
	 * experience the performance is better with emd_hat. for more on the difference
	 * between emd and emd_hat, see the paper:
	 * A Linear Time Histogram Metric for Improved SIFT Matching
	 * Ofir Pele, Michael Werman
	 * ECCV 2008
	 *
	 * To get shorter running time, set the ground distance function to
	 * be a thresholded distance. For example: min(L2, T). Where T is some threshold.
	 * Note that the running time is shorter with smaller T values. Note also that
	 * thresholding the distance will probably increase accuracy. Finally, a thresholded
	 * metric is also a metric. See paper:
	 * Fast and Robust Earth Mover's Distances
	 * Ofir Pele, Michael Werman
	 * ICCV 2009
	 *
	 * If you use this code, please cite the papers.
	 */

	
	public double distance(Signature signature1, Signature signature2, double extraMassPenalty) {

		final int n1 = signature1.getNumberOfFeatures();
		final int n2 = signature2.getNumberOfFeatures();
		
		double[] P = new double[n1 + n2];
		double[] Q = new double[n1 + n2];

		for (int i = 0; i < n1; i++) {
			P[i] = signature1.getWeights()[i];
		}

		for (int j = 0; j < n2; j++) {
			Q[j + n1] = signature2.getWeights()[j];
		}

		double[][] C = new double[P.length][P.length];
		
		Feature[] features1 = signature1.getFeatures();
		Feature[] features2 = signature2.getFeatures();
		
		for (int i = 0; i < n1; i++) {
			for (int j = 0; j < n2; j++) {
				double dist = features1[i].groundDist(features2[j]);
				assert (dist >= 0);
				C[i][j + n1] = dist;
				C[j + n1][i] = dist;
			}
		}

		return emdHat(P, Q, C, extraMassPenalty);
	}

 // ------------------------------------------------------------------
    	
    private long emdHatImpl(long[] Pc, long[] Qc, long[][] C, long extraMassPenalty) {
		final int N = Pc.length;
		assert (Qc.length == N);
		
		System.out.println("emdHatImpl: N = " + N);

		// Ensuring that the supplier - P, have more mass.
		// Note that we assume here that C is symmetric
		long sumP = 0;
		for (int i = 0; i < N; i++) {
			sumP += Pc[i];
		}
		
		long sumQ = 0;
		for (int i = 0; i < N; i++) {
			sumQ += Qc[i];
		}
		
		long[] P, Q;
		long absDiffSumPSumQ;
		if (sumQ > sumP) {
			P = Qc;
			Q = Pc;
			absDiffSumPSumQ = sumQ - sumP;
		} else {
			P = Pc;
			Q = Qc;
			absDiffSumPSumQ = sumP - sumQ;
		}

		// creating the b vector that contains all vertexes
//		Vector<Long> b = new Vector<Long>();
		long[] b = new long[2 * N + 2];
//		for (int i = 0; i < 2 * N + 2; i++) {
//			b.add(0l);
//		}
		
		int THRESHOLD_NODE = 2 * N;
		int ARTIFICIAL_NODE = 2 * N + 1; // need to be last !
		
		for (int i = 0; i < N; i++) {
//			b.set(i, P.get(i));
//			b.set(i, P[i]);
			b[i] = P[i];
		}
		for (int i = N; i < 2 * N; i++) {
//			b.set(i, Q.get(i - N));
//			b.set(i, Q[i - N]);
			b[i] = Q[i - N];
		}

		// remark*) I put here a deficit of the extra mass, as mass that flows
		// to the threshold node
		// can be absorbed from all sources with cost zero (this is in reverse
		// order from the paper,
		// where incoming edges to the threshold node had the cost of the
		// threshold and outgoing
		// edges had the cost of zero)
		// This also makes sum of b zero.
//		b.set(THRESHOLD_NODE, -absDiffSumPSumQ);
//		b.set(ARTIFICIAL_NODE, 0l);
		b[THRESHOLD_NODE] = -absDiffSumPSumQ;
		b[ARTIFICIAL_NODE] = 0;

		long maxC = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				assert (C[i][j] >= 0);
				if (C[i][j] > maxC) {
					maxC = C[i][j];
				}
			}
		}
		if (extraMassPenalty == -1)
			extraMassPenalty = maxC;

		Set<Integer> sourcesThatFlowNotOnlyToThresh = new HashSet<Integer>();
		Set<Integer> sinksThatGetFlowNotOnlyFromThresh = new HashSet<Integer>();
		long preFlowCost = 0;
		
//		System.out.println("emdHatImpl:  b.size() = " +  b.size());
		System.out.println("emdHatImpl:  b.length = " +  b.length);

		// regular edges between sinks and sources without threshold edges
//		Vector<List<Edge>> c = new Vector<List<Edge>>();
		@SuppressWarnings("unchecked")
		List<Edge>[] c = new LinkedList[b.length];
		for (int i = 0; i < b.length; i++) {	// for (int i = 0; i < b.size(); i++) {
			//c.add(new LinkedList<Edge>());
			c[i] = new LinkedList<Edge>();
		}
		for (int i = 0; i < N; i++) {
			if (b[i] == 0)	// if (b.get(i) == 0)
				continue;
			for (int j = 0; j < N; j++) {
				if (b[j + N] == 0)	// if (b.get(j + N) == 0)
					continue;
				if (C[i][j] == maxC)
					continue;
//				c.get(i).add(new Edge(j + N, C[i][j]));
				c[i].add(new Edge(j + N, C[i][j]));
			}
		}

		// checking which are not isolated
		for (int i = 0; i < N; i++) {
			if (b[i] == 0)	// if (b.get(i) == 0)
				continue;
			for (int j = 0; j < N; j++) {
				if (b[j + N] == 0)	// if (b.get(j + N) == 0)
					continue;
				if (C[i][j] == maxC)
					continue;
				sourcesThatFlowNotOnlyToThresh.add(i);
				sinksThatGetFlowNotOnlyFromThresh.add(j + N);
			}
		}

		// converting all sinks to negative
		for (int i = N; i < 2 * N; i++) {
//			b.set(i, -b.get(i));
			b[i] = -b[i];
		}

		// add edges from/to threshold node,
		// note that costs are reversed to the paper (see also remark* above)
		// It is important that it will be this way because of remark* above.
		for (int i = 0; i < N; ++i) {
//			c.get(i).add(new Edge(THRESHOLD_NODE, 0));
			c[i].add(new Edge(THRESHOLD_NODE, 0));
		}
		for (int j = 0; j < N; ++j) {
//			c.get(THRESHOLD_NODE).add(new Edge(j + N, maxC));
			c[THRESHOLD_NODE].add(new Edge(j + N, maxC));
		}

		// artificial arcs - Note the restriction that only one edge i,j is
		// artificial so I ignore it...
		for (int i = 0; i < ARTIFICIAL_NODE; i++) {
//			c.get(i).add(new Edge(ARTIFICIAL_NODE, maxC + 1));
			c[i].add(new Edge(ARTIFICIAL_NODE, maxC + 1));
//			c.get(ARTIFICIAL_NODE).add(new Edge(i, maxC + 1));
			c[ARTIFICIAL_NODE].add(new Edge(i, maxC + 1));
		}


		// Note here it should be vector<int> and not vector<int>
		// as I'm using -1 as a special flag !!!
		final int REMOVE_NODE_FLAG = -1;
		
//		Vector<Integer> nodesNewNames = new Vector<Integer>();
		int[] nodesNewNames = new int[b.length];
//		Vector<Integer> nodesOldNames = new Vector<Integer>();
		
		for (int i = 0; i < b.length; i++) {	// for (int i = 0; i < b.size(); i++) {
//			nodesNewNames.add(REMOVE_NODE_FLAG);
			nodesNewNames[i] = REMOVE_NODE_FLAG;
//			nodesOldNames.add(0);
		}
		
		// remove nodes with supply demand of 0
		// and vertexes that are connected only to the
		// threshold vertex
		int currentNodeName = 0;
		
		for (int i = 0; i < N * 2; i++) {
			if (b[i] != 0) {	// if (b.get(i) != 0) {
				if (sourcesThatFlowNotOnlyToThresh.contains(i)
						|| sinksThatGetFlowNotOnlyFromThresh.contains(i)) {
//					nodesNewNames.set(i, currentNodeName);
					nodesNewNames[i] = currentNodeName;
//					nodesOldNames.add(i);
					currentNodeName++;
				} else {
					if (i >= N) {
//						preFlowCost -= (b.get(i) * maxC);
						preFlowCost -= (b[i] * maxC);
					}
//					b.set(THRESHOLD_NODE, b.get(THRESHOLD_NODE) + b.get(i)); // add mass(i<N) or deficit (i>=N)
					b[THRESHOLD_NODE] = b[THRESHOLD_NODE] + b[i]; // add mass(i<N) or deficit (i>=N)
				}
			}
		}
//		nodesNewNames.set(THRESHOLD_NODE, currentNodeName);
		nodesNewNames[THRESHOLD_NODE] = currentNodeName;
//		nodesOldNames.add(THRESHOLD_NODE);
		currentNodeName++;
//		nodesNewNames.set(ARTIFICIAL_NODE, currentNodeName);
		nodesNewNames[ARTIFICIAL_NODE] = currentNodeName;
//		nodesOldNames.add(ARTIFICIAL_NODE);
		currentNodeName++;
		
		
//		System.out.println("emdHatImpl: nodesNewNames.size() = " +  nodesNewNames.size());
		System.out.println("emdHatImpl: nodesNewNames.length = " +  nodesNewNames.length);
//		System.out.println("emdHatImpl: nodesOldNames.size() = " +  nodesOldNames.size());
		System.out.println("emdHatImpl: currentNodeName = " +  currentNodeName);

//		Vector<Long> bb = new Vector<Long>();
		long[] bb = new long[currentNodeName];
//		for (int i = 0; i < currentNodeName; i++) {
//			bb.add(0l);
//		}
		int j = 0;
		for (int i = 0; i < b.length; i++) {	// for (int i = 0; i < b.size(); i++) {
			if (nodesNewNames[i] != REMOVE_NODE_FLAG) {	// if (nodesNewNames.get(i) != REMOVE_NODE_FLAG) {
//				bb.set(j, b.get(i));
//				bb.set(j, b[i]);
				bb[j] = b[i];
				j++;
			}
		}

		Vector<List<Edge>> cc = new Vector<List<Edge>>();
		for (int i = 0; i < bb.length; i++) {	// for (int i = 0; i < bb.size(); i++) {
			cc.add(new LinkedList<Edge>());
		}
		for (int i = 0; i < c.length; i++) {	// for (int i = 0; i < c.size(); i++) {
			if (nodesNewNames[i] == REMOVE_NODE_FLAG)	// if (nodesNewNames.get(i) == REMOVE_NODE_FLAG)
				continue;
			for (Edge it : c[i]) {	// for (Edge it : c.get(i)) {
//				if (nodesNewNames.get(it.to) != REMOVE_NODE_FLAG) {
//					cc.get(nodesNewNames.get(i)).add(
//							new Edge(nodesNewNames.get(it.to), it.cost));
//				}
				if (nodesNewNames[it.to] != REMOVE_NODE_FLAG) {
					cc.get(nodesNewNames[i]).add(
							new Edge(nodesNewNames[it.to], it.cost));
				}
			}
		}

//		long myDist;

//		Vector<List<Edge0>> flows = new Vector<List<Edge0>>(bb.size());
		Vector<List<Edge0>> flows = new Vector<List<Edge0>>(bb.length);
		for (int i = 0; i < bb.length; i++) {	// for (int i = 0; i < bb.size(); i++) {
			flows.add(new LinkedList<Edge0>());
		}
		
		MinCostFlow mcf = new MinCostFlow();

		long mcfDist = mcf.compute(bb, cc, flows);

		long myDist = preFlowCost + // pre-flowing on cases where it was possible
				mcfDist + // solution of the transportation problem
				(absDiffSumPSumQ * extraMassPenalty); // emd-hat extra mass penalty

		return myDist;
	}

//	private double emdHat(Vector<Double> P, Vector<Double> Q, Vector<Vector<Double>> C, double extraMassPenalty) {
//		// This condition should hold:
//		// ( 2^(sizeof(CONVERT_TO_T*8)) >= ( MULT_FACTOR^2 )
//		// Note that it can be problematic to check it because
//		// of overflow problems. I simply checked it with Linux calc
//		// which has arbitrary precision.
//		double MULT_FACTOR = 1000000;
//
//		// Constructing the input
//		int N = P.size();
//		Vector<Long> iP = new Vector<Long>();
//		Vector<Long> iQ = new Vector<Long>();
//		Vector<Vector<Long>> iC = new Vector<Vector<Long>>();
//		for (int i = 0; i < N; i++) {
//			iP.add(0l);
//			iQ.add(0l);
//			Vector<Long> vec = new Vector<Long>();
//			for (int j = 0; j < N; j++) {
//				vec.add(0l);
//			}
//			iC.add(vec);
//		}
//
//		// Converting to CONVERT_TO_T
//		double sumP = 0.0;
//		double sumQ = 0.0;
//		double maxC = C.get(0).get(0);
//		for (int i = 0; i < N; i++) {
//			sumP += P.get(i);
//			sumQ += Q.get(i);
//			for (int j = 0; j < N; j++) {
//				if (C.get(i).get(j) > maxC)
//					maxC = C.get(i).get(j);
//			}
//		}
//		double minSum = Math.min(sumP, sumQ);
//		double maxSum = Math.max(sumP, sumQ);
//		double PQnormFactor = MULT_FACTOR / maxSum;
//		double CnormFactor = MULT_FACTOR / maxC;
//		for (int i = 0; i < N; i++) {
//			iP.set(i, (long) (Math.floor(P.get(i) * PQnormFactor + 0.5)));
//			iQ.set(i, (long) (Math.floor(Q.get(i) * PQnormFactor + 0.5)));
//			for (int j = 0; j < N; j++) {
//				iC.get(i)
//				.set(j,
//						(long) (Math.floor(C.get(i).get(j)
//								* CnormFactor + 0.5)));
//			}
//		}
//
//		// computing distance without extra mass penalty
//		double dist = emdHatImpl(iP, iQ, iC, 0);
//		// unnormalize
//		dist = dist / PQnormFactor;
//		dist = dist / CnormFactor;
//
//		// adding extra mass penalty
//		if (extraMassPenalty == -1) {
//			extraMassPenalty = maxC;
//		}
//		dist += (maxSum - minSum) * extraMassPenalty;
//
//		return dist;
//	}
    
//    @Deprecated
//    private double emdHat(Vector<Double> P, Vector<Double> Q, Vector<Vector<Double>> C, double extraMassPenalty) {
//    	return emdHat(toArrayD(P), toArrayD(Q), toMatrixD(C), extraMassPenalty);
//    }
    
	private double emdHat(double[] P, double[] Q, double[][] C, double extraMassPenalty) {
		// This condition should hold:
		// ( 2^(sizeof(CONVERT_TO_T*8)) >= ( MULT_FACTOR^2 )
		// Note that it can be problematic to check it because
		// of overflow problems. I simply checked it with Linux calc
		// which has arbitrary precision.
		final double MULT_FACTOR = 1000000;

		// Constructing the input
//		final int N = P.size();
		final int N = P.length;
		
//		Vector<Long> iP = new Vector<Long>();
//		Vector<Long> iQ = new Vector<Long>();
		long[] iP = new long[N];
		long[] iQ = new long[N];
		
//		Vector<Vector<Long>> iC = new Vector<Vector<Long>>();
		long[][] iC = new long[N][N];
//		for (int i = 0; i < N; i++) {
//			iP[i] = 0l; // iP.add(0l);	// not needed!
//			iP[i] = 0l; // iQ.add(0l);
//			
//			Vector<Long> vec = new Vector<Long>();
//			for (int j = 0; j < N; j++) {
//				vec.add(0l);
//			}
//			iC.add(vec);
//		}

		// Converting to CONVERT_TO_T
		double sumP = 0.0;
		double sumQ = 0.0;
//		double maxC = C.get(0).get(0);
		double maxC = C[0][0];
		for (int i = 0; i < N; i++) {
//			sumP += P.get(i);
//			sumQ += Q.get(i);
			sumP += P[i];
			sumQ += Q[i];
			for (int j = 0; j < N; j++) {
//				if (C.get(i).get(j) > maxC) {
//					maxC = C.get(i).get(j);
//				}
				if (C[i][j] > maxC) {
					maxC = C[i][j];
				}
			}
		}
		double minSum = Math.min(sumP, sumQ);
		double maxSum = Math.max(sumP, sumQ);
		double PQnormFactor = MULT_FACTOR / maxSum;
		double CnormFactor = MULT_FACTOR / maxC;
		for (int i = 0; i < N; i++) {
//			iP.set(i, (long) (Math.floor(P.get(i) * PQnormFactor + 0.5)));
//			iQ.set(i, (long) (Math.floor(Q.get(i) * PQnormFactor + 0.5)));
//			iP[i] = (long) (Math.floor(P.get(i) * PQnormFactor + 0.5));
//			iQ[i] = (long) (Math.floor(Q.get(i) * PQnormFactor + 0.5));
			iP[i] = (long) (Math.floor(P[i] * PQnormFactor + 0.5));
			iQ[i] = (long) (Math.floor(Q[i] * PQnormFactor + 0.5));
			for (int j = 0; j < N; j++) {
//				iC.get(i).set(j, (long) (Math.floor(C.get(i).get(j) * CnormFactor + 0.5)));
//				iC[i][j] = (long) (Math.floor(C.get(i).get(j) * CnormFactor + 0.5));
				iC[i][j] = (long) (Math.floor(C[i][j] * CnormFactor + 0.5));
			}
		}

		// computing distance without extra mass penalty
		double dist = emdHatImpl(iP, iQ, iC, 0);
		// unnormalize
		dist = dist / PQnormFactor;
		dist = dist / CnormFactor;

		// adding extra mass penalty
		if (extraMassPenalty == -1) {
			extraMassPenalty = maxC;
		}
		dist += (maxSum - minSum) * extraMassPenalty;

		return dist;
	}
}