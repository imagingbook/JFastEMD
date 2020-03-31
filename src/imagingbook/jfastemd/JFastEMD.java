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
//import java.util.Vector;

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
		
//		System.out.println("emdHatImpl: N = " + N);

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
		long[] b = new long[2 * N + 2];
		int THRESHOLD_NODE = 2 * N;
		int ARTIFICIAL_NODE = 2 * N + 1; // need to be last !
		
		for (int i = 0; i < N; i++) {
			b[i] = P[i];
		}
		for (int i = N; i < 2 * N; i++) {
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

		// regular edges between sinks and sources without threshold edges
		@SuppressWarnings("unchecked")
		List<Edge>[] c = new LinkedList[b.length];
		for (int i = 0; i < b.length; i++) {
			c[i] = new LinkedList<Edge>();
		}
		for (int i = 0; i < N; i++) {
			if (b[i] == 0)
				continue;
			for (int j = 0; j < N; j++) {
				if (b[j + N] == 0)
					continue;
				if (C[i][j] == maxC)
					continue;
				c[i].add(new Edge(j + N, C[i][j]));
			}
		}

		// checking which are not isolated
		for (int i = 0; i < N; i++) {
			if (b[i] == 0)
				continue;
			for (int j = 0; j < N; j++) {
				if (b[j + N] == 0)
					continue;
				if (C[i][j] == maxC)
					continue;
				sourcesThatFlowNotOnlyToThresh.add(i);
				sinksThatGetFlowNotOnlyFromThresh.add(j + N);
			}
		}

		// converting all sinks to negative
		for (int i = N; i < 2 * N; i++) {
			b[i] = -b[i];
		}

		// add edges from/to threshold node,
		// note that costs are reversed to the paper (see also remark* above)
		// It is important that it will be this way because of remark* above.
		for (int i = 0; i < N; ++i) {
			c[i].add(new Edge(THRESHOLD_NODE, 0));
		}
		for (int j = 0; j < N; ++j) {
			c[THRESHOLD_NODE].add(new Edge(j + N, maxC));
		}

		// artificial arcs - Note the restriction that only one edge i,j is
		// artificial so I ignore it...
		for (int i = 0; i < ARTIFICIAL_NODE; i++) {
			c[i].add(new Edge(ARTIFICIAL_NODE, maxC + 1));
			c[ARTIFICIAL_NODE].add(new Edge(i, maxC + 1));
		}


		// Note here it should be vector<int> and not vector<int>
		// as I'm using -1 as a special flag !!!
		final int REMOVE_NODE_FLAG = -1;
		int[] nodesNewNames = new int[b.length];
		
		for (int i = 0; i < b.length; i++) {
			nodesNewNames[i] = REMOVE_NODE_FLAG;
		}
		
		// remove nodes with supply demand of 0
		// and vertexes that are connected only to the
		// threshold vertex
		int currentNodeName = 0;
		
		for (int i = 0; i < N * 2; i++) {
			if (b[i] != 0) {
				if (sourcesThatFlowNotOnlyToThresh.contains(i)
						|| sinksThatGetFlowNotOnlyFromThresh.contains(i)) {
					nodesNewNames[i] = currentNodeName;
					currentNodeName++;
				} else {
					if (i >= N) {
						preFlowCost -= (b[i] * maxC);
					}
					b[THRESHOLD_NODE] = b[THRESHOLD_NODE] + b[i]; // add mass(i<N) or deficit (i>=N)
				}
			}
		}

		nodesNewNames[THRESHOLD_NODE] = currentNodeName;
		currentNodeName++;
		nodesNewNames[ARTIFICIAL_NODE] = currentNodeName;
		currentNodeName++;

		long[] bb = new long[currentNodeName];
		
		int j = 0;
		for (int i = 0; i < b.length; i++) {
			if (nodesNewNames[i] != REMOVE_NODE_FLAG) {
				bb[j] = b[i];
				j++;
			}
		}

		@SuppressWarnings("unchecked")
		List<Edge>[] cc = new LinkedList[bb.length];
		for (int i = 0; i < bb.length; i++) {
			cc[i] = new LinkedList<Edge>();
		}
		for (int i = 0; i < c.length; i++) {
			if (nodesNewNames[i] == REMOVE_NODE_FLAG)
				continue;
			for (Edge it : c[i]) {
				if (nodesNewNames[it.to] != REMOVE_NODE_FLAG) {
					cc[nodesNewNames[i]].add(new Edge(nodesNewNames[it.to], it.cost));
				}
			}
		}

		@SuppressWarnings("unchecked")
		List<Edge0>[] flows = new LinkedList[bb.length];
		for (int i = 0; i < bb.length; i++) {
			flows[i] = new LinkedList<Edge0>();
		}
		
		MinCostFlow mcf = new MinCostFlow();
		long mcfDist = mcf.compute(bb, cc, flows);
		long myDist = preFlowCost + // pre-flowing on cases where it was possible
				mcfDist + // solution of the transportation problem
				(absDiffSumPSumQ * extraMassPenalty); // emd-hat extra mass penalty

		return myDist;
	}
    
	private double emdHat(double[] P, double[] Q, double[][] C, double extraMassPenalty) {
		// This condition should hold:
		// ( 2^(sizeof(CONVERT_TO_T*8)) >= ( MULT_FACTOR^2 )
		// Note that it can be problematic to check it because
		// of overflow problems. I simply checked it with Linux calc
		// which has arbitrary precision.
		final double MULT_FACTOR = 1000000;

		// Constructing the input
		final int N = P.length;
		long[] iP = new long[N];
		long[] iQ = new long[N];
		
		long[][] iC = new long[N][N];

		// Converting to CONVERT_TO_T
		double sumP = 0.0;
		double sumQ = 0.0;
		double maxC = C[0][0];
		for (int i = 0; i < N; i++) {
			sumP += P[i];
			sumQ += Q[i];
			for (int j = 0; j < N; j++) {
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
			iP[i] = Math.round(P[i] * PQnormFactor);
			iQ[i] = Math.round(Q[i] * PQnormFactor);
			for (int j = 0; j < N; j++) {
				iC[i][j] = Math.round(C[i][j] * CnormFactor);
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