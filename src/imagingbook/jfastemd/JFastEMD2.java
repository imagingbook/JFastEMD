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

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import imagingbook.jfastemd.Edges.Edge;
import imagingbook.jfastemd.Edges.EdgeWithFlow;

/**
 * This is a refactored version of the <strong>JFastEMD</strong> Java
 * implementation by Telmo Menezes, hosted at <a href=
 * "https://github.com/telmomenezes/JFastEMD">https://github.com/telmomenezes/JFastEMD</a>.
 * Most importantly, {@code Vector} data structures have been reduced to a
 * minimum and replaced by simple arrays wherever possible, which significantly
 * improves execution time (and memory requirements too). The original Java port
 * by Menezes was based on earlier implementations by Rubner (<a href=
 * "http://www.cs.duke.edu/~tomasi/software/emd.htm">http://www.cs.duke.edu/~tomasi/software/emd.htm</a>)
 * and Pele
 * (<a href="http://ofirpele.droppages.com/">http://ofirpele.droppages.com/</a>,
 * as described in the original note below.
 * 
 * <p>
 * <strong>Original Note:</strong>
 * </p>
 * <p>
 * This interface is similar to Rubner's interface. See:
 * http://www.cs.duke.edu/~tomasi/software/emd.htm
 * </p>
 * <p>
 * To get the same results as Rubner's code you should set extra_mass_penalty to
 * 0, and divide by the minimum of the sum of the two signature's weights.
 * However, I suggest not to do this as you lose the metric property and more
 * importantly, in my experience the performance is better with emd_hat. for
 * more on the difference between emd and emd_hat, see the paper: A Linear Time
 * Histogram Metric for Improved SIFT Matching by Ofir Pele, Michael Werman
 * (ECCV 2008).
 * </p>
 * <p>
 * To get shorter running time, set the ground distance function to be a
 * thresholded distance. For example: min(L2, T). Where T is some threshold.
 * Note that the running time is shorter with smaller T values. Note also that
 * thresholding the distance will probably increase accuracy. Finally, a
 * thresholded metric is also a metric. See paper: Fast and Robust Earth Mover's
 * Distances by Ofir Pele, Michael Werman (ICCV 2009).
 * </p>
 * <p>
 * If you use this code, please cite the papers.
 * </p>
 * <p>
 * <strong>Original License:</strong>
 * </p><hr>
 * <pre>
IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
By downloading, copying, installing or using the software you agree to this license.  If you do not
agree to this license, do not download, install, copy or use the software.

LICENSE CONDITIONS

Copyright (c) 2009-2012, Telmo Menezes and Ofir Pele
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
 * Neither the name of the The Hebrew University of Jerusalem nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * </pre>
 * <hr>
 * @author Telmo Menezes (telmo@telmomenezes.com)
 * @author Ofir Pele
 * @author W. Burger (refactored 2020)
 * 
 * @version 2020/03/31
 *
 */
public class JFastEMD2 {
	
	private static final double MULT_FACTOR = 1000000;
	private static final int REMOVE_NODE_FLAG = -1;	// as I'm using -1 as a special flag !!!

	private final double[] P, Q;
	private final double[][] C;
	private final double extraMassPenalty;

	public JFastEMD2(Signature signature1, Signature signature2) {
		this(signature1, signature2, -1);
	}

	// TODO: remove Signature, use original map (array) as input.
	@Deprecated
	public JFastEMD2(Signature signature1, Signature signature2, double extraMassPenalty) {
		int n1 = signature1.getNumberOfFeatures();
		int n2 = signature2.getNumberOfFeatures();
		this.P = new double[n1 + n2];
		this.Q = new double[n1 + n2];
		this.C = new double[P.length][P.length];
		this.extraMassPenalty = extraMassPenalty;
		init(signature1, signature2);
	}

	private void init(Signature signature1, Signature signature2) {
		int n1 = signature1.getNumberOfFeatures();
		int n2 = signature2.getNumberOfFeatures();
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
	}
	
	public JFastEMD2(double[] P, double[] Q, double[][] C, double extraMassPenalty) {
		this.P = P;
		this.Q = Q;
		this.C = C;
		this.extraMassPenalty = extraMassPenalty;
	}

	// ------------------------------------------------------------------

	public double getDistance() {
		return emdHat(P, Q, C);
	}

	// ------------------------------------------------------------------

	private double emdHat(double[] P, double[] Q, double[][] C) {
		// This condition should hold:
		// ( 2^(sizeof(CONVERT_TO_T*8)) >= ( MULT_FACTOR^2 )
		// Note that it can be problematic to check it because
		// of overflow problems. I simply checked it with Linux calc
		// which has arbitrary precision.
		
		// Constructing the input
		final int N = P.length;
		
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

		long[] iP = new long[N];
		long[] iQ = new long[N];
		long[][] iC = new long[N][N];
		
		double minSum = Math.min(sumP, sumQ);
		double maxSum = Math.max(sumP, sumQ);
		double PQnormFactor = MULT_FACTOR / maxSum;
		double CnormFactor =  MULT_FACTOR / maxC;
		for (int i = 0; i < N; i++) {
			iP[i] = Math.round(P[i] * PQnormFactor);
			iQ[i] = Math.round(Q[i] * PQnormFactor);
			for (int j = 0; j < N; j++) {
				iC[i][j] = Math.round(C[i][j] * CnormFactor);
			}
		}

		// computing distance without extra mass penalty
		double dist = emdHat(iP, iQ, iC, 0);
		// unnormalize
		dist = dist / PQnormFactor / CnormFactor;

		// adding extra mass penalty
		double emp = (extraMassPenalty < 0) ? maxC : extraMassPenalty;
		dist += (maxSum - minSum) * emp;

		return dist;
	}

	// ------------------------------------------------------------------
	/**
	 * EMD implementation using integer (long) quantities.
	 */
	private long emdHat(long[] Pc, long[] Qc, long[][] C, long extraMassPenalty) {
		if (Pc.length != Qc.length) {
			throw new IllegalArgumentException("Pc, Qc must be of same length!");
		}
		
		final int N = Pc.length;

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
		final long[] b = new long[2 * N + 2];
		final int THRESHOLD_NODE = 2 * N;
		final int ARTIFICIAL_NODE = 2 * N + 1; // need to be last !

//		System.arraycopy(P, 0, b, 0, N);	
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

		Set<Integer> sourcesThatFlowNotOnlyToThresh = new HashSet<>();
		Set<Integer> sinksThatGetFlowNotOnlyFromThresh = new HashSet<>();

		// regular edges between sinks and sources without threshold edges
		@SuppressWarnings("unchecked")
		List<Edge>[] c = new LinkedList[b.length];
		
		for (int i = 0; i < b.length; i++) {
			c[i] = new LinkedList<Edge>();
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
				c[i].add(new Edge(j + N, C[i][j]));
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
		for (int i = 0; i < N; i++) {
			c[i].add(new Edge(THRESHOLD_NODE, 0));
			c[THRESHOLD_NODE].add(new Edge(i + N, maxC));
		}
		
//		for (int j = 0; j < N; ++j) {
//			c[THRESHOLD_NODE].add(new Edge(j + N, maxC));
//		}

		// artificial arcs - Note the restriction that only one edge i,j is
		// artificial so I ignore it...
		for (int i = 0; i < ARTIFICIAL_NODE; i++) {
			c[i].add(new Edge(ARTIFICIAL_NODE, maxC + 1));
			c[ARTIFICIAL_NODE].add(new Edge(i, maxC + 1));
		}
		
		int[] nodesNewNames = new int[b.length];
		Arrays.fill(nodesNewNames, REMOVE_NODE_FLAG);
//		for (int i = 0; i < b.length; i++) {
//			nodesNewNames[i] = REMOVE_NODE_FLAG;
//		}

		// remove nodes with supply demand of 0
		// and vertexes that are connected only to the
		// threshold vertex
		int currentNodeName = 0;
		long preFlowCost = 0;

		for (int i = 0; i < N * 2; i++) {
			if (b[i] != 0) {
				if (sourcesThatFlowNotOnlyToThresh.contains(i) || sinksThatGetFlowNotOnlyFromThresh.contains(i)) {
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

		final long[] bb = new long[currentNodeName];

		int j = 0;
		for (int i = 0; i < b.length; i++) {
			if (nodesNewNames[i] != REMOVE_NODE_FLAG) {
				bb[j] = b[i];
				j++;
			}
		}

		@SuppressWarnings("unchecked")
		List<Edge>[] cc = new List[bb.length];
		
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
	
		long mcfDist = new MinCostFlow(bb, cc).getDistance();
		
		long emp = (extraMassPenalty == -1) ? maxC : extraMassPenalty;
		
		return preFlowCost + 			// pre-flowing on cases where it was possible
			mcfDist + 					// solution of the transportation problem
			(absDiffSumPSumQ * emp); 	// emd-hat extra mass penalty
	}

}