package imagingbook.jfastemd;

import java.util.ArrayList;
import java.util.List;

import imagingbook.jfastemd.Edges.Edge;
import imagingbook.jfastemd.Edges.EdgeReducedBackward;
import imagingbook.jfastemd.Edges.EdgeReducedForward;
import imagingbook.jfastemd.Edges.EdgeWithFlow;


class MinCostFlow {

	private final long[] e;					// e - supply(positive) and demand(negative).
	private final List<Edge>[] c;			// c[i] - edges that goes from node i. first is the second node
	private final List<EdgeWithFlow>[] x;	// x - the flow is returned in it
	
	private final int numNodes;
	private final int[] nodesToQ;
	private final long distance;

	@SuppressWarnings("unchecked")
	MinCostFlow(long[] e, List<Edge>[] c) {
		if (e.length != c.length)
			throw new IllegalArgumentException("Lengths of e, c do not match!");
		this.numNodes = e.length;
		this.e = e;
		this.c = c;
		this.x = new List[numNodes];
		this.nodesToQ = new int[numNodes];
		this.distance = compute();
	}

	// ----------------------------------------------------------------------------

	public long getDistance() {
		return this.distance;
	}

	public List<EdgeWithFlow>[] getFlow() {	// TODO: return flow as a matrix?
		return x;
	}

	// ----------------------------------------------------------------------------

	private long compute() {
		for (int i = 0; i < e.length; i++) {
			x[i] = new ArrayList<EdgeWithFlow>(e.length);
		}

		// init flow
		for (int from = 0; from < numNodes; from++) {
			for (Edge it : c[from]) {
				x[from].add(new EdgeWithFlow(it.to, it.cost, 0));
				x[it.to].add(new EdgeWithFlow(from, -it.cost, 0));
			}
		}

		// reduced costs for forward edges (c[i,j]-pi[i]+pi[j])
		// Note that for forward edges the residual capacity is infinity
		@SuppressWarnings("unchecked")
		List<EdgeReducedForward>[] rCostForward = new List[numNodes];
		@SuppressWarnings("unchecked")
		List<EdgeReducedBackward>[] rCostCapBackward = new List[numNodes];
		
		for (int i = 0; i < numNodes; i++) {
			int n = c[i].size();
			rCostForward[i]     = new ArrayList<EdgeReducedForward>(n);
			rCostCapBackward[i] = new ArrayList<EdgeReducedBackward>(n);
		}
		
		// reduced costs and capacity for backward edges
		// (c[j,i]-pi[j]+pi[i])
		// Since the flow at the beginning is 0, the residual capacity is
		// also zero
		for (int from = 0; from < numNodes; from++) {
			for (Edge it : c[from]) {
				rCostForward[from].add(new EdgeReducedForward(it.to, it.cost));
				rCostCapBackward[it.to].add(new EdgeReducedBackward(from, -it.cost, 0));
			}
		}

		// Max supply TODO:demand?, given U?, optimization-> min out of
		// demand,supply
		long U = 0;
		for (int i = 0; i < numNodes; i++) {
			if (e[i] > U) {
				U = e[i];
			}
		}

		final long[] d   = new long[numNodes];	// = 0L
		final int[] prev = new int[numNodes];	

		long delta = 1;
		while (true) { // until we break when S or T is empty
			long maxSupply = 0;
			int k = 0;
			for (int i = 0; i < numNodes; i++) {
				if (e[i] > 0) {
					if (maxSupply < e[i]) {
						maxSupply = e[i];
						k = i;
					}
				}
			}
			if (maxSupply == 0)
				break;
			delta = maxSupply;

			// delta= e[k];
			// if (-e[l]<delta) delta= e[k];
			// find delta (minimum on the path from k to l)
			int l = computeShortestPath(d, prev, k, rCostForward, rCostCapBackward, e);
			int to = l;
			do {
				int from = prev[to];	// assert (from != to);
				int itccb = 0;	// residual
				while ((itccb < rCostCapBackward[from].size()) && (rCostCapBackward[from].get(itccb).to != to)) {
					itccb++;
				}
				if (itccb < rCostCapBackward[from].size()) {
					if (rCostCapBackward[from].get(itccb).residual_capacity < delta)
						delta = rCostCapBackward[from].get(itccb).residual_capacity;
				}
				to = from;
			} while (to != k);

			// augment delta flow from k to l (backwards actually...)
			to = l;
			do {
				int from = prev[to];	// assert (from != to);
				// TODO - might do here O(n) can be done in O(1)
				int itx = 0;
				while (x[from].get(itx).to != to) {
					itx++;
				}
				x[from].get(itx).flow += delta;

				// update residual for backward edges
				int itccb = 0;
				while ((itccb < rCostCapBackward[to].size()) && (rCostCapBackward[to].get(itccb).to != from)) {
					itccb++;
				}
				if (itccb < rCostCapBackward[to].size()) {
					rCostCapBackward[to].get(itccb).residual_capacity += delta;
				}
				itccb = 0;
				while ((itccb < rCostCapBackward[from].size()) && (rCostCapBackward[from].get(itccb).to != to)) {
					itccb++;
				}
				if (itccb < rCostCapBackward[from].size()) {
					rCostCapBackward[from].get(itccb).residual_capacity -= delta;
				}

				// update e
				e[to] = e[to] + delta;
				e[from] = e[from] - delta;

				to = from;
			} while (to != k);
		}

		// compute distance from x
		long dist = 0;
		for (int from = 0; from < numNodes; from++) {
			for (EdgeWithFlow it : x[from]) {
				dist += (it.cost * it.flow);
			}
		}
		return dist;
	}

	// --------------------------------------------------------------------------------------

	private int computeShortestPath(long[] d, int[] prev, int from, 
			List<EdgeReducedForward>[] costForward, List<EdgeReducedBackward>[] costBackward, long[] e) {

		// Making heap (all inf except 0, so we are saving comparisons...)
		final List<Edge> Q = new ArrayList<Edge>(numNodes);
		for (int i = 0; i < numNodes; i++) {
			Q.add(new Edge());
		}

		Q.get(0).to = from;
		Q.get(0).cost = 0;
		nodesToQ[from] = 0;

		int j = 1;
		// TODO: both of these into a function?
		for (int i = 0; i < from; i++) {
			Q.get(j).to = i;
			Q.get(j).cost = Long.MAX_VALUE;
			nodesToQ[i] = j;
			j++;
		}

		for (int i = from + 1; i < numNodes; i++) {
			Q.get(j).to = i;
			Q.get(j).cost = Long.MAX_VALUE;
			nodesToQ[i] = j;
			j++;
		}

		final boolean[] finalNodesFlg = new boolean[numNodes];	// finalNodesFlg[i] = false;

		int l = 0;
		do {
			int u = Q.get(0).to;
			d[u] = Q.get(0).cost; // final distance
			finalNodesFlg[u] = true;
			if (e[u] < 0) {
				l = u;
				break;
			}

			heapRemoveFirst(Q, nodesToQ);

			// neighbors of u
			for (EdgeReducedForward it : costForward[u]) {
//				assert (it.cost >= 0);
				long alt = d[u] + it.cost;
				int v = it.to;
				if ((nodesToQ[v] < Q.size()) && (alt < Q.get(nodesToQ[v]).cost)) {
					heapDecreaseKey(Q, nodesToQ, v, alt);
					prev[v] = u;
				}
			}
			
			for (EdgeReducedBackward it : costBackward[u]) {
				if (it.residual_capacity > 0) {
//					assert (it.cost >= 0);
					long alt = d[u] + it.cost;
					int v = it.to;
					if ((nodesToQ[v] < Q.size()) && (alt < Q.get(nodesToQ[v]).cost)) {
						heapDecreaseKey(Q, nodesToQ, v, alt);
						prev[v] = u;
					}
				}
			}

		} while (!Q.isEmpty()); //(Q.size() > 0);

		for (int _from = 0; _from < numNodes; _from++) {
			for (EdgeReducedForward it : costForward[_from]) {
				if (finalNodesFlg[_from]) {
					it.cost += d[_from] - d[l];
				}
				if (finalNodesFlg[it.to]) {
					it.cost -= d[it.to] - d[l];
				}
			}
		}

		// reduced costs and capacity for backward edges
		// (c[j,i]-pi[j]+pi[i])
		for (int _from = 0; _from < numNodes; _from++) {
			for (EdgeReducedBackward it : costBackward[_from]) {
				if (finalNodesFlg[_from]) {
					it.cost += d[_from] - d[l];
				}
				if (finalNodesFlg[it.to]) {
					it.cost -= d[it.to] - d[l];
				}
			}
		}

		return l;
	}

	// --------------------------------------------------------------------------------------

	private void heapDecreaseKey(List<Edge> Q, int[] nodes_to_Q, int v, long alt) {
		int i = nodes_to_Q[v];
		Q.get(i).cost = alt;
		while (i > 0 && Q.get(getParent(i)).cost > Q.get(i).cost) {
			swapHeap(Q, nodes_to_Q, i, getParent(i));
			i = getParent(i);
		}
	}

	private void heapRemoveFirst(List<Edge> Q, int[] nodes_to_Q) {
		swapHeap(Q, nodes_to_Q, 0, Q.size() - 1);
		Q.remove(Q.size() - 1);
		heapify(Q, nodes_to_Q, 0);
	}	

	private void heapify(List<Edge> Q, int[] nodes_to_Q, int start) {
		int i = start;
		while (true) {
			// TODO: change to loop
			final int le = getLeft(i);
			final int re = getRight(i);
			int smallest = -1;
			if ((le < Q.size()) && (Q.get(le).cost < Q.get(i).cost)) {
				smallest = le;
			} else {
				smallest = i;
			}
			if ((re < Q.size()) && (Q.get(re).cost < Q.get(smallest).cost)) {
				smallest = re;
			}

			if (smallest == i)
				break;

			swapHeap(Q, nodes_to_Q, i, smallest);
			i = smallest;
		}
	}

	private void swapHeap(List<Edge> Q, int[] nodesToQ, int i, int j) {
		Edge tmp = Q.get(i);
		Q.set(i, Q.get(j));
		Q.set(j, tmp);
		nodesToQ[Q.get(j).to] = j;
		nodesToQ[Q.get(i).to] = i;
	}

	private int getLeft(int i) {
		return 2 * (i + 1) - 1;
	}

	private int getRight(int i) {
		return 2 * (i + 1); // 2 * (i + 1) + 1 - 1
	}

	private int getParent(int i) {
		return (i - 1) / 2;
	}
}