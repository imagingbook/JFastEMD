package imagingbook.jfastemd;

import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import imagingbook.jfastemd.Edges.Edge;
import imagingbook.jfastemd.Edges.Edge0;
import imagingbook.jfastemd.Edges.Edge1;
import imagingbook.jfastemd.Edges.Edge2;
import imagingbook.jfastemd.Edges.Edge3;



class MinCostFlow {
	private final int numNodes;
	private final int[] nodesToQ;
	
	MinCostFlow(int numNodes) {
		this.numNodes = numNodes;
		this.nodesToQ = new int[numNodes];
	}

	// e - supply(positive) and demand(negative).
	// c[i] - edges that goes from node i. first is the second node
	// x - the flow is returned in it
	long compute(long[] e, List<Edge>[] c, List<Edge0>[] x) {
		if (e.length != numNodes || c.length != numNodes || c.length != numNodes) {
			throw new IllegalArgumentException("Either e, c or x is not of expected length " + numNodes);
		}

		// init flow
		for (int from = 0; from < numNodes; from++) {
			for (Edge it : c[from]) {
				x[from].add(new Edge0(it.to, it.cost, 0));
				x[it.to].add(new Edge0(from, -it.cost, 0));
			}
		}

		// reduced costs for forward edges (c[i,j]-pi[i]+pi[j])
		// Note that for forward edges the residual capacity is infinity
		@SuppressWarnings("unchecked")
		List<Edge1>[] rCostForward = new List[numNodes];
		for (int i = 0; i < numNodes; i++) {
			rCostForward[i] = new LinkedList<Edge1>();
		}
		for (int from = 0; from < numNodes; from++) {
			for (Edge it : c[from]) {
				rCostForward[from].add(new Edge1(it.to, it.cost));
			}
		}

		// reduced costs and capacity for backward edges
		// (c[j,i]-pi[j]+pi[i])
		// Since the flow at the beginning is 0, the residual capacity is
		// also zero
		@SuppressWarnings("unchecked")
		List<Edge2>[] rCostCapBackward = new List[numNodes];
		for (int i = 0; i < numNodes; i++) {
			rCostCapBackward[i] = new LinkedList<Edge2>();
		}
		for (int from = 0; from < numNodes; from++) {
			for (Edge it : c[from]) {
				rCostCapBackward[it.to].add(new Edge2(from, -it.cost, 0));
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

		long[] d = new long[numNodes];
		int[] prev = new int[numNodes];

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
				int from = prev[to];
				assert (from != to);
				// residual
				int itccb = 0;
				while ((itccb < rCostCapBackward[from].size())
						&& (rCostCapBackward[from].get(itccb).to != to)) {
					itccb++;
				}
				if (itccb < rCostCapBackward[from].size()) {
					if (rCostCapBackward[from].get(itccb).residual_capacity < delta)
						delta = rCostCapBackward[from].get(itccb).residual_capacity;
				}
				to = from;
			} while (to != k);

			// augment delta flow from k to l (backwards actually...)
//			to = l[0];
			to = l;
			do {
				int from = prev[to];
				assert (from != to);

				// TODO - might do here O(n) can be done in O(1)
				int itx = 0;
				while (x[from].get(itx).to != to) {
					itx++;
				}
				x[from].get(itx).flow += delta;

				// update residual for backward edges
				int itccb = 0;
				while ((itccb < rCostCapBackward[to].size())
						&& (rCostCapBackward[to].get(itccb).to != from)) {
					itccb++;
				}
				if (itccb < rCostCapBackward[to].size()) {
					rCostCapBackward[to].get(itccb).residual_capacity += delta;
				}
				itccb = 0;
				while ((itccb < rCostCapBackward[from].size())
						&& (rCostCapBackward[from].get(itccb).to != to)) {
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
			for (Edge0 it : x[from]) {
				dist += (it.cost * it.flow);
			}
		}
		return dist;
	}

	// --------------------------------------------------------------------------------------

	int computeShortestPath(long[] d, int[] prev, int from, 
			List<Edge1>[] costForward, List<Edge2>[] costBackward, long[] e) {

		// Making heap (all inf except 0, so we are saving comparisons...)
		Vector<Edge3> Q = new Vector<Edge3>();
		for (int i = 0; i < numNodes; i++) {
			Q.add(new Edge3());
		}

		Q.get(0).to = from;
		nodesToQ[from] = 0;
		Q.get(0).cost = 0;

		int j = 1;
		// TODO: both of these into a function?
		for (int i = 0; i < from; i++) {
			Q.get(j).to = i;
			nodesToQ[i] = j;
			Q.get(j).cost = Long.MAX_VALUE;
			j++;
		}

		for (int i = from + 1; i < numNodes; i++) {
			Q.get(j).to = i;
			nodesToQ[i] = j;
			Q.get(j).cost = Long.MAX_VALUE;
			j++;
		}

		Vector<Boolean> finalNodesFlg = new Vector<Boolean>();
		for (int i = 0; i < numNodes; i++) {
			finalNodesFlg.add(false);
		}
		
		int l = 0;
		do {
			int u = Q.get(0).to;
			d[u] = Q.get(0).cost; // final distance
			finalNodesFlg.set(u, true);
			if (e[u] < 0) {
				l = u;
				break;
			}

			heapRemoveFirst(Q, nodesToQ);

			// neighbors of u
			for (Edge1 it : costForward[u]) {
				assert (it.cost >= 0);
				long alt = d[u] + it.cost;
				int v = it.to;
				if ((nodesToQ[v] < Q.size()) && (alt < Q.get(nodesToQ[v]).cost)) {
					heapDecreaseKey(Q, nodesToQ, v, alt);
					prev[v] = u;
				}
			}
			for (Edge2 it : costBackward[u]) {
				if (it.residual_capacity > 0) {
					assert (it.cost >= 0);
					long alt = d[u] + it.cost;
					int v = it.to;
					if ((nodesToQ[v] < Q.size()) && (alt < Q.get(nodesToQ[v]).cost)) {
						heapDecreaseKey(Q, nodesToQ, v, alt);
						prev[v] = u;
					}
				}
			}

		} while (Q.size() > 0);

		for (int _from = 0; _from < numNodes; _from++) {
			for (Edge1 it : costForward[_from]) {
				if (finalNodesFlg.get(_from)) {
					it.cost += d[_from] - d[l];
				}
				if (finalNodesFlg.get(it.to)) {
					it.cost -= d[it.to] - d[l];
				}
			}
		}

		// reduced costs and capacity for backward edges
		// (c[j,i]-pi[j]+pi[i])
		for (int _from = 0; _from < numNodes; _from++) {
			for (Edge2 it : costBackward[_from]) {
				if (finalNodesFlg.get(_from)) {
					it.cost += d[_from] - d[l];
				}
				if (finalNodesFlg.get(it.to)) {
					it.cost -= d[it.to] - d[l];
				}
			}
		}
		
		return l;
	}

	// --------------------------------------------------------------------------------------

	void heapDecreaseKey(Vector<Edge3> Q, int[] nodes_to_Q,
			int v, long alt) {
		int i = nodes_to_Q[v];
		Q.get(i).cost = alt;
		while (i > 0 && Q.get(PARENT(i)).cost > Q.get(i).cost) {
			swapHeap(Q, nodes_to_Q, i, PARENT(i));
			i = PARENT(i);
		}
	}

	void heapRemoveFirst(Vector<Edge3> Q, int[] nodes_to_Q) {
		swapHeap(Q, nodes_to_Q, 0, Q.size() - 1);
		Q.remove(Q.size() - 1);
		heapify(Q, nodes_to_Q, 0);
	}

	void heapify(Vector<Edge3> Q, int[] nodes_to_Q, int i) {
		do {
			// TODO: change to loop
			int l = LEFT(i);
			int r = RIGHT(i);
			int smallest;
			if ((l < Q.size()) && (Q.get(l).cost < Q.get(i).cost)) {
				smallest = l;
			} else {
				smallest = i;
			}
			if ((r < Q.size()) && (Q.get(r).cost < Q.get(smallest).cost)) {
				smallest = r;
			}

			if (smallest == i)
				return;

			swapHeap(Q, nodes_to_Q, i, smallest);
			i = smallest;

		} while (true);
	}

	void swapHeap(Vector<Edge3> Q, int[] nodesToQ, int i, int j) {
		Edge3 tmp = Q.get(i);
		Q.set(i, Q.get(j));
		Q.set(j, tmp);
		nodesToQ[Q.get(j).to] = j;
		nodesToQ[Q.get(i).to] = i;
	}

	int LEFT(int i) {
		return 2 * (i + 1) - 1;
	}

	int RIGHT(int i) {
		return 2 * (i + 1); // 2 * (i + 1) + 1 - 1
	}

	int PARENT(int i) {
		return (i - 1) / 2;
	}
}