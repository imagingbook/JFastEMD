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

	int numNodes;
	Vector<Integer> nodesToQ;
	
	@Deprecated
	long compute(long[] ea, List<Edge>[] ca, List<Edge0>[] xa) {
		
		Vector<Long> e = new Vector<>();
		for (long val : ea) {
			e.add(val);
		}
		
		Vector<List<Edge>> c = new Vector<List<Edge>>();
		for (List<Edge> le : ca) {
			c.add(le);
		}
		
		Vector<List<Edge0>> x = new Vector<List<Edge0>>();
		for (List<Edge0> le : xa) {
			x.add(le);
		}
		
//		return compute(e, c, x);
		return compute(ea, c, x);
	}
	
	

	// e - supply(positive) and demand(negative).
	// c[i] - edges that goes from node i. first is the second nod
	// x - the flow is returned in it
//	long compute(Vector<Long> e, Vector<List<Edge>> c, Vector<List<Edge0>> x) {
	long compute(long[] e, Vector<List<Edge>> c, Vector<List<Edge0>> x) {
//		assert (e.size() == c.size());
		assert (e.length == c.size());
		assert (x.size() == c.size());

//		numNodes = e.size();
		numNodes = e.length;
		nodesToQ = new Vector<Integer>();
		for (int i = 0; i < numNodes; i++) {
			nodesToQ.add(0);
		}

		// init flow
		for (int from = 0; from < numNodes; from++) {
			for (Edge it : c.get(from)) {
				x.get(from).add(new Edge0(it.to, it.cost, 0));
				x.get(it.to).add(new Edge0(from, -it.cost, 0));
			}
		}

		// reduced costs for forward edges (c[i,j]-pi[i]+pi[j])
		// Note that for forward edges the residual capacity is infinity
//		Vector<List<Edge1>> rCostForward = new Vector<List<Edge1>>();
		@SuppressWarnings("unchecked")
		List<Edge1>[] rCostForward = new List[numNodes];
		for (int i = 0; i < numNodes; i++) {
//			rCostForward.add(new LinkedList<Edge1>());
			rCostForward[i] = new LinkedList<Edge1>();
		}
		for (int from = 0; from < numNodes; from++) {
			for (Edge it : c.get(from)) {
//				rCostForward.get(from).add(new Edge1(it.to, it.cost));
				rCostForward[from].add(new Edge1(it));
			}
		}

		// reduced costs and capacity for backward edges
		// (c[j,i]-pi[j]+pi[i])
		// Since the flow at the beginning is 0, the residual capacity is
		// also zero
		Vector<List<Edge2>> rCostCapBackward = new Vector<List<Edge2>>();
		for (int i = 0; i < numNodes; i++) {
			rCostCapBackward.add(new LinkedList<Edge2>());
		}
		for (int from = 0; from < numNodes; from++) {
			for (Edge it : c.get(from)) {
				rCostCapBackward.get(it.to).add(
						new Edge2(from, -it.cost, 0));
			}
		}

		// Max supply TODO:demand?, given U?, optimization-> min out of
		// demand,supply
		long U = 0;
		for (int i = 0; i < numNodes; i++) {
//			if (e.get(i) > U) {
//				U = e.get(i);
//			}
			if (e[i] > U) {
				U = e[i];
			}
		}
//		long delta = (long) (Math.pow(2.0, Math.ceil(Math.log((double) (U)) / Math.log(2.0))));	// TODO: value never used

//		Vector<Long> d = new Vector<Long>();
//		Vector<Integer> prev = new Vector<Integer>();
		long[] d = new long[numNodes];		// TODO: scratch?
		int[] prev = new int[numNodes];
		
//		for (int i = 0; i < numNodes; i++) {
//			d.add(0l);
//			prev.add(0);
//		}
		
		long delta = 1;
		while (true) { // until we break when S or T is empty
			long maxSupply = 0;
			int k = 0;
			for (int i = 0; i < numNodes; i++) {
//				if (e.get(i) > 0) {
//					if (maxSupply < e.get(i)) {
//						maxSupply = e.get(i);
//						k = i;
//					}
//				}
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

			int[] l = new int[1];
			computeShortestPath(d, prev, k, rCostForward, rCostCapBackward, e, l);

			// find delta (minimum on the path from k to l)
			// delta= e[k];
			// if (-e[l]<delta) delta= e[k];
			int to = l[0];
			do {
//				int from = prev.get(to);
				int from = prev[to];
				assert (from != to);

				// residual
				int itccb = 0;
				while ((itccb < rCostCapBackward.get(from).size())
						&& (rCostCapBackward.get(from).get(itccb).to != to)) {
					itccb++;
				}
				if (itccb < rCostCapBackward.get(from).size()) {
					if (rCostCapBackward.get(from).get(itccb).residual_capacity < delta)
						delta = rCostCapBackward.get(from).get(itccb).residual_capacity;
				}

				to = from;
			} while (to != k);

			// augment delta flow from k to l (backwards actually...)
			to = l[0];
			do {
//				int from = prev.get(to);
				int from = prev[to];
				assert (from != to);

				// TODO - might do here O(n) can be done in O(1)
				int itx = 0;
				while (x.get(from).get(itx).to != to) {
					itx++;
				}
				x.get(from).get(itx).flow += delta;

				// update residual for backward edges
				int itccb = 0;
				while ((itccb < rCostCapBackward.get(to).size())
						&& (rCostCapBackward.get(to).get(itccb).to != from)) {
					itccb++;
				}
				if (itccb < rCostCapBackward.get(to).size()) {
					rCostCapBackward.get(to).get(itccb).residual_capacity += delta;
				}
				itccb = 0;
				while ((itccb < rCostCapBackward.get(from).size())
						&& (rCostCapBackward.get(from).get(itccb).to != to)) {
					itccb++;
				}
				if (itccb < rCostCapBackward.get(from).size()) {
					rCostCapBackward.get(from).get(itccb).residual_capacity -= delta;
				}

				// update e
//				e.set(to, e.get(to) + delta);
//				e.set(from, e.get(from) - delta);
				e[to] = e[to] + delta;
				e[from] = e[from] - delta;

				to = from;
			} while (to != k);
		}

		// compute distance from x
		long dist = 0;
		for (int from = 0; from < numNodes; from++) {
			for (Edge0 it : x.get(from)) {
				dist += (it.cost * it.flow);
			}
		}
		return dist;
	}
	
	// --------------------------------------------------------------------------------------

	void computeShortestPath(
			long[] d, // Vector<Long> d, // this should not be passed, since not used outside!
			int[] prev, // Vector<Integer> prev,
			int from, 
			List<Edge1>[] costForward,	// Vector<List<Edge1>> costForward,
			Vector<List<Edge2>> costBackward, 
			long[] e,  // Vector<Long> e, 
			int[] l) {
		// Making heap (all inf except 0, so we are saving comparisons...)
		Vector<Edge3> Q = new Vector<Edge3>();
		for (int i = 0; i < numNodes; i++) {
			Q.add(new Edge3());
		}

		Q.get(0).to = from;
		nodesToQ.set(from, 0);
		Q.get(0).cost = 0;

		int j = 1;
		// TODO: both of these into a function?
		for (int i = 0; i < from; ++i) {
			Q.get(j).to = i;
			nodesToQ.set(i, j);
			Q.get(j).cost = Long.MAX_VALUE;
			j++;
		}

		for (int i = from + 1; i < numNodes; i++) {
			Q.get(j).to = i;
			nodesToQ.set(i, j);
			Q.get(j).cost = Long.MAX_VALUE;
			j++;
		}

		Vector<Boolean> finalNodesFlg = new Vector<Boolean>();
		for (int i = 0; i < numNodes; i++) {
			finalNodesFlg.add(false);
		}
		do {
			int u = Q.get(0).to;

//			d.set(u, Q.get(0).cost); // final distance
			d[u] = Q.get(0).cost; // final distance
			finalNodesFlg.set(u, true);
			if (e[u] < 0) {	// if (e.get(u) < 0) {
				l[0] = u;
				break;
			}

			heapRemoveFirst(Q, nodesToQ);

			// neighbors of u
			for (Edge1 it : costForward[u]) {	// for (Edge1 it : costForward.get(u)) {
				assert (it.cost >= 0);
//				long alt = d.get(u) + it.cost;
				long alt = d[u] + it.cost;
				int v = it.to;
				if ((nodesToQ.get(v) < Q.size())
						&& (alt < Q.get(nodesToQ.get(v)).cost)) {
					heapDecreaseKey(Q, nodesToQ, v, alt);
//					prev.set(v, u);
					prev[v] = u;
				}
			}
			for (Edge2 it : costBackward.get(u)) {
				if (it.residual_capacity > 0) {
					assert (it.cost >= 0);
//					long alt = d.get(u) + it.cost;
					long alt = d[u] + it.cost;
					int v = it.to;
					if ((nodesToQ.get(v) < Q.size())
							&& (alt < Q.get(nodesToQ.get(v)).cost)) {
						heapDecreaseKey(Q, nodesToQ, v, alt);
//						prev.set(v, u);
						prev[v] = u;
					}
				}
			}

		} while (Q.size() > 0);

		for (int _from = 0; _from < numNodes; _from++) {
			for (Edge1 it : costForward[_from]) {	// for (Edge1 it : costForward.get(_from)) {
				if (finalNodesFlg.get(_from)) {
//					it.cost += d.get(_from) - d.get(l[0]);
					it.cost += d[_from] - d[l[0]];
				}
				if (finalNodesFlg.get(it.to)) {
//					it.cost -= d.get(it.to) - d.get(l[0]);
					it.cost -= d[it.to] - d[l[0]];
				}
			}
		}

		// reduced costs and capacity for backward edges
		// (c[j,i]-pi[j]+pi[i])
		for (int _from = 0; _from < numNodes; _from++) {
			for (Edge2 it : costBackward.get(_from)) {
				if (finalNodesFlg.get(_from)) {
//					it.cost += d.get(_from) - d.get(l[0]);
					it.cost += d[_from] - d[l[0]];
				}
				if (finalNodesFlg.get(it.to)) {
//					it.cost -= d.get(it.to) - d.get(l[0]);
					it.cost -= d[it.to] - d[l[0]];
				}
			}
		}
	}
	
	// --------------------------------------------------------------------------------------

	void heapDecreaseKey(Vector<Edge3> Q, Vector<Integer> nodes_to_Q,
			int v, long alt) {
		int i = nodes_to_Q.get(v);
		Q.get(i).cost = alt;
		while (i > 0 && Q.get(PARENT(i)).cost > Q.get(i).cost) {
			swapHeap(Q, nodes_to_Q, i, PARENT(i));
			i = PARENT(i);
		}
	}

	void heapRemoveFirst(Vector<Edge3> Q, Vector<Integer> nodes_to_Q) {
		swapHeap(Q, nodes_to_Q, 0, Q.size() - 1);
		Q.remove(Q.size() - 1);
		heapify(Q, nodes_to_Q, 0);
	}

	void heapify(Vector<Edge3> Q, Vector<Integer> nodes_to_Q, int i) {
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

	void swapHeap(Vector<Edge3> Q, Vector<Integer> nodesToQ, int i, int j) {
		Edge3 tmp = Q.get(i);
		Q.set(i, Q.get(j));
		Q.set(j, tmp);
		nodesToQ.set(Q.get(j).to, j);
		nodesToQ.set(Q.get(i).to, i);
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