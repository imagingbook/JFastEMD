package imagingbook.jfastemd;

/**
 * @author Telmo Menezes (telmo@telmomenezes.com)
 *
 */

class Edges {

	static class Edge {
		int to;
		long cost;

		Edge(int to, long cost) {
			this.to = to;
			this.cost = cost;
		}
	}

// -----------------------------------------------------------

	static class Edge0 extends Edge {
		long flow;

		Edge0(int to, long cost, long flow) {
			super(to, cost);
			this.flow = flow;
		}
	}

//-----------------------------------------------------------

	static class Edge1 extends Edge {

		Edge1(int to, long reduced_cost) {
			super(to, reduced_cost);
		}
	}

//-----------------------------------------------------------

	static class Edge2 extends Edge1 {
		long residual_capacity;

		Edge2(int to, long reduced_cost, long residual_capacity) {
			super(to, reduced_cost);
			this.residual_capacity = residual_capacity;
		}
	}

//-----------------------------------------------------------

	static class Edge3 extends Edge {
		Edge3() {
			super(0, 0);
		}

		Edge3(int to, long dist) {
			super(to, dist);
		}
	}

//-----------------------------------------------------------

}
