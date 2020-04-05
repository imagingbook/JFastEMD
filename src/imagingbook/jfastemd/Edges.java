package imagingbook.jfastemd;

/**
 * @author Telmo Menezes (telmo@telmomenezes.com)
 * @author Wilbur (refactored)
 */

abstract class Edges {

	static class Edge {
		int to;
		long cost;
		
		Edge() {
			this(0, 0);
		}

		Edge(int to, long cost) {
			this.to = to;
			this.cost = cost;
		}
	}

// -----------------------------------------------------------

	static class EdgeWithFlow extends Edge {
		long flow;

		EdgeWithFlow(int to, long cost, long flow) {
			super(to, cost);
			this.flow = flow;
		}
	}

//-----------------------------------------------------------

	static class EdgeReducedForward extends Edge {

		EdgeReducedForward(int to, long reduced_cost) {
			super(to, reduced_cost);
		}
	}

//-----------------------------------------------------------

	static class EdgeReducedBackward extends EdgeReducedForward {
		long residual_capacity;

		EdgeReducedBackward(int to, long reduced_cost, long residual_capacity) {
			super(to, reduced_cost);
			this.residual_capacity = residual_capacity;
		}
	}

//-----------------------------------------------------------

}
