package imagingbook.emd_opencv;

import imagingbook.jfastemd.signatures.Signature;

/**
 * This is a Java port of emd.cpp (3/14/98) taken from the OpenCV
 * release 4.
 * 
 * @author WB
 *
 */
public class EmdOpencv {
	
	private static final int MAX_ITERATIONS = 500;
	private static final double CV_EMD_INF = 1e20;
	private static final double CV_EMD_EPS = 1e-5;
	
	// state variables:
	double max_cost;	// max value in cost matrix
	int ssize, dsize;
	double[][] cost;
	boolean[][] is_x;
	CvNode1D[] u, v;	// CvNode1D
	CvNode2D enter_x;	// CvNode2D
	double weight;
	private int[] idx1;
	private int[] idx2;
	private double[] s;
	private double[] d;
	private boolean is_used;
	private double[][] delta;
	private CvNode2D[] loop;
	private Object rows_x;
	private CvNode2D[] cols_x;
	
	
	
	public double EMD(Signature _signature1, Signature _signature2,
			int distType, double[][] _cost) {
		
		int n1 = _signature1.getNumberOfFeatures();
		int n2 = _signature2.getNumberOfFeatures();
		
		double[] signature1 = _signature1.getWeights();
		double[] signature2 = _signature2.getWeights();
		
		double[][] flow = new double[n1][n2];
		
		return cvCalcEMD2(signature1, signature2, distType, _cost, flow, Double.NaN);
	}
	
	// ----------------------------------------------------------------------------------

	private double cvCalcEMD2(double[] signature1, double[] signature2, int distType, double[][] _cost,
			double[][] flow, double lower_bound) {
		
		double emd = 0;
		double total_cost = 0;
		int result = 0;
	    double eps, min_delta;
	    
	    int dims = 0; // signature1->cols - 1; // dimensionality of signature weights?
	    int size1 = signature1.length;
	    int size2 = signature2.length;
	    
	    Object dist_func = null;
	    Object user_param = null;
	    
	    result = icvInitEMD(signature1, signature2, dims, dist_func, user_param, _cost, lower_bound);
	    
	    if( result > 0 && Double.isFinite(lower_bound) ) {	// wilbur: simplify!
	    	emd = lower_bound;
	        return emd;
	    }
	    
	    eps = CV_EMD_EPS * this.max_cost;
	    
	    /* if ssize = 1 or dsize = 1 then we are done, else ... */
	    if( this.ssize > 1 && this.dsize > 1 ) {

	        for(int itr = 1; itr < MAX_ITERATIONS; itr++ ) {
	        	/* find basic variables */
	            result = icvFindBasicVariables( this.cost, this.is_x, this.u, this.v, this.ssize, this.dsize );
	            
	            if( result < 0 )
	                break;
	            
	            /* check for optimality */
	            min_delta = icvIsOptimal( this.cost, this.is_x,
	            		this.u, this.v,
	            		this.ssize, this.dsize, this.enter_x );
	            
	            if( min_delta == CV_EMD_INF )
	                throw new RuntimeException("did not converge (1)");
	            
	            /* if no negative deltamin, we found the optimal solution */
	            if( min_delta >= -eps )
	                break;
	            
	            /* improve solution */
	            if(!icvNewSolution())
	            	throw new RuntimeException("did not converge (2)");
	        }
	    }
	    
	    /* compute the total flow */
//	    for( xp = state._x; xp < state.end_x; xp++ ) {
//	    	
//	    }
		
	    emd = total_cost / this.weight;
		return emd;
	}
	
	/************************************************************************************\
	*          initialize structure, allocate buffers and generate initial solution      *
	\************************************************************************************/

	private int icvInitEMD(double[] signature1, double[] signature2, int dims, Object dist_func, Object user_param,
			double[][] cost, double lower_bound) {
		
		int size1 = signature1.length;
		int size2 = signature2.length;

		double s_sum = 0, d_sum = 0, diff;
	    //int i, j;
	    int ssize = 0, dsize = 0;
	    boolean equal_sums = true;
	    int buffer_size;
	    double max_cost = 0;
	    
	    this.idx1 = new int[size1 + 1];
	    this.idx2 = new int[size2 + 1];
	    
	    this.s = new double[size1 + 1];
	    this.d = new double[size2 + 1];
	    
	    /* sum up the supply and demand */
	    for(int i = 0; i < size1; i++ ) {
	        double weight = signature1[i];	// [i * (dims + 1)]
	        if (weight > 0) {
	            s_sum += weight;
	            this.s[ssize] = weight;
	            this.idx1[ssize++] = i;
	        }
	        else if( weight < 0 )
	            throw new IllegalArgumentException("signature1 must not contain negative weights");
	    }
	    
	    for(int i = 0; i < size2; i++ ) {
	        double weight = signature2[i];	// [i * (dims + 1)]
	        if (weight > 0) {
	            d_sum += weight;
	            this.d[dsize] = weight;
	            this.idx2[dsize++] = i;
	        }
	        else if( weight < 0 )
	        	throw new IllegalArgumentException("signature2 must not contain negative weights");
	    }
	    
	    if (ssize == 0)
	    	throw new IllegalArgumentException("signature1 must contain at least one non-zero value");
	    if (dsize == 0)
	    	throw new IllegalArgumentException("signature2 must contain at least one non-zero value");
	    
	    /* if supply different than the demand, add a zero-cost dummy cluster */
	    diff = s_sum - d_sum;
	    if (Math.abs(diff) >= CV_EMD_EPS * s_sum) {
	    	equal_sums = false;
	    	if (diff < 0) {
	            this.s[ssize] = -diff;
	            this.idx1[ssize++] = -1;
	        }
	    	else {
	            this.d[dsize] = diff;
	            this.idx2[dsize++] = -1;
	        }
	    }
	    
	    this.ssize = ssize;
	    this.dsize = dsize;
	    this.weight = s_sum > d_sum ? s_sum : d_sum;	// max?
	    
	    if (Double.isFinite(lower_bound) && equal_sums) {    /* check lower bound */
	    	int sz1 = size1 * (dims + 1);
	    	int sz2 = size2 * (dims + 1);
	    	double lb = 0;
	    	
	    	double[] xs = s;	// ??
	    	double[] xd = d;
	    	
			for (int j = 0; j < sz1; j += dims + 1) {
				double weight = signature1[j];
				for (int i = 0; i < dims; i++)
					xs[i] += signature1[j + i + 1] * weight; // ??
			}
			
			for (int j = 0; j < sz2; j += dims + 1) {
				double weight = signature2[j];
				for (int i = 0; i < dims; i++)
					xd[i] += signature2[j + i + 1] * weight;
			}
			
			lb = 0; // dist_func( xs, xd, user_param ) / state->weight;		// CHECK!!
			
			boolean ii = (lower_bound <= lb);	// ??
			lower_bound = lb;	// needs to be returned??
			if (ii) 
				return 1;
	    }
	    
	    /* assign pointers */
	    this.is_used = false;	// points to same address as delta matrix!
	    /* init delta matrix */
	    this.delta = new double[ssize][dsize];		// ssize * sizeof( float * );
	    
	    this.loop = new CvNode2D[ssize + dsize + 1];	// CHECK line 502
	    
	    /* init cost matrix */
	    this.cost = new double[ssize][dsize];
	    
	    /* compute the distance matrix */
		for (int i = 0; i < ssize; i++) {
			int ci = this.idx1[i];
			if (ci >= 0) {
				for(int j = 0; j < dsize; j++) {
					int cj = this.idx2[j];
					if (cj < 0)
	                    this.cost[i][j] = 0;
					else {
						double val;
						if (dist_func != null) {
							val = 1;	// ??
//									dist_func( signature1 + ci * (dims + 1) + 1,
//											signature2 + cj * (dims + 1) + 1,
//											user_param );
						}
						else {
							val = cost[ci][cj];
						}
						this.cost[i][j] = val;
						if (max_cost < val)
	                        max_cost = val;
					}
				}
			}

			else {
				for(int j = 0; j < dsize; j++ )	
					this.cost[i][j] = 0;		// why zero cost? needed at all?
			}
		}
		
		this.max_cost = max_cost;
		
		this.rows_x = new CvNode2D[ssize];
		this.cols_x = new CvNode2D[dsize];
		this.u = new CvNode1D[ssize];
		this.v = new CvNode1D[dsize];
		
		/* init is_x matrix */
		this.is_x = new boolean[ssize][dsize];
		
		icvRussel();
		
		// this.enter_x = (state->end_x)++; // what the f...?
	    
		return 0;
	}

	// ----------------------------------------------------------------------------------
	
	private int icvFindBasicVariables(double[][] cost, boolean[][] is_x, 
			CvNode1D[] u, CvNode1D[] v, int ssize, int dsize) {
		int i, j;
	    int u_cfound, v_cfound;
	    CvNode1D u0_head, u1_head; // *cur_u, *prev_u;
	    CvNode1D v0_head, v1_head; // *cur_v, *prev_v;
	    boolean found;
	    
	    // TODO: UNFINISHED!!
		return 0;
	}

	// ----------------------------------------------------------------------------------

	private double icvIsOptimal(double[][] cost, boolean[][] is_x, Object u2, Object v2, int ssize2, int dsize2,
			Object enter_x2) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	// ----------------------------------------------------------------------------------
	
	private boolean icvNewSolution() {
		// TODO Auto-generated method stub
		return false;
	}

	// ----------------------------------------------------------------------------------
	
	private void icvRussel() {
		// TODO Auto-generated method stub
		
	}


}
