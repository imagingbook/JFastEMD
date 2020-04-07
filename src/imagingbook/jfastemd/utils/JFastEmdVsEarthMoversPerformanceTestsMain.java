package imagingbook.jfastemd.utils;

import java.util.Vector;

import com.crtomirmajer.wmd4j.emd.EarthMovers;
import com.telmomenezes.jfastemd.JFastEMD;

import imagingbook.jfastemd.JFastEMD2;

/**
 * Created by Majer on 22. 09. 2016.
 * Modified by Wilbur 2020
 */
public class JFastEmdVsEarthMoversPerformanceTestsMain {
    
    public static void main(String[] args) {
    	
        int size = 20;
        
        double[] P = EarthMoversUtils.randomVector(size);	// the source weights
        double[] Q = EarthMoversUtils.randomVector(size);	// the target weights
        double[][] C = EarthMoversUtils.matrix(P, Q);		// the cost (distance) matrix
        
        Vector<Double> vecP = EarthMoversUtils.makeVector(P);
        Vector<Double> vecQ = EarthMoversUtils.makeVector(Q);
        Vector<Vector<Double>> vecC = EarthMoversUtils.makeMatrix(C);
        
        int repeats = 10000; // 10000;
        
        System.out.println("size = " + size);
        System.out.println("a = " + P.length);
        System.out.println("b = " + Q.length);
        System.out.println("m = " + C.length + " x " + C[0].length);
        
        //warm up
        System.out.println("warmup..." + repeats);
//        EarthMovers earthMovers = new EarthMovers();
        for(int i = 0 ; i < repeats ; i++) {
        	JFastEMD.emdHat(vecP, vecQ, vecC, 1);
        	new EarthMovers().distance(P, Q, C, 1);
        	new JFastEMD2(P, Q, C, 1).getDistance();
        }
        
        System.out.println("Results:");
        System.out.format("telmomenezes-jfastemd:           %.3f\n", JFastEMD.emdHat(vecP, vecQ, vecC, 1));
        System.out.format("crtomirmajer optimized-jfastemd: %.3f\n", new EarthMovers().distance(P, Q, C, 1));
        System.out.format("wilburs jfastemd:                %.3f\n", new JFastEMD2(P, Q, C, 1).getDistance());
        
        
        // ----------------------------------------------------------------------
//        System.out.println("starting telmomenezes-jfastemd...");
        long start = System.nanoTime();
        for(int i = 0 ; i < repeats ; i++) {
        	JFastEMD.emdHat(vecP, vecQ, vecC, 1);
            //EarthMoversUtils.jfastemd(P, Q, C, 1);
        }
        double time_telmomenezes = (System.nanoTime() - start) / 1E9;
        System.out.format("telmomenezes: %8.3fs = %5.1f%%\n", time_telmomenezes, 100.0);
        // ----------------------------------------------------------------------
//        System.out.println("starting crtomirmajer optimized-jfastemd...");
        start = System.nanoTime();
        for(int i = 0 ; i < repeats ; i++) {
        	new EarthMovers().distance(P, Q, C, 1);
        }
        double time_crtomirmajer = (System.nanoTime() - start) / 1E9;
        double ratio_crtomirmajer = (double) time_crtomirmajer / time_telmomenezes;
        System.out.format("crtomirmajer: %8.3fs = %5.1f%%\n", time_crtomirmajer, ratio_crtomirmajer * 100);
        // ----------------------------------------------------------------------
//        System.out.println("starting wilburs jfastemd...");
        start = System.nanoTime();
        for(int i = 0 ; i < repeats ; i++) {
            new JFastEMD2(P, Q, C, 1).getDistance();
        }
        double time_wilbur = (System.nanoTime() - start) / 1E9;
        double ratio_wilbur = (double) time_wilbur / time_telmomenezes;
        System.out.format("wilbur:       %8.3fs = %5.1f%%\n", time_wilbur, ratio_wilbur * 100);
        // ----------------------------------------------------------------------
        System.out.println("done.");
    }
}
/*
 * Results on i5 PC (for repeat = 10000):
telmomenezes:   15.567s = 100.0%
crtomirmajer:    7.592s =  48.8%
wilbur:          6.413s =  41.2%

 * Results on i7 Laptop (for repeat = 10000):
telmomenezes:   12.744s = 100.0%
crtomirmajer:    5.967s =  46.8%
wilbur:          5.019s =  39.4%
*/
