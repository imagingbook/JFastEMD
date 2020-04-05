package imagingbook.jfastemd.utils;

import com.crtomirmajer.wmd4j.emd.EarthMovers;

import imagingbook.jfastemd.JFastEMD2;

/**
 * Created by Majer on 22. 09. 2016.
 */
public class JFastEmdVsEarthMoversPerformanceTestsMain {
    
    public static void main(String[] args) {
    	
    	
        
        int size = 20;
        double[] a = EarthMoversUtils.randomVector(size);
        double[] b = EarthMoversUtils.randomVector(size);
        double[][] m = EarthMoversUtils.matrix(a, b);
        
        int repeats = 10000; // 10000;
        
        //warm up
        System.out.println("warmup...");
//        EarthMovers earthMovers = new EarthMovers();
        for(int i = 0 ; i < repeats ; i++) {
        	EarthMoversUtils.jfastemd(a, b, m, 1);
        	new EarthMovers().distance(a, b, m, 1);
        }
        
        System.out.println("Results:");
        System.out.format("telmomenezes-jfastemd:           %.3f\n", EarthMoversUtils.jfastemd(a, b, m, 1));
        System.out.format("crtomirmajer optimized-jfastemd: %.3f\n", new EarthMovers().distance(a, b, m, 1));
        System.out.format("wilburs jfastemd:                %.3f\n", new JFastEMD2(a, b, m, 1).getDistance());
        
        
        // ----------------------------------------------------------------------
//        System.out.println("starting telmomenezes-jfastemd...");
        long start = System.nanoTime();
        for(int i = 0 ; i < repeats ; i++) {
            EarthMoversUtils.jfastemd(a, b, m, 1);
        }
        double time_telmomenezes = (System.nanoTime() - start) / 1E9;
        System.out.format("telmomenezes: %8.3fs = %5.1f%%\n", time_telmomenezes, 100.0);
        // ----------------------------------------------------------------------
//        System.out.println("starting crtomirmajer optimized-jfastemd...");
        start = System.nanoTime();
        for(int i = 0 ; i < repeats ; i++) {
        	new EarthMovers().distance(a, b, m, 1);
        }
        double time_crtomirmajer = (System.nanoTime() - start) / 1E9;
        double ratio_crtomirmajer = (double) time_crtomirmajer / time_telmomenezes;
        System.out.format("crtomirmajer: %8.3fs = %5.1f%%\n", time_crtomirmajer, ratio_crtomirmajer * 100);
        // ----------------------------------------------------------------------
//        System.out.println("starting wilburs jfastemd...");
        start = System.nanoTime();
        for(int i = 0 ; i < repeats ; i++) {
            new JFastEMD2(a, b, m, 1).getDistance();
        }
        double time_wilbur = (System.nanoTime() - start) / 1E9;
        double ratio_wilbur = (double) time_wilbur / time_telmomenezes;
        System.out.format("wilbur:       %8.3fs = %5.1f%%\n", time_wilbur, ratio_wilbur * 100);
        // ----------------------------------------------------------------------
        System.out.println("done.");
    }
}
