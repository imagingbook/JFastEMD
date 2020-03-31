package imagingbook.jfastemd;

/**
 * @author Telmo Menezes (telmo@telmomenezes.com)
 *
 */
public class Test {

    static double getValue(double[] map, int x, int y, int bins) {
        return map[(y * bins) + x];
    }

    static Signature getSignature(double[] map, int bins) {
        // find number of entries in the sparse matrix
        int n = 0;
        for (int x = 0; x < bins; x++) {
            for (int y = 0; y < bins; y++) {
                if (getValue(map, x, y, bins) > 0) {
                    n++;
                }
            }
        }
        
        // compute features and weights
        Feature2D[] features = new Feature2D[n];
        double[] weights = new double[n];
        int i = 0;
        for (int x = 0; x < bins; x++) {
            for (int y = 0; y < bins; y++) {
                double val = getValue(map, x, y, bins);
                if (val > 0) {
                    Feature2D f = new Feature2D(x, y);
                    features[i] = f;
                    weights[i] = val;
                    i++;
                }
            }
        }

        Signature signature = new Signature();
        signature.setNumberOfFeatures(n);
        signature.setFeatures(features);
        signature.setWeights(weights);

        return signature;
    }

    static double emdDist(double[] map1, double[] map2, int bins) {
        Signature sig1 = getSignature(map1, bins);
        Signature sig2 = getSignature(map2, bins);

        double dist = new JFastEMD(sig1, sig2, -1).getDistance();

        return dist;
    }
    
    // ----------------------------------------------------------------------------------------------------------
    
	static double[] a0 = { 1.0, 0.0, 0.0, 0.0 };
	static double[] a1 = { 0.0, 1.0, 0.0, 0.0 };
	static double[] a2 = { 0.0, 1.0, 1.0, 0.0 };
	static double[] b0 = { 1.0, 0.31350830458876927, 0.475451529763324, 0.710099174235318, 0.8180547959863713,
			0.8501705482451378, 0.7091117393023645, 0.3421407576224318, 0.0, 0.0, 0.8648715755286225, 0.0, 0.0,
			0.2296406823738588, 0.32854154105225764, 0.41240916326716803, 0.2556550109727834, 0.0, 0.0, 0.0,
			0.9688367289608703, 0.0, 0.15675415229438464, 0.2556550109727834, 0.4473841617389759, 0.5493536070363795,
			0.4274423997306564, 0.0, 0.0, 0.0, 0.9897035895139639, 0.0, 0.0, 0.3659766656118061, 0.49450429339199403,
			0.6613006436919866, 0.5263432584090552, 0.0, 0.0, 0.0, 0.9501732176368471, 0.0989008586783988,
			0.0989008586783988, 0.3956034347135952, 0.5865550152810103, 0.710099174235318, 0.6355631023347673,
			0.2556550109727834, 0.0, 0.0, 0.8570091629463844, 0.0, 0.0, 0.4042535848157063, 0.5841965520250411,
			0.7250339741455023, 0.6900589756736945, 0.3421407576224318, 0.0, 0.0, 0.6999031392570128, 0.0,
			0.15675415229438464, 0.2776498124065264, 0.5399424749792293, 0.6234493289150748, 0.6160355170421022,
			0.31350830458876927, 0.0, 0.0, 0.434403964700911, 0.0, 0.0, 0.0, 0.3863948346682434, 0.434403964700911,
			0.2776498124065264, 0.0989008586783988, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	static double[] b1 = { 0.7690558279672058, 0.5431277903401206, 0.7400048870218827, 0.9431329651179853, 1.0,
			0.8608371807886335, 0.6014104886149799, 0.0, 0.0, 0.0, 0.3620851935600804, 0.0, 0.0, 0.0, 0.0,
			0.1810425967800402, 0.1810425967800402, 0.0, 0.0, 0.0, 0.8073470804011144, 0.2869457269295445,
			0.3620851935600804, 0.3620851935600804, 0.6014104886149799, 0.5082508251725361, 0.0, 0.0, 0.0, 0.0,
			0.8407357836698793, 0.1810425967800402, 0.5082508251725361, 0.6892934219525761, 0.7824530853950201,
			0.626304483621074, 0.42036789183493967, 0.0, 0.0, 0.0, 0.8189574032199598, 0.42036789183493967,
			0.649030920489625, 0.7241703871201608, 0.7241703871201608, 0.8300735172696652, 0.2869457269295445, 0.0, 0.0,
			0.0, 0.8703360187326165, 0.42036789183493967, 0.5082508251725361, 0.7073136187644842, 0.921047483801923,
			0.7549340506391291, 0.626304483621074, 0.0, 0.0, 0.0, 0.573891453859089, 0.3620851935600804,
			0.46798832370958465, 0.5082508251725361, 0.3620851935600804, 0.5431277903401206, 0.1810425967800402, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	static double[] b2 = { 0.8430654934475055, 0.6239426073894736, 0.7730847346217581, 0.9313448972142149, 1.0,
			0.8783820744571478, 0.5344908771065096, 0.1686130986895011, 0.0, 0.0, 0.6419699112188565, 0.0,
			0.5058392960685033, 0.604471635932257, 0.5833054849228712, 0.5058392960685033, 0.1686130986895011, 0.0, 0.0,
			0.0, 0.3372261973790022, 0.1686130986895011, 0.0, 0.5058392960685033, 0.43585853724275586,
			0.43585853724275586, 0.2672454385532548, 0.0, 0.0, 0.0, 0.7287336883921703, 0.43585853724275586,
			0.1686130986895011, 0.5833054849228712, 0.5601205897026693, 0.47335681252935546, 0.3372261973790022, 0.0,
			0.0, 0.0, 0.7730847346217581, 0.43585853724275586, 0.5058392960685033, 0.5344908771065096,
			0.6587529295664228, 0.6239426073894736, 0.1686130986895011, 0.0, 0.0, 0.0, 0.7830149820263362,
			0.3915074910131681, 0.43585853724275586, 0.7162562210501102, 0.6587529295664228, 0.6891997754414121,
			0.2672454385532548, 0.1686130986895011, 0.0, 0.0, 0.5344908771065096, 0.0, 0.3372261973790022,
			0.2672454385532548, 0.3915074910131681, 0.43585853724275586, 0.1686130986895011, 0.0, 0.0, 0.0,
			0.1686130986895011, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    public static void main(String[] args) {
        System.out.println("test 1: " + emdDist(a0, a0, 2) + " [expected: 0.0]");
        System.out.println("test 1: " + emdDist(a0, a1, 2) + " [expected: 1.0]");
        System.out.println("test 2: " + emdDist(a0, a2, 2) + " [expected: 2.0]");
        System.out.println("test 3: " + emdDist(b0, b1, 10) + " [expected: 19.1921]");
        System.out.println("test 4: " + emdDist(b0, b2, 10) + " [expected: 25.7637]");
    }
}

/*
test 1: 0.0 [expected: 0.0]
test 1: 1.0 [expected: 1.0]
test 2: 2.0 [expected: 2.0]
test 3: 19.19207181614855 [expected: 19.1921]
test 4: 25.763736318380168 [expected: 25.7637]
*/