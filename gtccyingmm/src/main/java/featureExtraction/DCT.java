/*
  Please feel free to use/modify this class. 
  If you give me credit by keeping this information or
  by sending me an email before using it or by reporting bugs , i will be happy.
  Email : gtiwari333@gmail.com,
  Blog : http://ganeshtiwaridotcomdotnp.blogspot.com/ 
 */
package featureExtraction;

import java.math.BigDecimal;
import java.util.ArrayList;

/**
 * performs Inverser Fourier Transform <br>
 * we used Dct because there is only real coeffs
 * 
 * @author Ganesh Tiwari
 */
public class DCT {

	/**
	 * number of GTCC coeffs
	 */
	int numCepstra;
	/**
	 * number of GTCC Filters
	 */
	int M;

	/**
	 * @param len
	 *            length of array, i.e., number of features
	 * @param M
	 *            numbe of Mel Filters
	 * @return
	 */
	public DCT(int numCepstra, int M) {
		this.numCepstra = numCepstra;
		this.M = M;
	}

	public ArrayList<Double> performDCT(ArrayList<Double> y) {
		double cepctemp =0;
		ArrayList<Double> cepc = new ArrayList<>();
		// perform DCT
		cepctemp += y.get(0) * Math.sqrt(1./64.);;
		for (int n = 1; n <= numCepstra; n++) {
			for (int i = 1; i <= M; i++) {
				cepctemp += y.get(i-1) * Math.cos((Math.PI * n * (i - 0.5))/M);
			}
			
			cepctemp=cepctemp*Math.sqrt(0.03125);
			//System.out.println((double)Math.sqrt((sqrt)));
			cepc.add(cepctemp);
		
		}
		return cepc;
	}
}
