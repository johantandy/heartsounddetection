package featureExtraction;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;

import javax.management.monitor.GaugeMonitorMBean;
import javax.sound.sampled.UnsupportedAudioFileException;
import javax.swing.JFrame;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

import util.Database;

public class GTCC {
	private static final double MAX_16_BIT = Short.MAX_VALUE;
	FFT fft = new FFT();
	private static int dataPerFrame;
	static ArrayList<Double> gain;
	static ArrayList<Double> erb_point = new ArrayList<>();

	// DC REMOVAL - 1st Step GTCC
	public double[] DCRemoval(double[] data) {
		double[] result = new double[data.length];
		double avarage = 0;
		for (int i = 0; i < data.length; i++) {
			avarage = avarage + data[i];
		}
		avarage = avarage / data.length;
		// System.out.println(avarage);
		for (int i = 0; i < data.length; i++) {
			result[i] = data[i] - avarage;
			// System.out.println("y["+i+"]="+data[i]+"-("+avarage+")="+result[i]);
			// System.out.println(result[i]);
		}

		return result;
	}

	// PREEMPHASIZE 2nd Step GTCC
	public double[] PreEmphasize(double[] data, double alpha) {
		double[] result = new double[data.length];
		result[0] = data[0];
		for (int i = 1; i < data.length; i++) {
			result[i] = data[i] - alpha * data[i - 1];
			// System.out.println(result[i]);
		}
		return result;
	}

	// FRAME BLOCKING 3nd Step GTCC
	public static ArrayList<double[]> FrameBlocking(double[] data, int frameSize, int overlap) {
		int length = data.length;
		ArrayList<double[]> result = new ArrayList<double[]>();
		int index = 0;
		int bound = length - overlap;
		while (index < bound) {

			double[] frameData = new double[frameSize];
			for (int i = 0; i < frameSize; i++) {
				int frameNumber = index + i;
				if (frameNumber >= length) {
					frameData[i] = 0;
				} else {
					frameData[i] = data[frameNumber];
				}
			}

			index = index + overlap;
			// System.out.println(frameData[1]);
			result.add(frameData);
		}

		/*
		 * for (int i = 0; i < result.size(); i++) { //System.out.println("split"); for
		 * (int j = 0; j < result.get(i).length; j++) {
		 * System.out.println(result.get(i)[j]); }
		 * 
		 * }
		 */

		return result;
	}

	// WINDOWING WITH HAMMING WINDOW 4nd Step GTCC
	public double[] Windowing(double[] data, int frameSize) {
		int length = data.length;
		double[] result = new double[length];
		double divider = frameSize - 1;
		double pi2 = 2 * Math.PI;
		for (int i = 0; i < length; i++) { // System.out.println(divider);
											// System.out.println((0.54-0.46*Math.cos(pi2*i/divider)));
			result[i] = data[i] * (0.54 - 0.46 * Math.cos(pi2 * i / divider));
		}
		return result;
	}

	// GAMMATONE FILTER BANK 5nd Step GTCC
	public static ArrayList<ArrayList<Double>> GammatoneFilterBank(double fs) {
		int low_freq = 100;
		int high_freq = 44100 / 4;
		int numfbank = 64;
		double ear_q = 9.26449; // # Glasberg and Moore Parameters
		double min_bw = 24.7;
		double T = 1 / fs;
		erb_point = new ArrayList<>();
		ArrayList<Double> erb = new ArrayList<>();
		ArrayList<Double> b = new ArrayList<>();
		ArrayList<Double> arg = new ArrayList<>();
		ArrayList<Complex> vec = new ArrayList<>();
		ArrayList<Double> b1 = new ArrayList<>();
		ArrayList<Double> b2 = new ArrayList<>();
		ArrayList<Double> k11 = new ArrayList<>();
		ArrayList<Double> k12 = new ArrayList<>();
		ArrayList<Double> k13 = new ArrayList<>();
		ArrayList<Double> k14 = new ArrayList<>();
		ArrayList<Double> A11 = new ArrayList<>();
		ArrayList<Double> A12 = new ArrayList<>();
		ArrayList<Double> A13 = new ArrayList<>();
		ArrayList<Double> A14 = new ArrayList<>();
		ArrayList<Double> common = new ArrayList<>();
		ArrayList<Complex> gain_arg_temp = new ArrayList<>();
		ArrayList<Complex> gain_arg_temp1 = new ArrayList<>();
		ArrayList<Complex> gain_arg = new ArrayList<>();
		ArrayList<Complex> c31 = new ArrayList<>();
		ArrayList<Complex> c32 = new ArrayList<>();
		ArrayList<Complex> c33 = new ArrayList<>();
		ArrayList<Complex> c34 = new ArrayList<>();
		ArrayList<Complex> c35 = new ArrayList<>();
		ArrayList<Complex> c36 = new ArrayList<>();
		ArrayList<Complex> c37 = new ArrayList<>();
		ArrayList<Double> gtResponse = new ArrayList<>();
		gain = new ArrayList<>();
		ArrayList<Integer> allfilts = new ArrayList<>();
		ArrayList<ArrayList<Double>> fcoefs = new ArrayList<>();

		Complex cdua = new Complex(0, 2);
		Complex csatu = new Complex(0, 1);
		double A0 = T;
		double A2 = 0;
		double B0 = 1;
		double rt_pos = Math.sqrt(3 + Math.pow(2, 1.5));
		// System.out.println(rt_pos);
		double rt_neg = Math.sqrt(3 - Math.pow(2, 1.5));

		// System.out.println(-ear_q*min_bw+Math.exp((1/numfbank)*(-Math.log(high_freq+ear_q*min_bw)+Math.log(low_freq+ear_q*min_bw)))*(high_freq+ear_q*min_bw));

		double width = 1.0;
		int order = 1;
		for (int i = 1; i <= numfbank; i++) {
			erb_point.add(-ear_q * min_bw + Math
					.exp(((double) i / numfbank)
							* (-Math.log(high_freq + ear_q * min_bw) + Math.log(low_freq + ear_q * min_bw)))
					* (high_freq + ear_q * min_bw));
			// System.out.println(erb_point.get(i-1));
			// System.out.println((double)i/numfbank);
			erb.add(Math.pow(width * (Math.pow((erb_point.get(i - 1) / ear_q), order) + Math.pow(min_bw, 1)),
					1 / order));
			// System.out.println(erb.get(i-1));
			b.add(1.019 * 2 * Math.PI * erb.get(i - 1));
			// System.out.println(b.get(i-1));
			// arg.add(2*T*Math.exp(-2*Math.PI*b.get(i)*T));
			arg.add((2 * erb_point.get(i - 1) * Math.PI * T));
			// System.out.println(arg.get(i-1));
			Complex c1 = new Complex(arg.get(i - 1), 0);
			vec.add(cdua.times(c1).exp());
			// System.out.println(vec.get(i-1));
			b1.add(-2 * Math.cos(arg.get(i - 1)) / Math.exp(b.get(i - 1) * T));
			// System.out.println(b1.get(i-1));
			b2.add(Math.exp(-2 * b.get(i - 1) * T));
			// System.out.println(b2.get(i-1));
			common.add(-T * Math.exp(-(b.get(i - 1) * T)));
			// System.out.println(common.get(i-1));

			k11.add(Math.cos(arg.get(i - 1)) + rt_pos * Math.sin(arg.get(i - 1)));
			// System.out.println(k11.get(i-1));
			k12.add(Math.cos(arg.get(i - 1)) - rt_pos * Math.sin(arg.get(i - 1)));
			// System.out.println(k12.get(i-1));
			k13.add(Math.cos(arg.get(i - 1)) + rt_neg * Math.sin(arg.get(i - 1)));
			// System.out.println(k13.get(i-1));
			k14.add(Math.cos(arg.get(i - 1)) - rt_neg * Math.sin(arg.get(i - 1)));
			// System.out.println(k14.get(i-1));

			A11.add(common.get(i - 1) * k11.get(i - 1));
			A12.add(common.get(i - 1) * k12.get(i - 1));
			A13.add(common.get(i - 1) * k13.get(i - 1));
			A14.add(common.get(i - 1) * k14.get(i - 1));

			Complex c = new Complex(arg.get(i - 1), 0);
			Complex c10 = new Complex(b.get(i - 1) * T, 0);

			gain_arg_temp.add(c);
			gain_arg_temp1.add(c10);

			gain_arg.add(csatu.times(gain_arg_temp.get(i - 1)).minus(gain_arg_temp1.get(i - 1)).exp());
			// System.out.println(gain_arg.get(i-1));
			Complex c3 = new Complex(k11.get(i - 1), 0);
			Complex c4 = new Complex(k12.get(i - 1), 0);
			Complex c5 = new Complex(k13.get(i - 1), 0);
			Complex c6 = new Complex(k14.get(i - 1), 0);
			Complex c7 = new Complex(T * Math.exp(b.get(i - 1) * T), 0);
			Complex c8 = new Complex(-1 / Math.exp(b.get(i - 1) * T) + 1, 0);
			Complex c9 = new Complex(1 - Math.exp(b.get(i - 1) * T), 0);

			c31.add(c3);
			c32.add(c4);
			c33.add(c5);
			c34.add(c6);
			c35.add(c7);
			c36.add(c8);
			c37.add(c9);

			gain.add(((vec.get(i - 1).minus(gain_arg.get(i - 1).times(c31.get(i - 1))))
					.times((vec.get(i - 1).minus(gain_arg.get(i - 1).times(c32.get(i - 1)))))
					.times((vec.get(i - 1).minus(gain_arg.get(i - 1).times(c33.get(i - 1)))))
					.times((vec.get(i - 1).minus(gain_arg.get(i - 1).times(c34.get(i - 1))))))
							.times((c35.get(i - 1))
									.divides((c36.get(i - 1).plus((vec.get(i - 1).times(c37.get(i - 1)))))).power(4))
							.abs());

			// System.out.println("split");

			allfilts.add(Collections.frequency(erb_point, erb_point.get(i - 1)));

			ArrayList<Double> fcoefstemp = new ArrayList<>();

			fcoefstemp.add(A0 * allfilts.get(i - 1));
			fcoefstemp.add(A11.get(i - 1));
			fcoefstemp.add(A12.get(i - 1));
			fcoefstemp.add(A13.get(i - 1));
			fcoefstemp.add(A14.get(i - 1));
			fcoefstemp.add(A2 * allfilts.get(i - 1));
			fcoefstemp.add(B0 * allfilts.get(i - 1));
			fcoefstemp.add(b1.get(i - 1));
			fcoefstemp.add(b2.get(i - 1));
			fcoefstemp.add(gain.get(i - 1));

			fcoefs.add(fcoefstemp);

		}
		// counting
		// System.out.println(gain.size());
		/*
		 * double sumOfX = 0; for (int j = 0; j < 1; j++) {
		 * //System.out.println(source.get(j)); for (int j2 = 0; j2 < gain.size(); j2++)
		 * {
		 * 
		 * 
		 * //System.out.println(source.size()); sumOfX+=(source.get(j)*gain.get(j2)); }
		 * //System.out.println(sumOfX); gtResponse.add(sumOfX); sumOfX = 0;
		 * //System.out.println(gtResponse.size()); //System.out.println(gtResponse);
		 * //System.out.println("space"); }
		 */

		/*
		 * for (int j = 0; j < fcoefs.size(); j++) { for (int j2 = 0; j2 <
		 * fcoefs.get(0).size(); j2++) { System.out.println(fcoefs.get(j).get(j2));
		 * 
		 * } System.out.println("space");
		 * 
		 * }
		 */

		// System.out.println(erb_point.get(63));
		/*
		 * for (int j = 0; j < gain.size(); j++) { System.out.println(gain.get(j)); }
		 */

		/*
		 * for (int j = 0; j < fcoefs.size(); j++) { for (int j2 = 0; j2 <
		 * fcoefs.get(0).size(); j2++) { System.out.println(fcoefs.get(j).get(j2));
		 * 
		 * } System.out.println("space");
		 * 
		 * }
		 */
		return fcoefs;

	}

	public ArrayList<Double> countXm(ArrayList<Double> frequencyValue, ArrayList<ArrayList<Double>> GTFIlterBank) {
		double sumOfX = 0.;
		ArrayList<Double> xm = new ArrayList<>();
		for (int j2 = 0; j2 < frequencyValue.size() / 2 - 1; j2++) {
			for (int j = 0; j < GTFIlterBank.size(); j++) {
				for (int i = 0; i < GTFIlterBank.get(j).size(); i++) {

					sumOfX += (Math.pow(frequencyValue.get(j2), 2) * GTFIlterBank.get(j).get(i));
				}
			}
			// System.out.println(source.get(j));

			xm.add(sumOfX);
			sumOfX = 0.;
			// System.out.println(sumOfX);
			// System.out.println(gtResponse.size());
			// System.out.println(gtResponse);
			// System.out.println("space");
		}
		return xm;
	}

	public ArrayList<Double> EqualLoudnessCurve(ArrayList<Double> f, ArrayList<Double> erb) {
		ArrayList<Double> e = new ArrayList<>();
		ArrayList<Double> eql = new ArrayList<>();
		for (int j = 0; j < erb.size(); j++) {
			double wtemp = 2 * Math.PI * (double) erb.get(j);
			e.add(((Math.pow(wtemp, 2) + 56.8 * Math.pow(10, 6) * Math.pow(wtemp, 4))
					/ Math.pow((Math.pow(wtemp, 2) + 6.3 * Math.pow(10, 6)), 2)
					* (Math.pow(wtemp, 2) + 0.38 * Math.pow(10, 9)) * (Math.pow(wtemp, 6) + 9.58 * Math.pow(10, 26))));
			// System.out.println(e.get(j));
			eql.add(e.get(j) * f.get(j));
			// System.out.println(eql.get(j));
		}
		return eql;

	}

	public ArrayList<Double> LogaritmhmicCompression(ArrayList<Double> xme) {
		ArrayList<Double> xln = new ArrayList<>();
		for (int i = 0; i < xme.size(); i++) {
			xln.add(Math.log10(xme.get(i)));
			// System.out.println(xln.get(i));
		}

		return xln;

	}

	public double[] GetFeatureVector(double[] data, double alpha, int size, int overlap) {
		double[] dcRemoval = DCRemoval(data);
		double[] preEmphasized = PreEmphasize(dcRemoval, alpha);
		// ArrayList<double[]> frame = FrameBlocking(preEmphasized,frameSize);
		ArrayList<double[]> frame = FrameBlocking(preEmphasized, size, overlap);

		ArrayList<double[]> result = new ArrayList<double[]>();
		double[][] gtccFeature = new double[frame.size()][];
		double[][] framedSignal = new double[frame.size()][];
		ArrayList<ArrayList<Double>> GTFIlterBank = GammatoneFilterBank(16000);
		ArrayList<ArrayList<Double>> xm = new ArrayList<>();
		ArrayList<ArrayList<Double>> equalloudness = new ArrayList<>();
		ArrayList<ArrayList<Double>> logarithmicCompression = new ArrayList<>();
		ArrayList<ArrayList<Double>> discreteCosineTransform = new ArrayList<>();
		System.out.println(frame.size());
		DCT dct = new DCT(12, 64);
		for (int i = 0; i < frame.size(); i++) {

			double[] window = Windowing(frame.get(i), size);
			framedSignal[i] = window;

			/*
			 * for (int j = 0; j < window.length; j++) { System.out.println(window[j]); }
			 */
			Complex[] signal = new Complex[window.length];
			// System.out.println(window.length);
			for (int x = 0; x < window.length; x++) {
				Complex c = new Complex(window[x], 0);
				signal[x] = c;
			}

			Complex[] hasil = FFT.fft1D(signal);
			/*
			 * for (int j = 0; j < hasil.length; j++) { System.out.println(hasil[j]); }
			 * 
			 * System.out.println("split");
			 */

			ArrayList<Double> frequencyValue = new ArrayList<>();

			for (int k = 0; k < window.length; k++) {
				double nilaiReal = hasil[k].re();
				double nilaiImag = hasil[k].im();

				frequencyValue.add(Math.sqrt(nilaiReal * nilaiReal + nilaiImag * nilaiImag));
				// System.out.println(frequencyValue.si);
				// System.out.println(frequencyValue.get(k));
			}

			xm.add(countXm(frequencyValue, GTFIlterBank));

			equalloudness.add(EqualLoudnessCurve(xm.get(i), erb_point));
			logarithmicCompression.add(LogaritmhmicCompression(equalloudness.get(i)));
			discreteCosineTransform.add(dct.performDCT(logarithmicCompression.get(i)));
			// System.out.println(discreteCosineTransform.get(i));
			// System.out.println(" ");

		}

		/*System.out.println("reduce dimension");
		double[][] datafeature = new double[discreteCosineTransform.size()][discreteCosineTransform.get(0).size()];

		for (int i = 0; i < datafeature.length; i++) {
			for (int j = 0; j < datafeature[0].length; j++) {
				datafeature[i][j] = discreteCosineTransform.get(i).get(j);
			}
		}
		Database db = new Database();
		int lengthlast=0;
		try {
			lengthlast = db.selectlength("lengtha");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int multip = 0, itera = 0;
		int tempitera=0;
		System.out.println(lengthlast);
		while (multip < lengthlast) {
			tempitera=itera;
			multip = 12 * itera;
			itera++;
		}
		System.out.println("itera "+tempitera);
		itera=tempitera;

		SingularValueDecomposition svd = new SingularValueDecomposition(MatrixUtils.createRealMatrix(datafeature).transpose());
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();
		double[][] news = new double[12][itera];
		double[][] newv = new double[itera][itera];
		for (int j = 0; j < news.length; j++) {
			for (int j2 = 0; j2 < itera; j2++) {
				news[j][j2] = sar[j][j2];
			}
		}
		
		for (int i = 0; i < itera; i++) {
			for (int j = 0; j < itera; j++) {
				newv[i][j]=vtar[i][j];
			}
		}
		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);

		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		double[][] featurereduced = featuresvd.transpose().getData();
		double[] featuresvdreduc = new double[featurereduced.length*featurereduced[0].length];
		int iterators = 0;
		for (int i = 0; i < featurereduced.length; i++) {
			for (int j = 0; j < featurereduced[0].length; j++) {
				featuresvdreduc[iterators] = featurereduced[i][j];
				iterators++;
			}
		}*/
		
		double[] gtccFeatures = new double[discreteCosineTransform.get(0).size() * discreteCosineTransform.size()];

		int iterator = 0;
		for (int i = 0; i < frame.size(); i++) {
			for (int j = 0; j < discreteCosineTransform.get(i).size(); j++) {
				// System.out.println();
				gtccFeatures[iterator] = discreteCosineTransform.get(i).get(j);
				iterator++;
			}
		}

		//System.out.println(gtccFeatures.length);
		//System.out.println("length svdfeature "+featuresvdreduc.length );
		return gtccFeatures;
	}

	public static void write(PrintStream printStream, String... values) throws Exception {

		try {
			for (String value : values) {
				printStream.println(value);
			}
			printStream.flush();
		} catch (Exception e) {
			// handling exception
		}
	}

	public static void main(String[] args) {
		// GammatoneFilterBank(null, 16000);
		String Sound1 = "C:\\Users\\extre\\Desktop\\heart audio\\Atraining_murmur\\Atraining_murmur\\201108222223.wav";
		double[] data = StdAudio.read(Sound1);
		double[] datasample = new double[301];
		for (int i = 0; i <= 300; i++) {
			datasample[i] = data[i];
			// System.out.println(datasample[i]);
		}

		/*
		 * PrintStream printStream; try { printStream = new PrintStream(new
		 * File("output3.txt")); // hook for closing the stream
		 * Runtime.getRuntime().addShutdownHook(new Thread(printStream::close)); //
		 * writing write(printStream,Arrays.toString(datasample)); } catch
		 * (FileNotFoundException e) { // TODO Auto-generated catch block
		 * e.printStackTrace(); } catch (Exception e) { // TODO Auto-generated catch
		 * block e.printStackTrace(); }
		 */
		/*
		 * Yin yin = new Yin(16000); try { yin.main(Sound1); } catch
		 * (UnsupportedAudioFileException | IOException e) { // TODO Auto-generated
		 * catch block e.printStackTrace(); }
		 */

		// System.out.println("oke");
		GTCC gtcc = new GTCC();

		// gtcc.GammatoneFilterBank(16000);
		double[] dataa = gtcc.GetFeatureVector(data, 0.9, 4000, 16000);

		/*
		 * ArrayList<Float> pitch1 = new ArrayList<Float>(); ArrayList<Float> timesuji =
		 * new ArrayList<Float>();
		 */
		/*
		 * Yin yin = new Yin(16000); try { yin.main(Sound1); } catch
		 * (UnsupportedAudioFileException e) { // TODO Auto-generated catch block
		 * e.printStackTrace(); } catch (IOException e) { // TODO Auto-generated catch
		 * block e.printStackTrace(); } pitch1=yin.getPitchs(); timesuji=yin.getTimes();
		 * System.out.println(pitch1.size()); System.out.println("MULAI"); for (int i =
		 * 0; i < pitch1.size(); i++) { if (pitch1.get(i)!=-1.0) {
		 * 
		 * System.out.println(timesuji.get(i)+"  "+pitch1.get(i)); } }
		 * System.out.println("selesai");
		 */
		/*
		 * for (int i = 0; i < dataa.length; i++) { System.out.println(dataa[i]); }
		 */

		// GMM gmm = new GMM();
		// gmm.main(gtcc.GetFeatureVector(data, 0.9, 4000, 16000));

		/*
		 * for (int i = 1; i <=195; i++) { Syst em.out.println(i); }
		 */
		/*
		 * XYSeries xys = new XYSeries("plot audio"); for (int i = 0; i < data.length;
		 * i++) { xys.add(i, data[i]); } XYSeriesCollection dataset = new
		 * XYSeriesCollection(); dataset.addSeries(xys); XyChart.dataSet=dataset;
		 * XyChart chart = new XyChart("plot audio", "plot audio"); chart.pack();
		 * RefineryUtilities.centerFrameOnScreen(chart);
		 * chart.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		 * chart.setVisible(true);
		 */
		/*
		 * GTCC gtcc = new GTCC(); gtcc.GetFeatureVector(data, 0.9, 4000, 16000);
		 */

		/*
		 * Complex c = new Complex(0, 2); Complex c1 = new Complex(1.1, 0);
		 * 
		 * System.out.println((c.times(c1).exp()));
		 */

		// ArrayList<double[]> hasilframeblocking = FrameBlocking(data, 4000,1600);

	}

}
