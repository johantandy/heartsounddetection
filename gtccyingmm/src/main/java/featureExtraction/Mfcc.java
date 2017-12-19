package featureExtraction;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.JFrame;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;



public class Mfcc {
	FFT fft = new FFT();
	DCTs dct;
	public int noOfFrames;
	public int numCepstra=12;
	private static int dataPerFrame;
	static String fileName;
	public  Mfcc()
	{
	}
	
	public double[] CeptralLiftering(double[] data)
	{
		int n = data.length;
		double[] result = new double[data.length];
		for (int i=0;i<data.length;i++)
		{
			result[i] = data[i] *(1+(n/2)*Math.sin(i*Math.PI/n));
		}
		return result;
	}
	
	public double[] DCRemoval (double[] data)
	{
		double[] result = new double[data.length];
		double avarage=0;
		for (int i=0;i<data.length;i++)
		{
			avarage = avarage + data[i];
		}
		avarage = avarage/data.length;
		for (int i=0;i<data.length;i++)
		{
			result[i] = data[i]-avarage;
		}
		return result;
	}
	
	/*public double[] DCT(double[] source)
	{
		int data = source.length;
		double[] result = new double[data];
		for (int i=0;i<data;i++)
		{
			double a = i*Math.PI/data;
			for (int j=0;j<data;j++)
			{
				result[i] = result[i] + (source[j]*Math.cos(a*(j-0.5)));
			}
		}
		return result;
	}*/
	
	public ArrayList<double[]> FrameBlocking (double[] data, int frameSize, int overlap)
	{
		int length = data.length;
		ArrayList<double[]> result = new ArrayList<double[]>();
		int index = 0;
		int bound = length-overlap;
		//XYSeriesCollection dataset = new XYSeriesCollection();
		//int va = index;
		while (index<bound)
		{
			double[] frameData = new double[frameSize];
			for (int i=0;i<frameSize;i++)
			{
				int frameNumber = index+i;
				if (frameNumber>=length)
				{
					frameData[i]=0;
				}
				else
				{
					frameData[i] = data[frameNumber];
				}
			}
			//va = index;
			index = index + overlap;
			result.add(frameData);
			/*XYSeries frames = new XYSeries ("frame "+result.size());
			for (int i=0;i<frameSize;i++)
			{
				frames.add(va, frameData[i]);
				va++;
			}
			dataset.addSeries(frames);*/
		}
		/*XyChart.dataSet = dataset;
		XyChart chart2 = new XyChart("frameChart "+fileName,"framing tetap Chart "+fileName);
		chart2.pack( );          
	      RefineryUtilities.centerFrameOnScreen( chart2 );
	      chart2.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	      chart2.setVisible( true );*/
		return result;
	}
	
	//experiment
	public static ArrayList<double[]> FrameBlocking(double[] data, int totalFrame)
	{
		int dataLength = data.length;
		dataPerFrame = (int) (Math.ceil((double) dataLength / (double) totalFrame)*1.25);
		int overlap = (int) Math.ceil(0.25 * (double)dataPerFrame);
		/*System.out.println("parameter");
		System.out.println("panjang data = "+dataLength);
		System.out.println("panjang frame = "+dataPerFrame);
		System.out.println("overlap = "+overlap);*/
		ArrayList<double[]> result = new ArrayList<double[]>();
		int pointer = 0;
		//XYSeriesCollection dataset = new XYSeriesCollection();
		//int te = pointer;
		while (result.size() < totalFrame)
		{
			double[] temp = new double[dataPerFrame];
			for (int i=0;i<dataPerFrame;i++)
			{
				if (pointer < data.length)
				{
					temp[i] = data[pointer];
				}
				else
				{
					temp[i] = 0;
				}
				pointer++;
			}
			result.add(temp);
			/*XYSeries frames = new XYSeries ("frame "+result.size());
			te = pointer;
			for (int i=0;i<dataPerFrame;i++)
			{
				frames.add(te, temp[i]);
				te++;
			}
			dataset.addSeries(frames);*/
			pointer = pointer - overlap;
		}
		/*XyChart.dataSet = dataset;
		XyChart chart2 = new XyChart("frameChart "+fileName,"framing fleksibel Chart "+fileName);
		chart2.pack( );          
	      RefineryUtilities.centerFrameOnScreen( chart2 );
	      chart2.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	      chart2.setVisible( true );*/
		return result;
	}
	
	public double[] MelFrequencyWrapping(double[] source, double fs)
	{
		double melLow = 2595*Math.log10(1+(300/700));
		double melHigh = 2595*Math.log10(1+(8000/700));
		int p = 26;
		int data = source.length;
		double divider = (melHigh-melLow)/(p+1);
		double[] fb = new double[p+2];
		for (int i=0;i<=p+1;i++)
		{
			double mel = melLow+i*divider;
			double inverseMel = 700*(Math.pow(10, mel/2595)-1);
			fb[i] = data/fs*inverseMel;
		}
		double[] x = new double[p];
		for (int i = 1; i<= p;i++)
		{
			for (int k=0;k<data;k++)
			{
				if (k<fb[i-1])
				{
				}
				else if ((fb[i-1]<k) && (k<=fb[i]))
				{
					x[i-1] = x[i-1] + (source[k]*((k-fb[i-1])/(fb[i]-fb[i-1])));
				}
				else if ((fb[i]<k) && (k<=fb[i+1]))
				{
					x[i-1] = x[i-1] + (source[k]*((fb[i+1]-k)/(fb[i+1]-fb[i])));
				}
			}
		}
		return x;
	}
	
	public double[] PreEmphasize (double[] data, double alpha)
	{
		double[] result = new double[data.length];
		result[0] = data[0];
		for (int i=1;i<data.length;i++)
		{
			result[i] = data[i] - alpha*data[i-1];
		}
		return result;
	}
	
	public double[] Windowing (double[] data, int frameSize)
	{
		int length = data.length;
		double[] result = new double[length];
		double divider = frameSize-1;
		double pi2 = 2*Math.PI;
		for (int i=0;i<length;i++)
		{
			result[i] = data[i] * (0.54-0.46*Math.cos(pi2*i/divider));
		}
		return result;
	}
	
	public double[][] doCepstralMeanNormalization (double[][] mfccFeature)
	{
		double sum;
		double mean;
		double[][] mCeps = new double[mfccFeature.length][numCepstra];
		for (int i=0;i<numCepstra;i++)
		{
			sum = 0.0;
			for (int j=0;j<mfccFeature.length;j++)
			{
				sum = sum + mfccFeature[j][i];
			}
			mean = sum / noOfFrames;
			for (int j=0;j<mfccFeature.length;j++)
			{
				mCeps[j][i] = mfccFeature[j][i] - mean;
			}
		}
		return mCeps;
	}
	
	public double[] GetFeatureVector (double[] data, double alpha, int size, int overlap)
	{
		double[] dcRemoval_1 = DCRemoval(data);
		//System.out.println("Selesai DCRemoval");
		double[] preEmphasized = PreEmphasize(dcRemoval_1,alpha);
		//System.out.println("Selesai PReEmphasize");
		ArrayList<double[]> frame = FrameBlocking(preEmphasized,size,overlap);
		//System.out.println("Selesai FrameBlocking");
		ArrayList<double[]> result = new ArrayList<double[]>();
		dct = new DCTs(numCepstra,26);
		double[][] mfccFeature = new double[frame.size()][];
		double[][] framedSignal = new double[frame.size()][];
		//System.out.println ("Frame size = " + frame.size());
		System.out.println(frame.size());
		for (int i=0;i<frame.size();i++)
		{
			double[] window = Windowing (frame.get(i),size);
			//System.out.println("Selesai windowing");
			framedSignal[i] = window;
			/*if (i==0)
			{
				System.out.println ("frame 0 window");
				System.out.print("[");
				for (int j=0;j<window.length;j++)
				{
					System.out.print(window[j]+" , ");
				}
				System.out.println ("]");
			}*/
			Complex[] signal = new Complex[window.length];
			for (int x=0;x<window.length;x++)
			{
				Complex c = new Complex (window[x],0);
				signal[x] = c;
			}
			Complex[] hasil = fft.fft1D(signal);
			//System.out.println("Selesai fft");
			double[] frequencyValue = new double[window.length];
			for (int k = 0; k < window.length; k++) {
				double nilaiReal = hasil[k].re();
				double nilaiImag = hasil[k].im();
				frequencyValue[k] = Math.sqrt(nilaiReal * nilaiReal + nilaiImag * nilaiImag);
			}
			double[] melFrequency = MelFrequencyWrapping(frequencyValue,16000);
			//System.out.println("Selesai MelWrapping");
			double[] cepstrum = dct.performDCT(melFrequency);
			//System.out.println("Selesai DCT");
			double[] cepstral = CeptralLiftering(cepstrum);
			//System.out.println("Selesai Liftering");
			mfccFeature[i] = cepstral;
		}
		
		noOfFrames = frame.size();
		double[][] normalCeps = doCepstralMeanNormalization(mfccFeature);
		double[][] featureVector = new double[frame.size()][3*numCepstra+3];
		for (int i=0;i<frame.size();i++)
		{
			for (int j = 0; j < numCepstra; j++) {
				featureVector[i][j] = mfccFeature[i][j];
			}
		}
		
		ArrayList<Double> store = new ArrayList<>();
		int iter=0;
		for (int i = 0; i < mfccFeature.length; i++) {
			for (int j = 0; j < mfccFeature[i].length; j++) {
				store.add(mfccFeature[i][j]);
				
				}
		}
		
		
		
		double []features=new double[store.size()];
		for (int i = 0; i < store.size(); i++) {
			features[i]=store.get(i);
		}
		System.out.println(Arrays.deepToString(mfccFeature));
		return features;
	}
	
	
	public double[][] GetFeatureVector (double[] data, double alpha, int frameSize)
	{
		double[] dcRemoval_1 = DCRemoval(data);
		//System.out.println("Selesai DCRemoval");
		double[] preEmphasized = PreEmphasize(dcRemoval_1,alpha);
		//System.out.println("Selesai PReEmphasize");
		ArrayList<double[]> frame = FrameBlocking(preEmphasized,frameSize);
		//System.out.println("Selesai FrameBlocking");
		ArrayList<double[]> result = new ArrayList<double[]>();
		dct = new DCTs(numCepstra,26);
		double[][] mfccFeature = new double[frame.size()][];
		double[][] framedSignal = new double[frame.size()][];
		//System.out.println ("Frame size = " + frame.size());
		for (int i=0;i<frame.size();i++)
		{
			double[] window = Windowing (frame.get(i),dataPerFrame);
			//System.out.println("Selesai windowing");
			framedSignal[i] = window;
			/*if (i==0)
			{
				System.out.println ("frame 0 window");
				System.out.print("[");
				for (int j=0;j<window.length;j++)
				{
					System.out.print(window[j]+" , ");
				}
				System.out.println ("]");
			}*/
			Complex[] signal = new Complex[window.length];
			for (int x=0;x<window.length;x++)
			{
				Complex c = new Complex (window[x],0);
				signal[x] = c;
			}
			Complex[] hasil = fft.fft1D(signal);
			//System.out.println("Selesai fft");
			double[] frequencyValue = new double[window.length];
			for (int k = 0; k < window.length; k++) {
				double nilaiReal = hasil[k].re();
				double nilaiImag = hasil[k].im();
				frequencyValue[k] = Math.sqrt(nilaiReal * nilaiReal + nilaiImag * nilaiImag);
			}
			double[] melFrequency = MelFrequencyWrapping(frequencyValue,16000);
			//System.out.println("Selesai MelWrapping");
			double[] cepstrum = dct.performDCT(melFrequency);
			//System.out.println("Selesai DCT");
			double[] cepstral = CeptralLiftering(cepstrum);
			//System.out.println("Selesai Liftering");
			mfccFeature[i] = cepstral;
		}

		noOfFrames = frame.size();
		double[][] normalCeps = doCepstralMeanNormalization(mfccFeature);
	
		double[][] featureVector = new double[frame.size()][3*numCepstra+3];
		for (int i=0;i<frame.size();i++)
		{
 			for (int j = 0; j < numCepstra; j++) {
				featureVector[i][j] = mfccFeature[i][j];
				//System.out.println(featureVector[i][j]);
			}
			
		}
		
		double[] features = new double[featureVector.length*featureVector[0].length];
		for (int i = 0; i < featureVector.length; i++) {
			for (int j = 0; j < featureVector[i].length; j++) {
				features[j]=featureVector[i][j];
			}
		}
		
		System.out.println(Arrays.toString(features));
		return featureVector;
	}
	
	
	public static void main(String[] args) {
		/*String sound1 = "C:\\Users\\extre\\Desktop\\heart audio\\UNSUPERVISED\\201012172012.wav";
		String sound2 = "C:\\Users\\extre\\Desktop\\heart audio\\UNSUPERVISED\\201101051104.wav";
		double[] sound1data = StdAudio.read(sound1);
		double[] sound1data1 = StdAudio.read(sound2);
		Mfcc mfcc = new Mfcc();
		System.out.println(Arrays.deepToString(mfcc.GetFeatureVector(sound1data, 0.9, 4000, 1600)));
		mfcc = new Mfcc();
		System.out.println(Arrays.deepToString(mfcc.GetFeatureVector(sound1data1, 0.9, 4000, 1600)));*/
		
		String path = "C:\\Users\\extre\\Desktop\\heart audio\\UNSUPERVISED";
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		double[] duration = new double[directoryListing.length];
		double[][] fctore = new double[directoryListing.length][];
		int iter = 0;
		if (directoryListing != null) {
			for (File child : directoryListing) {
				double[] sound1data = StdAudio.read(child.getAbsolutePath());
				Mfcc mfcc = new Mfcc();
				//mfcc.GetFeatureVector(sound1data, 0.9, 4000, 1600);
				//System.out.println(Arrays.deepToString(mfcc.GetFeatureVector(sound1data, 0.9, 4000, 1600)));
		} 
	}
	}}
	/*public static void main (String[] args)
	{
		String sound1 = "C:\\Users\\hobert\\workspace\\dataLatihSuara\\siapa\\sample_2_siapa.wav";
		String sound2 = "C:\\Users\\hobert\\workspace\\data suara panjang\\panjang_siapa.wav";
		double[] sound1data = StdAudio.read(sound1);
		double[] sound2data = StdAudio.read(sound2);
		fileName = "sample_2_siapa.wav";
		ArrayList<double[]> frameSound1tetap = FrameBlocking(sound1data,4000,1600);
		ArrayList<double[]> frameSound1fleksibel = FrameBlocking (sound1data,13);
		fileName = "siapa panjang.wav";
		ArrayList<double[]> frameSound2tetap = FrameBlocking (sound2data,4000,1600);
		ArrayList<double[]> frameSound2fleksibel = FrameBlocking(sound2data,13);
		
	}*/
