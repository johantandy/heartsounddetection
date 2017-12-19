package gtccyingmm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.chart.Effect3D;

import com.github.habernal.confusionmatrix.ConfusionMatrix;
import com.google.gson.Gson;

import Jama.Matrix;
import featureExtraction.GTCC;
import featureExtraction.GTCCMatlab;
import featureExtraction.StdAudio;
import machineLearning.GaussianMixtureModel;
import machineLearning.gmm1;
import util.CSVHelper;
import util.Database;
import util.ObjectFeatureGtcc;
import util.ObjectParamGMM;
import yin.Yin;

public class Controller {
	/*
	 * String databasea = "gtccaset"; String databaseb = "gtccbset"; String
	 * databaseayin = "yinaset"; String databasebyin = "yinbset"; String
	 * databasegmma = "gtccgmmaset"; String databasegmmb = "gtccgmmbset"; String
	 * databasegmmyina = "yingmmaset"; String databasegmmyinb = "yingmmbset";
	 */

	String databasea = "gtccasetnormalless";
	String databaseb = "gtccbsetnormalless";
	String databaseayin = "yinasetnormalless";
	String databasebyin = "yinbsetnormalless";
	static String databasegmma = "gtccgmmasetnormalless";
	static String databasegmmb = "gtccgmmbsetnormalless";
	static String databasegmmyina = "yingmmasetnormalless";
	static String databasegmmyinb = "yingmmbsetnormalless";

	public void toJson(double [][][] cov,double [][] means , double[][] featurearray) {
		Map<String, double[][][]> jscov = new HashMap<>();
		jscov.put("cov", cov);
		Map<String, double[][]> jsmean = new HashMap<>();
		jsmean.put("mean",means);
		Map<String, double[][]> jsx = new HashMap<>();
		jsx.put("x", featurearray);
		Gson gson = new Gson();
		
		try {
			Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("cova.json"),"utf-8"));
			writer.write(gson.toJson(jscov));
			writer.close();
			Writer writer1 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("means.json"),"utf-8"));
			writer1.write(gson.toJson(jsmean));
			writer1.close();
			Writer writer2 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("x.json"),"utf-8"));
			writer2.write(gson.toJson(jsx));
			writer2.close();
		} catch (UnsupportedEncodingException e) 	{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void toJsonPredictTrain(ArrayList<Integer> px) {
		Map<String, ArrayList<Integer>> jspx = new HashMap<>();
		jspx.put("px", px);
		Gson gson = new Gson();
		
		try {
			Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("pxtrain.json"),"utf-8"));
			writer.write(gson.toJson(jspx));
			writer.close();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void toJsonPredict(ArrayList<Integer> px) {
		Map<String, ArrayList<Integer>> jspx = new HashMap<>();
		jspx.put("px", px);
		Gson gson = new Gson();
		
		try {
			Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("pxtest.json"),"utf-8"));
			writer.write(gson.toJson(jspx));
			writer.close();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public String SaveFeatureDatasetA(String Path, String classy, String namefile) {
		double[] feature = null;
		try {
			GTCC gtcc = new GTCC();
			double[] data = StdAudio.read(Path);
			feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
		} catch (Exception e) {
			// TODO: handle exception
			return "Audio can't be processed, errorcode = " + e.toString();
		}

		try {
			Database db = new Database();
			db.insertGtccDataset(classy, feature,namefile, databasea);
			
		} catch (Exception e) {
			// TODO: handle exception
			return "Feature insert failed, errorcode= " + e.toString();
		}

		return "Feature from "+namefile+" inserted succesfully";
	}

	public String SaveFeatureDatasetAYin(String Path, String classy, String namefile) {
		double[] feature = null;
		try {
			Yin yin = new Yin(16000);
			yin.main(Path);
			ArrayList<Float> yins = yin.getPitchs();
			feature = new double[yins.size()];
			for (int i = 0; i < yins.size(); i++) {
				feature[i] = Double.valueOf(yins.get(i).toString());
			}

		} catch (Exception e) {
			// TODO: handle exception
			return "Audio can't be processed, errorcode = " + e.toString();
		}

		try {
			Database db = new Database();
			db.insertGtccDataset(classy, feature, namefile, databaseayin);
		} catch (Exception e) {
			// TODO: handle exception
			return "Feature insert failed, errorcode= " + e.toString();
		}

		return "Feature from "+namefile+" inserted succesfully";
	}

	public String SaveFeatureDatasetBYin(String Path, String classy, String namefile) {
		double[] feature = null;
		try {
			Yin yin = new Yin(16000);
			yin.main(Path);
			ArrayList<Float> yins = yin.getPitchs();
			feature = new double[yins.size()];
			for (int i = 0; i < yins.size(); i++) {
				feature[i] = Double.valueOf(yins.get(i).toString());
			}

		} catch (Exception e) {
			// TODO: handle exception
			return "Audio can't be processed, errorcode = " + e.toString();
		}

		try {
			Database db = new Database();
			db.insertGtccDataset(classy, feature, namefile, databasebyin);
		} catch (Exception e) {
			// TODO: handle exception
			return "Feature insert failed, errorcode= " + e.toString();
		}

		return "Feature from "+namefile+" inserted succesfully";
	}

	public String SaveFeatureDatasetB(String Path, String classy, String namefile) {
		double[] feature = null;
		try {
			GTCC gtcc = new GTCC();
			double[] data = StdAudio.read(Path);
			feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
		} catch (Exception e) {
			// TODO: handle exception
			return "Audio can't be processed, errorcode = " + e.toString();
		}

		try {
			Database db = new Database();
			db.insertGtccDataset(classy, feature, namefile, databaseb);
		} catch (Exception e) {
			// TODO: handle exception
			return "Feature insert failed, errorcode= " + e.toString();
		}

		return "Feature from "+namefile+" inserted succesfully";
	}
	
	public String TrainAsetWorkNoArtifact() {
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[3][];
		double[][][] covariances = new double[3][][];
		double[] phi = new double[3];
		int[] sizedata = new int[3];
		Database db = new Database();
		int[][] index = new int[3][];

		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databasea);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databasea);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */

		}

		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}

		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.;
				}

			}
		}

		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();

		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}

		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);

		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			//System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}
		//GaussianMixtureModel gmm = new GaussianMixtureModel(featurearray, 4, means);
		 GaussianMixtureModel gmm = new GaussianMixtureModel(means,covariances,phi,3);
		// gmm.EMApache(featurearray);
		gmm.EMKmeansWork(featurearray);
		double accuracy = predictaTrainNoArtifact();

		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmma);

		return String.valueOf(gmm.getLoglikelihoods());

	}
	
	public Map<String, String> TrainAsetWeb() {
		/*ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[4][];
		double[][][] covariances = new double[4][][];
		double[] phi = new double[4];
		int[] sizedata = new int[4];
		Database db = new Database();
		int[][] index = new int[4][];
		Main.status="Collecting data from DB";
		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databasea);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databasea);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 

		}

		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}

		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.;
				}

			}
		}

		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();

		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}

		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);

		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			//System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}
		Main.status="Starting training process ... ";
		//GaussianMixtureModel gmm = new GaussianMixtureModel(featurearray, 4, means);
		GaussianMixtureModel gmm;
		try {
			gmm = new GaussianMixtureModel(means,covariances,phi,4);
			gmm.EMKmeansWork(featurearray);
		} catch (Exception e) {
			// TODO: handle exception
			gmm = new GaussianMixtureModel(featurearray, 4, means);
			gmm.EMKmeansWork(featurearray);
		}*/
		
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[4][];
		double[][][] covariances = new double[4][][];
		double[] phi = new double[4];
		int[] sizedata = new int[4];
		Database db = new Database();
		int[][] index = new int[4][];
	
		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databasea);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databasea);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);
	
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();
	
			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */
	
		}
	
		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}
	
		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.;
				}
	
			}
		}
	
		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();
	
		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}
	
		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);
	
		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();
	
		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {
	
				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;
	
		}
	
		featurearray = newdata;
	
		//System.out.println(Arrays.deepToString(featurearray));
	
		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */
	
		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);
	
						// System.out.println(k2);
						sums += featurearray[j2][i];
					}
	
				}
				meanstemp.add(sums / sizedata[j]);
			}
	
			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}
	
			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}
	
		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));
	
			covariances[j] = cova.getArray();
			//System.out.println(Arrays.deepToString(covariances[j]));
		}
	
		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);
	
		}
		GaussianMixtureModel gmm;
		try {
			gmm = new GaussianMixtureModel(means,covariances,phi,4);
			gmm.EMKmeansWork(featurearray);
		} catch (Exception e) {
			// TODO: handle exception
			gmm = new GaussianMixtureModel(featurearray, 4, means);
			gmm.EMKmeansWork(featurearray);
		}
		
		double accuracy = predictaTrainWeb();

		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmma);
		Map<String, String> model = new HashMap<>();
		for (int i = 0; i < sizedata.length; i++) {
			model.put("class"+i, String.valueOf(sizedata[i]));
		}
		model.put("accuracy", String.valueOf(accuracy));
		model.put("likelihood", String.valueOf(gmm.getLoglikelihoods()));
		return model;

	}
	

	public String TrainAsetWork() {
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[4][];
		double[][][] covariances = new double[4][][];
		double[] phi = new double[4];
		int[] sizedata = new int[4];
		Database db = new Database();
		int[][] index = new int[4][];
	
		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databasea);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databasea);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);
	
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();
	
			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */
	
		}
	
		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}
	
		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.;
				}
	
			}
		}
	
		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();
	
		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}
	
		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);
	
		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();
	
		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {
	
				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;
	
		}
	
		featurearray = newdata;
	
		//System.out.println(Arrays.deepToString(featurearray));
	
		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */
	
		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);
	
						// System.out.println(k2);
						sums += featurearray[j2][i];
					}
	
				}
				meanstemp.add(sums / sizedata[j]);
			}
	
			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}
	
			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}
	
		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));
	
			covariances[j] = cova.getArray();
			//System.out.println(Arrays.deepToString(covariances[j]));
		}
	
		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);
	
		}
		GaussianMixtureModel gmm;
		try {
			gmm = new GaussianMixtureModel(means,covariances,phi,4);
			gmm.EMKmeansWork(featurearray);
		} catch (Exception e) {
			// TODO: handle exception
			gmm = new GaussianMixtureModel(featurearray, 4, means);
			gmm.EMKmeansWork(featurearray);
		}
		//GaussianMixtureModel gmm = new GaussianMixtureModel(featurearray, 4, means);
		// gmm = new GaussianMixtureModel(means,covariances,phi,4);
		// gmm.EMApache(featurearray);
		//gmm.EMKmeansWork(featurearray);
		double accuracy = predictaTrain();
		
		toJson(covariances, means, featurearray);
		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmma);
	
		return String.valueOf(gmm.getLoglikelihoods());
	
	}

	public String TrainAsetWorkYINNoArtifact() {
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[3][];
		double[][][] covariances = new double[3][][];
		double[] phi = new double[3];
		int[] sizedata = new int[3];
		Database db = new Database();
		int[][] index = new int[3][];

		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databaseayin);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databaseayin);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */

		}

		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}

		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.;
				}

			}
		}

		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();

		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}

		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);

		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}
		GaussianMixtureModel gmm = new GaussianMixtureModel(featurearray, 3, means);
		// GaussianMixtureModel gmm = new GaussianMixtureModel(means,covariances,phi,4);
		// gmm.EMApache(featurearray);
		gmm.EMKmeansWork(featurearray);
		// gmm.fit(featurearray,4);
		double accuracy = predictayinTrainNoArtifact();

		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmmyina);
		return String.valueOf(gmm.getLoglikelihoods());

	}
	
	public Map<String, String> TrainAsetWebYIN() {
		/*ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[4][];
		double[][][] covariances = new double[4][][];
		double[] phi = new double[4];
		int[] sizedata = new int[4];
		Database db = new Database();
		int[][] index = new int[4][];

		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databaseayin);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databaseayin);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 

		}

		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}

		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.;
				}

			}
		}

		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();

		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}

		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);

		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}
		GaussianMixtureModel gmm;
		try {
			System.out.println("try");
			 gmm= new GaussianMixtureModel(featurearray, 4, means);
			
			 gmm.EMKmeansWork(featurearray);
		} catch (Exception e) {
			// TODO: handle exception
			System.out.println("catch");
			gmm = new GaussianMixtureModel(means,covariances,phi,4);
			 gmm.EMKmeansWork(featurearray);
		}*/
		
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[4][];
		double[][][] covariances = new double[4][][];
		double[] phi = new double[4];
		int[] sizedata = new int[4];
		Database db = new Database();
		int[][] index = new int[4][];

		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databaseayin);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databaseayin);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */

		}

		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}

		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.+Math.random();
				}

			}
		}

		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();

		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}

		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);

		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}
		GaussianMixtureModel gmm = null ;
		try {
			gmm = new GaussianMixtureModel(means,covariances,phi,4);
			
			gmm.EMKmeansWork(featurearray);
		} catch (Exception e) {
			// TODO: handle exception
			try {
				System.out.println("catch");
				gmm = new GaussianMixtureModel(featurearray, 4, means);
				gmm.EMKmeansWork(featurearray);
			} catch (Exception e2) {
				// TODO: handle exception
				TrainAsetWebYIN();
			}
			
		}
		
		// gmm.fit(featurearray,4);
		double accuracy = predictayinTrainWeb();

		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmmyina);
		Map<String, String> model = new HashMap<>();
		for (int i = 0; i < sizedata.length; i++) {
			model.put("class"+i, String.valueOf(sizedata[i]));
		}
		model.put("accuracy", String.valueOf(accuracy));
		model.put("likelihood", String.valueOf(gmm.getLoglikelihoods()));
		return model;

	}
	
	
	public String TrainAsetWorkYIN() {
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[4][];
		double[][][] covariances = new double[4][][];
		double[] phi = new double[4];
		int[] sizedata = new int[4];
		Database db = new Database();
		int[][] index = new int[4][];

		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databaseayin);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databaseayin);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */

		}

		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}

		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.+Math.random();
				}

			}
		}

		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();

		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}

		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);

		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}
		GaussianMixtureModel gmm ;
		try {
			gmm = new GaussianMixtureModel(means,covariances,phi,4);
			
			gmm.EMKmeansWork(featurearray);
		} catch (Exception e) {
			// TODO: handle exception
			System.out.println("catch");
			gmm = new GaussianMixtureModel(featurearray, 4, means);
			gmm.EMKmeansWork(featurearray);
		}
		
		// GaussianMixtureModel gmm = new GaussianMixtureModel(means,covariances,phi,4);
		// gmm.EMApache(featurearray);
		
		// gmm.fit(featurearray,4);
		double accuracy = predictayinTrain();
		toJson(covariances, means, featurearray);
		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmmyina);
		return String.valueOf(gmm.getLoglikelihoods());

	}


	public Map<String, String> TrainBsetWebYIN() {
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[3][];
		double[][][] covariances = new double[3][][];
		double[] phi = new double[3];
		int[] sizedata = new int[3];
		Database db = new Database();
		int[][] index = new int[3][];

		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databasebyin);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databasebyin);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */

		}

		// find max
		System.out.println("find max ");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}
max+=5;
		System.out.println("fitting with max  "+max);
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.+Math.random();
				}

			}
		}
if (max>featurearray.length) {
	System.out.println("Reduce dimension");
	RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
	SingularValueDecomposition svd = new SingularValueDecomposition(rm);
	RealMatrix u = svd.getU();
	RealMatrix s = svd.getS();
	RealMatrix vt = svd.getVT();
	double[][] sar = s.getData();
	double[][] vtar = vt.getData();

	double[][] news = new double[featurearray.length][featurearray.length];
	double[][] newv = new double[featurearray.length][featurearray.length];
	for (int j = 0; j < featurearray.length; j++) {
		for (int j2 = 0; j2 < featurearray.length; j2++) {
			news[j][j2] = sar[j][j2];
			newv[j][j2] = vtar[j][j2];
		}
	}

	s = MatrixUtils.createRealMatrix(news);
	vt = MatrixUtils.createRealMatrix(newv);

	RealMatrix featuresvd = u.multiply(s).multiply(vt);
	int featlenght = featurearray.length;
	featurearray = new double[featlenght][featlenght];
	featurearray = featuresvd.getData();
}
else {
	System.out.println("No need Reduction");
}
		

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}GaussianMixtureModel gmm;
		try {
			System.out.println("try");
			 gmm= new GaussianMixtureModel(featurearray, 3, means);
			
			 gmm.EMKmeansWork(featurearray);
		} catch (Exception e) {
			// TODO: handle exception
			System.out.println("catch");
			gmm = new GaussianMixtureModel(means,covariances,phi,3);
			 gmm.EMKmeansWork(featurearray);
		}
		//
		
		// gmm.EMApache(featurearray);
		
		// gmm.fit(featurearray,4);
		double accuracy = predictbyinTrainWeb();

		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmmyinb);

		Map<String, String> model = new HashMap<>();
		for (int i = 0; i < sizedata.length; i++) {
			model.put("class"+i, String.valueOf(sizedata[i]));
		}
		model.put("accuracy", String.valueOf(accuracy));
		model.put("likelihood", String.valueOf(gmm.getLoglikelihoods()));
		return model;

	}

	public String TrainBsetWorkYIN() {
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[3][];
		double[][][] covariances = new double[3][][];
		double[] phi = new double[3];
		int[] sizedata = new int[3];
		Database db = new Database();
		int[][] index = new int[3][];

		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databasebyin);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databasebyin);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */

		}

		// find max
		System.out.println("find max ");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}
		System.out.println("fitting with max  "+max);
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.+Math.random();
				}

			}
		}
if (max>featurearray.length) {
	System.out.println("Reduce dimension");
	RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
	SingularValueDecomposition svd = new SingularValueDecomposition(rm);
	RealMatrix u = svd.getU();
	RealMatrix s = svd.getS();
	RealMatrix vt = svd.getVT();
	double[][] sar = s.getData();
	double[][] vtar = vt.getData();

	double[][] news = new double[featurearray.length][featurearray.length];
	double[][] newv = new double[featurearray.length][featurearray.length];
	for (int j = 0; j < featurearray.length; j++) {
		for (int j2 = 0; j2 < featurearray.length; j2++) {
			news[j][j2] = sar[j][j2];
			newv[j][j2] = vtar[j][j2];
		}
	}

	s = MatrixUtils.createRealMatrix(news);
	vt = MatrixUtils.createRealMatrix(newv);

	RealMatrix featuresvd = u.multiply(s).multiply(vt);
	int featlenght = featurearray.length;
	featurearray = new double[featlenght][featlenght];
	featurearray = featuresvd.getData();
}
else {
	System.out.println("No need Reduction");
}
		

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}GaussianMixtureModel gmm;
		try {
			System.out.println("try");
			 gmm= new GaussianMixtureModel(featurearray, 3, means);
			
			 gmm.EMKmeansWork(featurearray);
		} catch (Exception e) {
			// TODO: handle exception
			System.out.println("catch");
			gmm = new GaussianMixtureModel(means,covariances,phi,3);
			 gmm.EMKmeansWork(featurearray);
		}
		//
		
		// gmm.EMApache(featurearray);
		
		// gmm.fit(featurearray,4);
		double accuracy = predictbyinTrain();
		toJson(covariances, means, featurearray);
		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmmyinb);

		return String.valueOf(gmm.getLoglikelihoods());

	}
	public Map<String, String> TrainBsetWeb() {
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[3][];
		double[][][] covariances = new double[3][][];
		double[] phi = new double[3];
		int[] sizedata = new int[3];
		Database db = new Database();
		int[][] index = new int[3][];

		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databaseb);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databaseb);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */

		}

		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}

		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.;
				}

			}
		}

		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();

		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}

		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);

		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}
		GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances, phi, 3);
		// GaussianMixtureModel gmm = new GaussianMixtureModel(featurearray, 3,means);
		// GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances,
		// phi,4);
		// gmm.EMKmeans1(featurearray,4);
		 //gmm.EMApache(featurearray);
		gmm.EMKmeansWork(featurearray);
		// gmm.fit(featurearray,4);
		double accuracy = predictbTrainWeb();

		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmmb);

		Map<String, String> model = new HashMap<>();
		for (int i = 0; i < sizedata.length; i++) {
			model.put("class"+i, String.valueOf(sizedata[i]));
		}
		model.put("accuracy", String.valueOf(accuracy));
		model.put("likelihood", String.valueOf(gmm.getLoglikelihoods()));
		return model;

	}
	
	
	public String TrainBsetWork() {
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		ArrayList<Integer> label = new ArrayList<>();
		double[][] means = new double[3][];
		double[][][] covariances = new double[3][][];
		double[] phi = new double[3];
		int[] sizedata = new int[3];
		Database db = new Database();
		int[][] index = new int[3][];

		System.out.println("Collecting data from DB");
		for (int i = 0; i < index.length; i++) {
			try {
				index[i] = db.selectGtccDatasetIndex(i, databaseb);
				for (int j = 0; j < index[i].length; j++) {
					ArrayList<Double> featuretemp = new ArrayList<>();
					ObjectFeatureGtcc of = new ObjectFeatureGtcc();
					db.selectGtccDataset(index[i][j], databaseb);
					try {
						FileInputStream fileIn = new FileInputStream("Out.ser");
						ObjectInputStream in = new ObjectInputStream(fileIn);
						of = new ObjectFeatureGtcc();
						of = (ObjectFeatureGtcc) in.readObject();
						in.close();
						fileIn.close();
					} catch (Exception e) {
						// TODO: handle exception
					}
					for (int k = 0; k < of.data.length; k++) {
						featuretemp.add(of.data[k]);
					}
					feature.add(featuretemp);
					label.add(i);

				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// sizedata[i] = feature.get(i).size();

			/*
			 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
			 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
			 * feature.get(i1).size(); } }
			 */

		}

		// find max
		System.out.println("find max");
		int max = feature.get(0).size();
		for (int j = 0; j < feature.size(); j++) {
			if (feature.get(j).size() > max) {
				max = feature.get(j).size();
			}
		}

		System.out.println("fitting with max");
		// Convert to 2d array & fitting with max
		double[][] featurearray = new double[feature.size()][max];
		for (int j1 = 0; j1 < feature.size(); j1++) {
			for (int j2 = 0; j2 < max; j2++) {
				if (j2 < feature.get(j1).size()) {
					featurearray[j1][j2] = feature.get(j1).get(j2);
				} else {
					featurearray[j1][j2] = 0.;
				}

			}
		}

		System.out.println("Reduce dimension");
		RealMatrix rm = MatrixUtils.createRealMatrix(featurearray);
		SingularValueDecomposition svd = new SingularValueDecomposition(rm);
		RealMatrix u = svd.getU();
		RealMatrix s = svd.getS();
		RealMatrix vt = svd.getVT();
		double[][] sar = s.getData();
		double[][] vtar = vt.getData();

		double[][] news = new double[featurearray.length][featurearray.length];
		double[][] newv = new double[featurearray.length][featurearray.length];
		for (int j = 0; j < featurearray.length; j++) {
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				news[j][j2] = sar[j][j2];
				newv[j][j2] = vtar[j][j2];
			}
		}

		s = MatrixUtils.createRealMatrix(news);
		vt = MatrixUtils.createRealMatrix(newv);

		RealMatrix featuresvd = u.multiply(s).multiply(vt);
		int featlenght = featurearray.length;
		featurearray = new double[featlenght][featlenght];
		featurearray = featuresvd.getData();

		ArrayList<Double> datas = new ArrayList<>();
		for (int i = 0; i < featurearray.length; i++) {
			for (int j = 0; j < featurearray[i].length; j++) {
				datas.add(featurearray[i][j]);
			}
		}
		double[] arraytempdata = new double[datas.size()];
		for (int i = 0; i < datas.size(); i++) {
			arraytempdata[i] = datas.get(i);
		}
		DescriptiveStatistics ds = new DescriptiveStatistics(arraytempdata);
		double mean = ds.getMean();
		double stdev = ds.getStandardDeviation();
		double[][] newdata = new double[featurearray.length][featurearray[0].length];
		for (int i = 0; i < featurearray.length; i++) {
			double[] temp = new double[featurearray[i].length];
			for (int j = 0; j < featurearray[0].length; j++) {

				temp[j] = (featurearray[i][j] - mean) / stdev;
			}
			newdata[i] = temp;

		}

		featurearray = newdata;

		System.out.println(Arrays.deepToString(featurearray));

		// find occurence
		for (int i = 0; i < index.length; i++) {
			int count = 0;
			for (int j = 0; j < featurearray.length; j++) {
				if (label.get(j) == i) {
					count++;
				}
			}
			sizedata[i] = count;
		}
		/*
		 * System.out.println("Compute Means"); // compute means ArrayList<Double>
		 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
		 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
		 * (label.get(j2) == j) { //System.out.println(j2 + " " +
		 * featurearray[j2].length + " " + sizedata[j]); meanstemp = new ArrayList<>();
		 * double sum = 0.; for (int k = 0; k < featurearray[j2].length; k++) { for (int
		 * k2 = 0; k2 < sizedata[j]; k2++) { //System.out.println(k2); sum +=
		 * featurearray[k2][k]; } meanstemp.add(sum/sizedata[j]); }
		 * 
		 * }
		 * 
		 * }
		 */

		System.out.println("Compute Means");
		// compute means
		ArrayList<Double> meanstemp = new ArrayList<>();
		double[] meansarray = null;
		double sums = 0;
		for (int j = 0; j < index.length; j++) {
			meanstemp = new ArrayList<>();
			for (int i = 0; i < featurearray[0].length; i++) {
				sums = 0.;
				for (int j2 = 0; j2 < featurearray.length; j2++) {
					if (label.get(j2) == j) {
						// System.out.println(j2 + " " + featurearray[j2].length + " " + sizedata[j]);

						// System.out.println(k2);
						sums += featurearray[j2][i];
					}

				}
				meanstemp.add(sums / sizedata[j]);
			}

			meansarray = new double[meanstemp.size()];
			for (int j1 = 0; j1 < meanstemp.size(); j1++) {
				meansarray[j1] = meanstemp.get(j1);
			}

			means[j] = meansarray;
			System.out.println(Arrays.toString(means[j]));
		}

		System.out.println("Compute covariances");
		// compute covariances
		for (int j = 0; j < index.length; j++) {
			double[][] datatemp = new double[sizedata[j]][];
			int iter = 0;
			for (int j2 = 0; j2 < featurearray.length; j2++) {
				if (label.get(j2) == j) {
					datatemp[iter] = featurearray[j2];
					iter++;
				}
			}
			/*
			 * RealMatrix cov = MatrixUtils.createRealMatrix(datatemp); cov = new
			 * Covariance(cov).getCovarianceMatrix(); Matrix cova = new
			 * Matrix(cov.getData()); cova = cova.plus(Matrix.identity(featlenght,
			 * featlenght).times(0.1)); covariances[j] = cova.getArray();
			 * System.out.println(Arrays.deepToString(covariances[j])); EigenDecomposition
			 * ed = new EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
			 * System.out.println(ed.getDeterminant());
			 */
			Matrix cova;
			cova = new Matrix(datatemp, sizedata[j], datatemp[0].length);
			cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]);
			cova = cova.plus(Matrix.identity(datatemp[0].length, datatemp[0].length).times(0.1));

			covariances[j] = cova.getArray();
			System.out.println(Arrays.deepToString(covariances[j]));
		}

		// compute phi
		double sum = 0;
		for (int j = 0; j < sizedata.length; j++) {
			sum += sizedata[j];
		}
		for (int j = 0; j < phi.length; j++) {
			phi[j] = sizedata[j] / sum;
			System.out.println(phi[j] + " = " + sizedata[j] + " / " + sum);

		}
		GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances, phi, 3);
		// GaussianMixtureModel gmm = new GaussianMixtureModel(featurearray, 3,means);
		// GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances,
		// phi,4);
		// gmm.EMKmeans1(featurearray,4);
		 //gmm.EMApache(featurearray);
		gmm.EMKmeansWork(featurearray);
		// gmm.fit(featurearray,4);
		double accuracy = predictbTrain();

		// gmm.fit(featurearray,4);
		db.insertGMMParam(gmm.covariances, gmm.means, gmm.pi, accuracy, databasegmmb);
		toJson(covariances, means, featurearray);
		return String.valueOf(gmm.getLoglikelihoods());

	}
	
	public static int PredictAsetTrain(String path) {
		// get latest parameter
		int result = 0;
		try {
			/*Database db = new Database();
			db.selectGtccParam(databasegmma);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();
*/
			double [][] means = GaussianMixtureModel.means;
			double[][][] covariances = GaussianMixtureModel.covariances;
			double[] pi = GaussianMixtureModel.pi;
			GTCC gtcc = new GTCC();
			double[] data = StdAudio.read(path);
			double[] feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
			
			/*GTCCMatlab gm = new GTCCMatlab();
			double [] feature =gm.extractFeature(path);
			gm.closeMatlab();
			*/

			GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances, pi, 4);
			result = gmm.predict(feature, means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}
	
	public static int PredictAsetTrainNoArtifact(String path) {
		// get latest parameter
		int result = 0;
		try {
			/*Database db = new Database();
			db.selectGtccParam(databasegmma);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();
*/
			double [][] means = GaussianMixtureModel.means;
			double[][][] covariances = GaussianMixtureModel.covariances;
			double[] pi = GaussianMixtureModel.pi;
			GTCC gtcc = new GTCC();
			double[] data = StdAudio.read(path);
			double[] feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
			
			/*GTCCMatlab gm = new GTCCMatlab();
			double [] feature =gm.extractFeature(path);
			gm.closeMatlab();
			*/

			GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances, pi, 3);
			result = gmm.predict(feature, means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}

	public static int PredictAset(String path) {
		// get latest parameter
		int result = 0;
		try {
			Database db = new Database();
			db.selectGtccParam(databasegmma);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();

			GTCC gtcc = new GTCC();
			double[] data = StdAudio.read(path);
			StdAudio.close();
			double[] feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);

			GaussianMixtureModel gmm = new GaussianMixtureModel(of.means, of.covariances, of.pi, 4);
			result = gmm.predict(feature, of.means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}
	
	public static int PredictAsetNoArtifact(String path) {
		// get latest parameter
		int result = 0;
		try {
			Database db = new Database();
			db.selectGtccParam(databasegmma);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();

			GTCC gtcc = new GTCC();
			double[] data = StdAudio.read(path);
			double[] feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);

			GaussianMixtureModel gmm = new GaussianMixtureModel(of.means, of.covariances, of.pi, 3);
			result = gmm.predict(feature, of.means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}
	
	public double predictaTrain() {

		int truelabel = 0, falselabel = 0;
		int count = 0;
		ArrayList<Integer> list = new ArrayList<>();
		ConfusionMatrix cm = new ConfusionMatrix();
		for (int i = 0; i < 4; i++) {
			int nol = 0,satu = 0,dua=0,tiga=0;
			//String path = "0Norm\\" + i;
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\a\\" + i;
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\a\\pisah\\"+i;
			//String path = "0Real\\" + i;
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\atest\\"+i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = PredictAsetTrain(child.getAbsolutePath());
					list.add(result);
					if (i == result) {
						truelabel++;

					} else {
						falselabel++;
					}
					if (result==0) {
						nol++;
					}
					else if (result==1) {
						satu++;
					}
					else if (result==2) {
						dua++;
					}
					else {
						tiga++;
					}
					// System.out.println(" file " + child.getName() + " prediksi = " + result);

				}
			}
			
			if (i==0) {
				System.out.println("nol");
				cm.increaseValue("Normal", "Normal",nol);
				cm.increaseValue("Normal", "Murmur",satu);
				cm.increaseValue("Normal", "Extrasystole",dua);
				cm.increaseValue("Normal", "Artifact",tiga);
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
			if (i==1) {
				System.out.println("satu");
				cm.increaseValue("Murmur", "Normal",nol);
				cm.increaseValue("Murmur", "Murmur",satu);
				cm.increaseValue("Murmur", "Extrasystole",dua);
				cm.increaseValue("Murmur", "Artifact",tiga);
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
			if (i==2) {
				System.out.println("dua");
				cm.increaseValue("Extrasystole", "Normal",nol);
				cm.increaseValue("Extrasystole", "Murmur",satu);
				cm.increaseValue("Extrasystole", "Extrasystole",dua);
				cm.increaseValue("Extrasystole", "Artifact",tiga);
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
			if (i==3) {
				System.out.println("tiga");
				cm.increaseValue("Artifact", "Normal",nol);
				cm.increaseValue("Artifact", "Murmur",satu);
				cm.increaseValue("Artifact", "Extrasystole",dua);
				cm.increaseValue("Artifact", "Artifact",tiga);
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
		}
		
		System.out.println(cm);
		System.out.println(cm.printLabelPrecRecFm());
		System.out.println(cm.getPrecisionForLabels());
		System.out.println(cm.getRecallForLabels());
		System.out.println(cm.printNiceResults());
		
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
		toJsonPredictTrain(list);
		return accuracy;

	}
	
	public double predictaTrainWeb() {

		int truelabel = 0, falselabel = 0;
		int count = 0;
		for (int i = 0; i < 4; i++) {
			//String path = "0Real\\" + i;
			String path = "upload\\0\\"+i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = PredictAsetTrain(child.getAbsolutePath());
					if (i == result) {
						truelabel++;

					} else {
						falselabel++;
					}
					System.out.println(" file " + child.getName() + " real ="+i+" prediksi = " + result);
				}
			}
		}
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);

		return accuracy;

	}
	
	public double predictaTrainNoArtifact() {

		int truelabel = 0, falselabel = 0;
		int count = 0;
		for (int i = 0; i < 3; i++) {
			String path = "0Real\\" + i;
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\atest\\"+i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = PredictAsetTrainNoArtifact(child.getAbsolutePath());
					if (i == result) {
						truelabel++;

					} else {
						falselabel++;
					}
					// System.out.println(" file " + child.getName() + " prediksi = " + result);

				}
			}
		}
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);

		return accuracy;

	}

	public double predicta() {

		int truelabel = 0, falselabel = 0;
		int count = 0;
		for (int i = 0; i < 4; i++) {
			String path = "0\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = PredictAset(child.getAbsolutePath());
					if (i == result) {
						truelabel++;

					} else {
						falselabel++;
					}
					// System.out.println(" file " + child.getName() + " prediksi = " + result);

				}
			}
		}
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);

		return accuracy;

	}
	
	public double predictbTrain() {

		int truelabel = 0, falselabel = 0;
		int count = 0;
		ArrayList<Integer> list = new ArrayList<>();
		ConfusionMatrix cm = new ConfusionMatrix();
		for (int i = 0; i < 3; i++) {
			int nol = 0,satu = 0,dua=0,tiga=0;
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\btest\\" + i;
			//String path = "1\\" + i;
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\b\\" + i;
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\b\\pisah\\"+i;
			//String path = "150\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = PredictBsetTrain(child.getAbsolutePath());
					list.add(result);
					if (i == result) {
						truelabel++;

					} else {
						falselabel++;
					}if (result==0) {
						nol++;
					}
					else if (result==1) {
						satu++;
					}
					else if (result==2) {
						dua++;
					}
					
					// System.out.println(" file " + child.getName() + " prediksi = " + result);

				}
			}
			
			if (i==0) {
				System.out.println("nol");
				cm.increaseValue("Normal", "Normal",nol);
				cm.increaseValue("Normal", "Murmur",satu);
				cm.increaseValue("Normal", "Extrasystole",dua);
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
			if (i==1) {
				System.out.println("satu");
				cm.increaseValue("Murmur", "Normal",nol);
				cm.increaseValue("Murmur", "Murmur",satu);
				cm.increaseValue("Murmur", "Extrasystole",dua);
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
			if (i==2) {
				System.out.println("dua");
				cm.increaseValue("Extrasystole", "Normal",nol);
				cm.increaseValue("Extrasystole", "Murmur",satu);
				cm.increaseValue("Extrasystole", "Extrasystole",dua);
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
			
		}
		
		System.out.println(cm);
		System.out.println(cm.printLabelPrecRecFm());
		System.out.println(cm.getPrecisionForLabels());
		System.out.println(cm.getRecallForLabels());
		System.out.println(cm.printNiceResults());
		toJsonPredictTrain(list);
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);

		return accuracy;

	}
	
	public double predictbTrainWeb() {

		int truelabel = 0, falselabel = 0;
		int count = 0;
		for (int i = 0; i < 3; i++) {
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\btest\\" + i;
			String path = "upload\\1\\" + i;
			//String path = "150\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = PredictBsetTrain(child.getAbsolutePath());
					if (i == result) {
						truelabel++;

					} else {
						falselabel++;
					}
					System.out.println(" file " + child.getName() + " real ="+i+" prediksi = " + result);

				}
			}
		}
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);

		return accuracy;

	}
	
	
public static double predictb() {
		
		int truelabel = 0, falselabel = 0;
		int count = 0;
		for (int i = 0; i < 3; i++) {
			String path = "1Real//" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = PredictBset(child.getAbsolutePath());
					if (i == result) {
						truelabel++;

					} else {
						falselabel++;
					}
					//System.out.println(" file " + child.getName() + " prediksi = " + result);

				}
			}
		}
		double accuracy = ((double )truelabel / (double)count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
		return accuracy;
	}


public static double predictayin() {
	
	int truelabel = 0, falselabel = 0;
	int count = 0;
	for (int i = 0; i < 4; i++) {
		String path = "0//" + i;
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			for (File child : directoryListing) {
				count++;
				int result = PredictAsetYin(child.getAbsolutePath());
				if (i == result) {
					truelabel++;

				} else {
					falselabel++;
				}
				//System.out.println(" file " + child.getName() + " prediksi = " + result);

			}
		}
	}
	double accuracy = ((double )truelabel / (double)count) * 100;
	System.out.println(
			"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
	return accuracy;
}

public static double predictayinTrain() {
	
	int truelabel = 0, falselabel = 0;
	int count = 0;
	ArrayList<Integer> list = new ArrayList<>();
	ConfusionMatrix cm = new ConfusionMatrix();
	for (int i = 0; i < 4; i++) {
		int nol = 0,satu = 0,dua=0,tiga=0;
		//String path = "0Real//" + i;
		//String path = "0Norm//" + i;
		String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\a\\pisah\\"+i;
		//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\a\\"+i;
		
		//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\atest\\"+i;
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			for (File child : directoryListing) {
				count++;
				int result = PredictAsetYinTrain(child.getAbsolutePath());
				list.add(i);
				if (i == result) {
					truelabel++;

				} else {
					falselabel++;
				}
				if (result==0) {
					nol++;
				}
				else if (result==1) {
					satu++;
				}
				else if (result==2) {
					dua++;
				}
				else {
					tiga++;
				}
				//System.out.println(" file " + child.getName() + " prediksi = " + result);

			}
		}if (i==0) {
			System.out.println("nol");
			cm.increaseValue("Normal", "Normal",nol);
			cm.increaseValue("Normal", "Murmur",satu);
			cm.increaseValue("Normal", "Extrasystole",dua);
			cm.increaseValue("Normal", "Artifact",tiga);
			nol=0;
			satu=0;
			dua=0;
			tiga=0;
		}
		if (i==1) {
			System.out.println("satu");
			cm.increaseValue("Murmur", "Normal",nol);
			cm.increaseValue("Murmur", "Murmur",satu);
			cm.increaseValue("Murmur", "Extrasystole",dua);
			cm.increaseValue("Murmur", "Artifact",tiga);
			nol=0;
			satu=0;
			dua=0;
			tiga=0;
		}
		if (i==2) {
			System.out.println("dua");
			cm.increaseValue("Extrasystole", "Normal",nol);
			cm.increaseValue("Extrasystole", "Murmur",satu);
			cm.increaseValue("Extrasystole", "Extrasystole",dua);
			cm.increaseValue("Extrasystole", "Artifact",tiga);
			nol=0;
			satu=0;
			dua=0;
			tiga=0;
		}
		if (i==3) {
			System.out.println("tiga");
			cm.increaseValue("Artifact", "Normal",nol);
			cm.increaseValue("Artifact", "Murmur",satu);
			cm.increaseValue("Artifact", "Extrasystole",dua);
			cm.increaseValue("Artifact", "Artifact",tiga);
			nol=0;
			satu=0;
			dua=0;
			tiga=0;
		}
	}
	
	System.out.println(cm);
	System.out.println(cm.printLabelPrecRecFm());
	System.out.println(cm.getPrecisionForLabels());
	System.out.println(cm.getRecallForLabels());
	System.out.println(cm.printNiceResults());
	
	double accuracy = ((double )truelabel / (double)count) * 100;
	System.out.println(
			"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
	toJsonPredictTrain(list);
	return accuracy;
}

public static double predictayinTrainWeb() {
	
	int truelabel = 0, falselabel = 0;
	int count = 0;
	for (int i = 0; i < 4; i++) {
		//String path = "0//" + i;
		String path = "upload\\0\\"+i;
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			for (File child : directoryListing) {
				count++;
				int result = PredictAsetYinTrain(child.getAbsolutePath());
				if (i == result) {
					truelabel++;

				} else {
					falselabel++;
				}
				System.out.println(" file " + child.getName() + " real ="+i+" prediksi = " + result);
			}
		}
	}
	double accuracy = ((double )truelabel / (double)count) * 100;
	System.out.println(
			"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
	return accuracy;
}


public static double predictayinTrainNoArtifact() {
	
	int truelabel = 0, falselabel = 0;
	int count = 0;
	for (int i = 0; i < 3; i++) {
		String path = "0//" + i;
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			for (File child : directoryListing) {
				count++;
				int result = PredictAsetYinTrainNoArtifact(child.getAbsolutePath());
				if (i == result) {
					truelabel++;

				} else {
					falselabel++;
				}
				//System.out.println(" file " + child.getName() + " prediksi = " + result);

			}
		}
	}
	double accuracy = ((double )truelabel / (double)count) * 100;
	System.out.println(
			"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
	return accuracy;
}



public static double predictbyin() {
	
	int truelabel = 0, falselabel = 0;
	int count = 0;
	for (int i = 0; i < 3; i++) {
		String path = "1Real//" + i;
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			for (File child : directoryListing) {
				count++;
				int result = PredictBsetYin(child.getAbsolutePath());
				if (i == result) {
					truelabel++;

				} else {
					falselabel++;
				}
				//System.out.println(" file " + child.getName() + " prediksi = " + result);

			}
		}
	}
	double accuracy = ((double )truelabel / (double)count) * 100;
	System.out.println(
			"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
	return accuracy;
}

public static double predictbyinTrain() {
	
	int truelabel = 0, falselabel = 0;
	int count = 0;
	ArrayList<Integer> list = new ArrayList<>();
	ConfusionMatrix cm = new ConfusionMatrix();
	for (int i = 0; i < 3; i++) {
		int nol = 0,satu = 0,dua=0,tiga=0;
		//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\btest\\" + i;
		//String path = "150//" + i;
		//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\b\\"+i;
		String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\b\\pisah//" + i;
		
		//String path = "1//" + i;
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			for (File child : directoryListing) {
				count++;
				int result = PredictBsetYinTrain(child.getAbsolutePath());
				list.add(i);
				if (i == result) {
					truelabel++;

				} else {
					falselabel++;
				}
				if (result==0) {
					nol++;
				}
				else if (result==1) {
					satu++;
				}
				else if (result==2) {
					dua++;
				}
				
				//System.out.println(" file " + child.getName() + " prediksi = " + result);

			}
		}if (i==0) {
			System.out.println("nol");
			cm.increaseValue("Normal", "Normal",nol);
			cm.increaseValue("Normal", "Murmur",satu);
			cm.increaseValue("Normal", "Extrasystole",dua);
			nol=0;
			satu=0;
			dua=0;
			tiga=0;
		}
		if (i==1) {
			System.out.println("satu");
			cm.increaseValue("Murmur", "Normal",nol);
			cm.increaseValue("Murmur", "Murmur",satu);
			cm.increaseValue("Murmur", "Extrasystole",dua);
			nol=0;
			satu=0;
			dua=0;
			tiga=0;
		}
		if (i==2) {
			System.out.println("dua");
			cm.increaseValue("Extrasystole", "Normal",nol);
			cm.increaseValue("Extrasystole", "Murmur",satu);
			cm.increaseValue("Extrasystole", "Extrasystole",dua);
			nol=0;
			satu=0;
			dua=0;
			tiga=0;
		}
		
	}
	
	System.out.println(cm);
	System.out.println(cm.printLabelPrecRecFm());
	System.out.println(cm.getPrecisionForLabels());
	System.out.println(cm.getRecallForLabels());
	System.out.println(cm.printNiceResults());
	double accuracy = ((double )truelabel / (double)count) * 100;
	System.out.println(
			"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
	toJsonPredictTrain(list);
	return accuracy;
}

public static double predictbyinTrainWeb() {
	
	int truelabel = 0, falselabel = 0;
	int count = 0;
	for (int i = 0; i < 3; i++) {
		//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\btest\\" + i;
		String path = "upload//1//" + i;
		//String path = "1//" + i;
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			for (File child : directoryListing) {
				count++;
				int result = PredictBsetYinTrain(child.getAbsolutePath());
				if (i == result) {
					truelabel++;

				} else {
					falselabel++;
				}
				System.out.println(" file " + child.getName() + " real ="+i+" prediksi = " + result);
			}
		}
	}
	double accuracy = ((double )truelabel / (double)count) * 100;
	System.out.println(
			"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
	return accuracy;
}

public static int PredictBsetTrain(String path) {
	// get latest parameter
	int result = 0;
	try {
		/*Database db = new Database();
		db.selectGtccParam(databasegmma);
		ObjectParamGMM of;
		FileInputStream fileIn = new FileInputStream("Out.ser");
		ObjectInputStream in = new ObjectInputStream(fileIn);
		of = new ObjectParamGMM();
		of = (ObjectParamGMM) in.readObject();
		in.close();
		fileIn.close();
*/
		double [][] means = GaussianMixtureModel.means;
		double[][][] covariances = GaussianMixtureModel.covariances;
		double[] pi = GaussianMixtureModel.pi;
		GTCC gtcc = new GTCC();
		double[] data = StdAudio.read(path);
		double[] feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);

		GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances, pi, 3);
		result = gmm.predict(feature, means[0].length);
	} catch (Exception e) {
		// TODO: handle exception
		e.printStackTrace();
	}
	return result;

}

	public static int PredictAsetYin(String path) {
		// get latest parameter
		int result = 0;
		try {
			Database db = new Database();
			db.selectGtccParam(databasegmmyina);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();

			double[] feature = null;
			try {
				Yin yin = new Yin(16000);
				yin.main(path);
				ArrayList<Float> yins = yin.getPitchs();
				feature = new double[yins.size()];
				for (int i1 = 0; i1 < yins.size(); i1++) {
					feature[i1] = Double.valueOf(yins.get(i1).toString());
				}

			} catch (Exception e) {
				// TODO: handle exception
				// return "Audio tidak dapat diproses, dengan errorcode = " + e.toString();
			}
			GaussianMixtureModel gmm = new GaussianMixtureModel(of.means, of.covariances, of.pi, 4);
			result = gmm.predict(feature, of.means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}
	
	public static int PredictAsetYinna(String path) {
		// get latest parameter
		int result = 0;
		try {
			Database db = new Database();
			db.selectGtccParam(databasegmmyina);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();

			double[] feature = null;
			try {
				Yin yin = new Yin(16000);
				yin.main(path);
				ArrayList<Float> yins = yin.getPitchs();
				feature = new double[yins.size()];
				for (int i1 = 0; i1 < yins.size(); i1++) {
					feature[i1] = Double.valueOf(yins.get(i1).toString());
				}

			} catch (Exception e) {
				// TODO: handle exception
				// return "Audio tidak dapat diproses, dengan errorcode = " + e.toString();
			}
			GaussianMixtureModel gmm = new GaussianMixtureModel(of.means, of.covariances, of.pi, 3);
			result = gmm.predict(feature, of.means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}
	
	public static int PredictAsetYinTrain(String path) {
		// get latest parameter
		int result = 0;
		try {
/*			Database db = new Database();
			db.selectGtccParam(databasegmmyina);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();*/
			
			double [][] means = GaussianMixtureModel.means;
			double[][][] covariances = GaussianMixtureModel.covariances;
			double[] pi = GaussianMixtureModel.pi;

			double[] feature = null;
			try {
				Yin yin = new Yin(16000);
				yin.main(path);
				ArrayList<Float> yins = yin.getPitchs();
				feature = new double[yins.size()];
				for (int i1 = 0; i1 < yins.size(); i1++) {
					feature[i1] = Double.valueOf(yins.get(i1).toString());
				}

			} catch (Exception e) {
				// TODO: handle exception
				// return "Audio tidak dapat diproses, dengan errorcode = " + e.toString();
			}
			GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances, pi, 4);
			result = gmm.predict(feature, means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}
	
	public static int PredictAsetYinTrainNoArtifact(String path) {
		// get latest parameter
		int result = 0;
		try {
/*			Database db = new Database();
			db.selectGtccParam(databasegmmyina);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();*/
			
			double [][] means = GaussianMixtureModel.means;
			double[][][] covariances = GaussianMixtureModel.covariances;
			double[] pi = GaussianMixtureModel.pi;

			double[] feature = null;
			try {
				Yin yin = new Yin(16000);
				yin.main(path);
				ArrayList<Float> yins = yin.getPitchs();
				feature = new double[yins.size()];
				for (int i1 = 0; i1 < yins.size(); i1++) {
					feature[i1] = Double.valueOf(yins.get(i1).toString());
				}

			} catch (Exception e) {
				// TODO: handle exception
				// return "Audio tidak dapat diproses, dengan errorcode = " + e.toString();
			}
			GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances, pi, 3);
			result = gmm.predict(feature, means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}
	
	
	public static int PredictBsetYinTrain(String path) {
		// get latest parameter
		int result = 0;
		try {
/*			Database db = new Database();
			db.selectGtccParam(databasegmmyina);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();*/
			
			double [][] means = GaussianMixtureModel.means;
			double[][][] covariances = GaussianMixtureModel.covariances;
			double[] pi = GaussianMixtureModel.pi;

			double[] feature = null;
			try {
				Yin yin = new Yin(16000);
				yin.main(path);
				ArrayList<Float> yins = yin.getPitchs();
				feature = new double[yins.size()];
				for (int i1 = 0; i1 < yins.size(); i1++) {
					feature[i1] = Double.valueOf(yins.get(i1).toString());
				}

			} catch (Exception e) {
				// TODO: handle exception
				// return "Audio tidak dapat diproses, dengan errorcode = " + e.toString();
			}
			GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances, pi, 4);
			result = gmm.predict(feature, means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}

	public int[] PredictAsetArrays(double[][] data) {
		// get latest parameter
		int[] result = new int[161];
		try {
			Database db = new Database();
			db.selectGtccParam(databasegmma);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();

			for (int i = 0; i < data.length; i++) {
				GaussianMixtureModel gmm = new GaussianMixtureModel(of.means, of.covariances, of.pi, 4);
				result[i] = gmm.predict(data[i], of.means[0].length);
			}

		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}

	public static int PredictBset(String path) {
		// get latest parameter
		int result = 0;
		try {
			Database db = new Database();
			db.selectGtccParam(databasegmmb);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();
			double[] feature=null;
			try {
				GTCC gtcc = new GTCC();
				double[] data = StdAudio.read(path);
				feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
			} catch (Exception e) {
				// TODO: handle exception
			 /*File file = new File(path);
			 file.deleteOnExit();*/
			}
			

			GaussianMixtureModel gmm = new GaussianMixtureModel(of.means, of.covariances, of.pi, 4);
			result = gmm.predict(feature, of.means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;
	}

	public static int PredictBsetYin(String path) {
		// get latest parameter
		int result = 0;
		try {
			Database db = new Database();
			db.selectGtccParam(databasegmmyinb);
			ObjectParamGMM of;
			FileInputStream fileIn = new FileInputStream("Out.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			of = new ObjectParamGMM();
			of = (ObjectParamGMM) in.readObject();
			in.close();
			fileIn.close();

			double[] feature = null;
			try {
				Yin yin = new Yin(16000);
				yin.main(path);
				ArrayList<Float> yins = yin.getPitchs();
				feature = new double[yins.size()];
				for (int i1 = 0; i1 < yins.size(); i1++) {
					feature[i1] = Double.valueOf(yins.get(i1).toString());
				}

			} catch (Exception e) {
				// TODO: handle exception
				// return "Audio tidak dapat diproses, dengan errorcode = " + e.toString();
			}
			GaussianMixtureModel gmm = new GaussianMixtureModel(of.means, of.covariances, of.pi, 4);
			result = gmm.predict(feature, of.means[0].length);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return result;

	}

	

}

/*
 * public String TrainAset() {
 * 
 * ArrayList<ArrayList<Double>> feature = new ArrayList<>(); ArrayList<Integer>
 * label = new ArrayList<>(); double[][] means = new double[4][]; double[][][]
 * covariances = new double[4][][]; double[] phi = new double[4]; int[] sizedata
 * = new int[4]; Database db = new Database(); int[][] index = new int[4][];
 * 
 * System.out.println("Collecting data from DB"); for (int i = 0; i <
 * index.length; i++) { try { index[i] = db.selectGtccDatasetIndex(i,
 * databasea); for (int j = 0; j < index[i].length; j++) { ArrayList<Double>
 * featuretemp = new ArrayList<>(); ObjectFeatureGtcc of = new
 * ObjectFeatureGtcc(); db.selectGtccDataset(index[i][j], databasea); try {
 * FileInputStream fileIn = new FileInputStream("Out.ser"); ObjectInputStream in
 * = new ObjectInputStream(fileIn); of = new ObjectFeatureGtcc(); of =
 * (ObjectFeatureGtcc) in.readObject(); in.close(); fileIn.close(); } catch
 * (Exception e) { // TODO: handle exception } for (int k = 0; k <
 * of.data.length; k++) { featuretemp.add(of.data[k]);
 * 
 * 
 * } label.add(Integer.valueOf(of.classy));
 * System.out.println("class "+of.classy +" "+of.name);
 * feature.add(featuretemp);
 * 
 * 
 * } } catch (Exception e) { // TODO Auto-generated catch block
 * e.printStackTrace(); } // sizedata[i] = feature.get(i).size();
 * 
 * 
 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
 * feature.get(i1).size(); } }
 * 
 * 
 * }
 * 
 * // find max System.out.println("find max"); int max = feature.get(0).size();
 * for (int j = 0; j < feature.size(); j++) { if (feature.get(j).size() > max) {
 * max = feature.get(j).size(); } }
 * 
 * System.out.println("fitting with max"); // Convert to 2d array & fitting with
 * max double[][] featurearray = new double[feature.size()][max]; for (int j1 =
 * 0; j1 < feature.size(); j1++) { for (int j2 = 0; j2 < max; j2++) { if (j2 <
 * feature.get(j1).size()) { featurearray[j1][j2] = feature.get(j1).get(j2); }
 * else { featurearray[j1][j2] = 0.; }
 * 
 * } }
 * 
 * System.out.println("Reduce dimension"); RealMatrix rm =
 * MatrixUtils.createRealMatrix(featurearray); SingularValueDecomposition svd =
 * new SingularValueDecomposition(rm); RealMatrix u = svd.getU(); RealMatrix s =
 * svd.getS(); RealMatrix vt = svd.getVT(); double[][] sar = s.getData();
 * double[][] vtar = vt.getData();
 * 
 * double[][] news = new double[featurearray.length][featurearray.length];
 * double[][] newv = new double[featurearray.length][featurearray.length]; for
 * (int j = 0; j < featurearray.length; j++) { for (int j2 = 0; j2 <
 * featurearray.length; j2++) { news[j][j2] = sar[j][j2]; newv[j][j2] =
 * vtar[j][j2]; } }
 * 
 * s = MatrixUtils.createRealMatrix(news); vt =
 * MatrixUtils.createRealMatrix(newv);
 * 
 * RealMatrix featuresvd = u.multiply(s).multiply(vt); int featlenght =
 * featurearray.length; featurearray = new double[featlenght][featlenght];
 * featurearray = featuresvd.getData();
 * 
 * 
 * ArrayList<Double> datas =new ArrayList<>(); for (int i = 0; i <
 * featurearray.length; i++) { for (int j = 0; j < featurearray[i].length; j++)
 * { datas.add(featurearray[i][j]); } } double[] arraytempdata = new
 * double[datas.size()]; for (int i = 0; i < datas.size(); i++) {
 * arraytempdata[i]=datas.get(i); } DescriptiveStatistics ds = new
 * DescriptiveStatistics(arraytempdata); double mean = ds.getMean(); double
 * stdev = ds.getStandardDeviation(); double [][] newdata = new
 * double[featurearray.length][featurearray[0].length]; for (int i = 0; i <
 * featurearray.length; i++) { double[] temp = new
 * double[featurearray[i].length]; for (int j = 0; j < featurearray[0].length;
 * j++) {
 * 
 * temp[j]=(featurearray[i][j]-mean)/stdev; } newdata[i]=temp;
 * 
 * }
 * 
 * featurearray=newdata;
 * 
 * System.out.println(Arrays.deepToString(featurearray));
 * 
 * // find occurence for (int i = 0; i < index.length; i++) { int count = 0; for
 * (int j = 0; j < featurearray.length; j++) { if (label.get(j) == i) { count++;
 * } } sizedata[i] = count;
 * System.out.println("kelas "+i+" data count = "+sizedata[i]); }
 * 
 * System.out.println("Compute Means"); // compute means ArrayList<Double>
 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
 * (label.get(j2) == j) { //System.out.println(j2 + " " +
 * featurearray[j2].length + " " + sizedata[j]); System.out.println(j2);
 * meanstemp = new ArrayList<>(); double sum = 0.; for (int k = 0; k <
 * featurearray[j2].length; k++) { for (int k2 = 0; k2 < sizedata[j]; k2++) {
 * //System.out.println(k2); sum += featurearray[k2][k]; }
 * meanstemp.add(sum/sizedata[j]); }
 * 
 * }
 * 
 * } System.out.println("spaces"); meansarray = new double[meanstemp.size()];
 * for (int j1 = 0; j1 < meanstemp.size(); j1++) { meansarray[j1] =
 * meanstemp.get(j1); }
 * 
 * means[j] = meansarray; System.out.println(Arrays.toString(means[j])); }
 * 
 * System.out.println("Compute covariances"); // compute covariances for (int j
 * = 0; j < index.length; j++) { double[][] datatemp = new
 * double[sizedata[j]][]; int iter = 0; for (int j2 = 0; j2 <
 * featurearray.length; j2++) { if (label.get(j2) == j) { datatemp[iter] =
 * featurearray[j2]; iter++; } }
 * 
 * Matrix cova = new Matrix(cov.getData()); cova =
 * cova.plus(Matrix.identity(featlenght, featlenght).times(0.1)); covariances[j]
 * = cova.getArray(); System.out.println(Arrays.deepToString(covariances[j]));
 * EigenDecomposition ed = new
 * EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
 * System.out.println(ed.getDeterminant());
 * System.out.println("ke "+j+" length data = "+datatemp.length); Matrix cova;
 * cova = new Matrix(datatemp, sizedata[j], datatemp[0].length); //cova1 = new
 * Matrix(cov.getData(), sizedata[j], datatemp[0].length);
 * 
 * cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]); cova =
 * cova.plus(Matrix.identity(datatemp[0].length,
 * datatemp[0].length).times(0.1));
 * 
 * covariances[j]=cova.getArray(); for (int i = 0; i < covariances[j].length;
 * i++) { System.out.println(Arrays.toString(covariances[j][i]));
 * 
 * } System.out.println("space");
 * 
 * }
 * 
 * 
 * 
 * // compute phi double sum = 0; for (int j = 0; j < sizedata.length; j++) {
 * sum += sizedata[j]; } System.out.println("sum "+sum); for (int j = 0; j <
 * phi.length; j++) { for (int j1 = 0; j1 < sizedata.length; j1++) {
 * System.out.println("size "+sizedata[j]); phi[j] = sizedata[j] / sum; }
 * 
 * System.out.println("pi"+phi[j]); } //GaussianMixtureModel gmm = new
 * GaussianMixtureModel(featurearray, 4,means); GaussianMixtureModel gmm = new
 * GaussianMixtureModel(means, covariances, phi,4);
 * //gmm.EMKmeans1(featurearray,4,means); //gmm.EMKmeans1(featurearray,4);
 * gmm.EMKmeans1(featurearray); //gmm.fit(featurearray,4);
 * db.insertGMMParam(GaussianMixtureModel.covariances,
 * GaussianMixtureModel.means, GaussianMixtureModel.pi, databasegmma);
 * 
 * return String.valueOf(gmm.getLoglikelihoods());
 * 
 * }
 * 
 * public String TrainBset() { ArrayList<ArrayList<Double>> feature = new
 * ArrayList<>(); ArrayList<Integer> label = new ArrayList<>(); double[][] means
 * = new double[3][]; double[][][] covariances = new double[3][][]; double[] phi
 * = new double[3]; int[] sizedata = new int[3]; Database db = new Database();
 * int[][] index = new int[3][];
 * 
 * System.out.println("Collecting data from DB"); for (int i = 0; i <
 * index.length; i++) { try { index[i] = db.selectGtccDatasetIndex(i,
 * databaseb); for (int j = 0; j < index[i].length; j++) { ArrayList<Double>
 * featuretemp = new ArrayList<>(); ObjectFeatureGtcc of = new
 * ObjectFeatureGtcc(); db.selectGtccDataset(index[i][j], databaseb); try {
 * FileInputStream fileIn = new FileInputStream("Out.ser"); ObjectInputStream in
 * = new ObjectInputStream(fileIn); of = new ObjectFeatureGtcc(); of =
 * (ObjectFeatureGtcc) in.readObject(); in.close(); fileIn.close(); } catch
 * (Exception e) { // TODO: handle exception } for (int k = 0; k <
 * of.data.length; k++) { featuretemp.add(of.data[k]);
 * 
 * 
 * } label.add(Integer.valueOf(of.classy));
 * //System.out.println("class "+of.classy +" "+of.name);
 * feature.add(featuretemp);
 * 
 * 
 * } } catch (Exception e) { // TODO Auto-generated catch block
 * e.printStackTrace(); } // sizedata[i] = feature.get(i).size();
 * 
 * 
 * //find max int max = feature.get(0).size(); for (int i1 = 1; i1 <
 * feature.size(); i1++) { if (feature.get(i1).size() < max) { max =
 * feature.get(i1).size(); } }
 * 
 * 
 * }
 * 
 * // find max System.out.println("find max"); int max = feature.get(0).size();
 * for (int j = 0; j < feature.size(); j++) { if (feature.get(j).size() > max) {
 * max = feature.get(j).size(); } }
 * 
 * System.out.println("fitting with max "+max); // Convert to 2d array & fitting
 * with max double[][] featurearray = new double[feature.size()][max]; for (int
 * j1 = 0; j1 < feature.size(); j1++) { for (int j2 = 0; j2 < max; j2++) { if
 * (j2 < feature.get(j1).size()) { featurearray[j1][j2] =
 * feature.get(j1).get(j2); } else { featurearray[j1][j2] = 0.; }
 * 
 * } } int featlenght = featurearray.length;
 * System.out.println("Reduce dimension"); RealMatrix rm =
 * MatrixUtils.createRealMatrix(featurearray); SingularValueDecomposition svd =
 * new SingularValueDecomposition(rm); RealMatrix u = svd.getU(); RealMatrix s =
 * svd.getS(); RealMatrix vt = svd.getVT(); double[][] sar = s.getData();
 * double[][] vtar = vt.getData();
 * 
 * double[][] news = new double[featurearray.length][featurearray.length];
 * double[][] newv = new double[featurearray.length][featurearray.length]; for
 * (int j = 0; j < featurearray.length; j++) { for (int j2 = 0; j2 <
 * featurearray.length; j2++) { news[j][j2] = sar[j][j2]; newv[j][j2] =
 * vtar[j][j2]; } }
 * 
 * s = MatrixUtils.createRealMatrix(news); vt =
 * MatrixUtils.createRealMatrix(newv);
 * 
 * RealMatrix featuresvd = u.multiply(s).multiply(vt);
 * 
 * featurearray = new double[featlenght][featlenght]; featurearray =
 * featuresvd.getData();
 * 
 * 
 * ArrayList<Double> datas =new ArrayList<>(); for (int i = 0; i <
 * featurearray.length; i++) { for (int j = 0; j < featurearray[i].length; j++)
 * { datas.add(featurearray[i][j]); } } double[] arraytempdata = new
 * double[datas.size()]; for (int i = 0; i < datas.size(); i++) {
 * arraytempdata[i]=datas.get(i); } DescriptiveStatistics ds = new
 * DescriptiveStatistics(arraytempdata); double mean = ds.getMean(); double
 * stdev = ds.getStandardDeviation(); double [][] newdata = new
 * double[featurearray.length][featurearray[0].length]; for (int i = 0; i <
 * featurearray.length; i++) { double[] temp = new
 * double[featurearray[i].length]; for (int j = 0; j < featurearray[0].length;
 * j++) {
 * 
 * temp[j]=(featurearray[i][j]-mean)/stdev; } newdata[i]=temp;
 * 
 * }
 * 
 * featurearray=newdata;
 * 
 * //System.out.println(Arrays.deepToString(featurearray));
 * 
 * // find occurence for (int i = 0; i < index.length; i++) { int count = 0; for
 * (int j = 0; j < featurearray.length; j++) { if (label.get(j) == i) { count++;
 * } } sizedata[i] = count;
 * System.out.println("kelas "+i+" data count = "+sizedata[i]); }
 * 
 * System.out.println("Compute Means"); // compute means ArrayList<Double>
 * meanstemp = new ArrayList<>(); double[] meansarray = null; for (int j = 0; j
 * < index.length; j++) { for (int j2 = 0; j2 < featurearray.length; j2++) { if
 * (label.get(j2) == j) { //System.out.println(j2 + " " +
 * featurearray[j2].length + " " + sizedata[j]); //System.out.println(j2);
 * meanstemp = new ArrayList<>(); double sum = 0.; for (int k = 0; k <
 * featurearray[j2].length; k++) { for (int k2 = 0; k2 < sizedata[j]; k2++) {
 * //System.out.println(k2); sum += featurearray[k2][k]; }
 * meanstemp.add(sum/sizedata[j]); }
 * 
 * }
 * 
 * } //System.out.println("spaces"); meansarray = new double[meanstemp.size()];
 * for (int j1 = 0; j1 < meanstemp.size(); j1++) { meansarray[j1] =
 * meanstemp.get(j1); }
 * 
 * means[j] = meansarray; //System.out.println(Arrays.toString(means[j])); }
 * 
 * System.out.println("Compute covariances"); // compute covariances for (int j
 * = 0; j < index.length; j++) { double[][] datatemp = new
 * double[sizedata[j]][]; int iter = 0; for (int j2 = 0; j2 <
 * featurearray.length; j2++) { if (label.get(j2) == j) { datatemp[iter] =
 * featurearray[j2]; iter++; } }
 * 
 * Matrix cova = new Matrix(cov.getData()); cova =
 * cova.plus(Matrix.identity(featlenght, featlenght).times(0.1)); covariances[j]
 * = cova.getArray(); System.out.println(Arrays.deepToString(covariances[j]));
 * EigenDecomposition ed = new
 * EigenDecomposition(MatrixUtils.createRealMatrix(covariances[j]));
 * System.out.println(ed.getDeterminant());
 * //System.out.println("ke "+j+" length data = "+datatemp.length); Matrix cova;
 * cova = new Matrix(datatemp, sizedata[j], datatemp[0].length); //cova1 = new
 * Matrix(cov.getData(), sizedata[j], datatemp[0].length);
 * 
 * cova = (cova.transpose().times(cova)).times(1.0d / sizedata[j]); cova =
 * cova.plus(Matrix.identity(datatemp[0].length,
 * datatemp[0].length).times(0.1));
 * 
 * covariances[j]=cova.getArray(); for (int i = 0; i < covariances[j].length;
 * i++) { System.out.println(Arrays.toString(covariances[j][i]));
 * 
 * } System.out.println("space");
 * 
 * }
 * 
 * ArrayList<ArrayList<String>> list= new ArrayList<>(); for (int i1 = 0; i1 <
 * covariances[2].length; i1++) {
 * 
 * ArrayList<String> ap = new ArrayList<>();
 * 
 * for (int j = 0; j < covariances[2][i1].length; j++) {
 * 
 * 
 * ap.add(String.valueOf(newdata[i1][j])); }
 * 
 * list.add(ap); }
 * 
 * CSVHelper csv = new CSVHelper(); File csvfile = new
 * File("C:\\\\Users\\\\extre\\\\Desktop\\\\heart audio\\\\covariances.csv");
 * FileOutputStream fos; try { fos = new FileOutputStream(csvfile); Writer w =
 * new OutputStreamWriter(fos,"UTF-8");
 * 
 * 
 * for (ArrayList<String> h : list) { try { csv.writeLine(w, h); } catch
 * (Exception e) { // TODO Auto-generated catch block e.printStackTrace(); } }
 * 
 * w.flush(); w.close(); } catch (FileNotFoundException e1) { // TODO
 * Auto-generated catch block e1.printStackTrace(); } catch
 * (UnsupportedEncodingException e1) { // TODO Auto-generated catch block
 * e1.printStackTrace(); } catch (IOException e) { // TODO Auto-generated catch
 * block e.printStackTrace(); }
 * 
 * 
 * ArrayList<ArrayList<String>> list1= new ArrayList<>(); for (int i1 = 0; i1 <
 * featurearray.length; i1++) {
 * 
 * ArrayList<String> ap = new ArrayList<>();
 * 
 * for (int j = 0; j < featurearray[i1].length; j++) {
 * 
 * 
 * ap.add(String.valueOf(featurearray[i1][j])); }
 * 
 * list1.add(ap); }
 * 
 * CSVHelper csvs = new CSVHelper(); File csvfiles = new
 * File("C:\\\\Users\\\\extre\\\\Desktop\\\\heart audio\\\\data.csv");
 * FileOutputStream foss; try { foss = new FileOutputStream(csvfiles); Writer ww
 * = new OutputStreamWriter(foss,"UTF-8");
 * 
 * 
 * for (ArrayList<String> h : list1) { try { csvs.writeLine(ww, h); } catch
 * (Exception e) { // TODO Auto-generated catch block e.printStackTrace(); } }
 * 
 * ww.flush(); ww.close(); } catch (FileNotFoundException e1) { // TODO
 * Auto-generated catch block e1.printStackTrace(); } catch
 * (UnsupportedEncodingException e1) { // TODO Auto-generated catch block
 * e1.printStackTrace(); } catch (IOException e) { // TODO Auto-generated catch
 * block e.printStackTrace(); }
 * 
 * 
 * 
 * 
 * // compute phi double sum = 0; for (int j = 0; j < sizedata.length; j++) {
 * sum += sizedata[j]; } //System.out.println("sum "+sum); for (int j = 0; j <
 * phi.length; j++) { for (int j1 = 0; j1 < sizedata.length; j1++) {
 * //System.out.println("size "+sizedata[j]); phi[j] = sizedata[j] / sum; }
 * 
 * //System.out.println("pi"+phi[j]); } //GaussianMixtureModel gmm = new
 * GaussianMixtureModel(featurearray, 4); //GaussianMixtureModel gmm = new
 * GaussianMixtureModel(featurearray,means, covariances, phi);
 * GaussianMixtureModel gmm = new GaussianMixtureModel(means, covariances, phi,
 * 3); gmm.EMKmeans1(featurearray); //gmm.EMKmeans1(featurearray,4);
 * //gmm.EMKmeans(featurearray); //gmm.fit(featurearray,4);
 * db.insertGMMParam(GaussianMixtureModel.covariances,
 * GaussianMixtureModel.means, GaussianMixtureModel.pi, databasegmmb);
 * 
 * return String.valueOf(gmm.getLoglikelihoods()); }
 */
