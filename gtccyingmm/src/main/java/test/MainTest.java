package test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import javax.swing.JOptionPane;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import com.github.habernal.confusionmatrix.ConfusionMatrix;
import com.google.gson.Gson;
/*import com.mathworks.engine.EngineException;
import com.mathworks.engine.MatlabEngine;*/

import featureExtraction.GTCC;
import featureExtraction.GTCCMatlab;
import featureExtraction.Mfcc;
import featureExtraction.StdAudio;
import gtccyingmm.Controller;
import machineLearning.GaussianMixtureModel;
import util.DB;
import util.Database;
import util.ObjectDB;
import util.ObjectFeatureGtcc;
import yin.Yin;

public class MainTest {
	/*static String databasea = "gtccaset";
	static String databaseb = "gtccbset";
	static String databaseayin = "yinaset";
	static String databasebyin = "yinbset";*/
	

	static String databasea = "gtccasetnormalless";
	static String databaseb = "gtccbsetnormalless";
	static String databaseamat = "gtccasetmat";
	static String databasebmat = "gtccbsetnormalless";

	static String databaseayin = "yinasetnormalless";
	static String databasebyin = "yinbsetnormalless";
	static String databasegmma = "gtccgmmasetnormalless";
	static String databasegmmb = "gtccgmmbsetnormalless";
	static String databasegmmyina = "yingmmasetnormalless";
	static String databasegmmyinb = "yingmmbsetnormalless";

	static Controller c = new Controller();
/*	static GTCCMatlab gm;
public static void initMatlab() throws EngineException, IllegalArgumentException, IllegalStateException, InterruptedException {
	 gm= new GTCCMatlab();
}*/
	
public static void readProp() {
		
		Properties prop = new Properties();
		InputStream input = null;
		OutputStream output = null;

		try {
			String filename = "config.properties";
			//URI uri = Main.class.getProtectionDomain().getCodeSource().getLocation().toURI();
			//input = Main.class.getClassLoader().getResourceAsStream(filename);
			input = new FileInputStream(new File(filename));
			// load a properties file
			prop.load(input);
			ObjectDB odb = new ObjectDB();
			odb.url = prop.getProperty("database");
			odb.username = prop.getProperty("dbuser");
			odb.password = prop.getProperty("dbpassword");
			odb.webport = prop.getProperty("webport");
			odb.port = prop.getProperty("dbport");
			// get the property value and print it out
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
			
	}
	
	public static void traina() {

		c.TrainAsetWork();

	}
	
	public static void trainana() {

		c.TrainAsetWorkNoArtifact();

	}


	public static void trainayin() {

		c.TrainAsetWorkYIN();

	}
	
	public static void trainayinna() {

		c.TrainAsetWorkYINNoArtifact();

	}

	public static void trainb() {
		c.TrainBsetWork();
	}

	public static void trainbyin() {
		c.TrainBsetWorkYIN();
	}
	
	public static void toJson(double[][] featurearray,ArrayList<Integer> px) {
		Map<String, ArrayList<Integer>> jpx = new HashMap<>();
		jpx.put("px", px);
		Map<String, double[][]> jsx = new HashMap<>();
		jsx.put("x", featurearray);
		Gson gson = new Gson();
		
		try {
			Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("px.json"),"utf-8"));
			writer.write(gson.toJson(jpx));
			writer.close();
			Writer writer2 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("xtest.json"),"utf-8"));
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

	public static void predicta() {
		
		
		int truelabel = 0, falselabel = 0;
		int count = 0;
		ArrayList<Integer> res = new ArrayList<>();
		ArrayList<Integer> real = new ArrayList<>();
		GaussianMixtureModel.featuretest = new ArrayList<>();
		ConfusionMatrix cm = new ConfusionMatrix();
		System.out.println("mulai");
		for (int i = 0; i < 4; i++) {
			int nol = 0,satu = 0,dua=0,tiga=0;
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\atest\\" + i;
			//String path = "0Real\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = c.PredictAset(child.getAbsolutePath());
					res.add(result);
					real.add(i);
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
					 /*logger.info(" file " + child.getName() + " prediksi = " + result);
					 logger1.info(child.getName());
					 logger2.info(String.valueOf(result));
					 logger3.info(String.valueOf(i));*/
					//System.out.println(child.getName());
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
		System.out.println("Predict");
		for (int i = 0; i < res.size(); i++) {
			System.out.println(res.get(i));
		}
		System.out.println("Real");
		for (int i = 0; i < real.size(); i++) {
			System.out.println(real.get(i));
		}
		
		System.out.println(cm);
		System.out.println(cm.printLabelPrecRecFm());
		System.out.println(cm.getPrecisionForLabels());
		System.out.println(cm.getRecallForLabels());
		System.out.println(cm.printNiceResults());
		double[][] xtest = new double[GaussianMixtureModel.featuretest.size()][GaussianMixtureModel.featuretest.get(0).size()];
		for (int i = 0; i < GaussianMixtureModel.featuretest.size(); i++) {
			for (int j = 0; j < GaussianMixtureModel.featuretest.get(i).size(); j++) {
				xtest[i][j]=GaussianMixtureModel.featuretest.get(i).get(j);
			}
			
		}
		toJson(xtest,res);
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);

	}
	
	public static void predictana() {

		int truelabel = 0, falselabel = 0;
		int count = 0;
		for (int i = 0; i < 3; i++) {
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\atest\\" + i;
			//String path = "0Real\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = c.PredictAsetNoArtifact(child.getAbsolutePath());
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

	}

	public static void predictayisn() {
		for (int i = 0; i < 4; i++) {
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\atest\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					int result = c.PredictAsetYin(child.getAbsolutePath());
					System.out.println("folder " + i + " file " + child.getName() + " prediksi = " + result);

				}
			}
		}
	}
	
	
public static void predictayinna() {
		
		
		int truelabel = 0, falselabel = 0;
		int count = 0;
		for (int i = 0; i < 3; i++) {
			//String path = "C:\\\\Users\\\\Lenovo\\\\Documents\\\\Johan\\\\datauji\\\\atest\\\\\\" + i;
			String path = "0Real\\" + i;
			
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = c.PredictAsetYinna(child.getAbsolutePath());
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

		
	/*	
		int nol = 0, satu = 0, dua = 0;
		for (int i = 0; i < 4; i++) {
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datatrain\\b\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					int result = c.PredictBsetYin(child.getAbsolutePath());

					System.out.println("folder " + i + " file " + child.getName() + " prediksi = " + result);

				}
			}
		}*/
	}

public static void predictayin() {
	
	
	int truelabel = 0, falselabel = 0;
	int count = 0;
	ArrayList<Integer> res = new ArrayList<>();
	ArrayList<Integer> real = new ArrayList<>();
	ArrayList<Integer> list = new ArrayList<>();
	GaussianMixtureModel.featuretest = new ArrayList<>();
	ConfusionMatrix cm = new ConfusionMatrix();
	System.out.println("mulai");
	for (int i = 0; i < 4; i++) {
		int nol = 0,satu = 0,dua=0,tiga=0;
		String path = "C:\\\\Users\\\\Lenovo\\\\Documents\\\\Johan\\\\datauji\\\\atest\\\\\\" + i;
		//String path = "0Real\\" + i;
		
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			for (File child : directoryListing) {
				count++;
				int result = c.PredictAsetYin(child.getAbsolutePath());
				list.add(result);
				res.add(result);
				real.add(i);
				if (i == result) {
					truelabel++;

				} else {
					falselabel++;
				}/*
				// System.out.println(" file " + child.getName() + " prediksi = " + result);
				System.out.println(child.getName());
			}
		}
	}*/
	
	System.out.println("Predict");
	for (int i1 = 0; i1 < res.size(); i1++) {
		System.out.println(res.get(i1));
	}
	System.out.println("Real");
	for (int i1 = 0; i1 < real.size(); i1++) {
		System.out.println(real.get(i1));
	}
	/*System.out.println("number");
	for (int i1 = 1; i1 <= 52; i1++) {
		System.out.println(i1);
	}*/
	
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
				 /*logger.info(" file " + child.getName() + " prediksi = " + result);
				 logger1.info(child.getName());
				 logger2.info(String.valueOf(result));
				 logger3.info(String.valueOf(i));*/
				//System.out.println(child.getName());
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
/*	System.out.println("Predict");
	for (int i = 0; i < res.size(); i++) {
		System.out.println(res.get(i));
	}
	System.out.println("Real");
	for (int i = 0; i < real.size(); i++) {
		System.out.println(real.get(i));
	}*/
	
	System.out.println(cm);
	System.out.println(cm.printLabelPrecRecFm());
	System.out.println(cm.getPrecisionForLabels());
	System.out.println(cm.getRecallForLabels());
	System.out.println(cm.printNiceResults());			
	//Controller.toJsonPredict(list);
	double[][] xtest = new double[GaussianMixtureModel.featuretest.size()][GaussianMixtureModel.featuretest.get(0).size()];
	for (int i = 0; i < GaussianMixtureModel.featuretest.size(); i++) {
		for (int j = 0; j < GaussianMixtureModel.featuretest.get(i).size(); j++) {
			xtest[i][j]=GaussianMixtureModel.featuretest.get(i).get(j);
		}
		
	}
	toJson(xtest,res);
	double accuracy = ((double) truelabel / (double) count) * 100;
	System.out.println(
			"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);

	
/*	
	int nol = 0, satu = 0, dua = 0;
	for (int i = 0; i < 4; i++) {
		String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datatrain\\b\\" + i;
		File dir = new File(path);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			for (File child : directoryListing) {
				int result = c.PredictBsetYin(child.getAbsolutePath());

				System.out.println("folder " + i + " file " + child.getName() + " prediksi = " + result);

			}
		}
	}*/
}

	public static void predictb() {

		int truelabel = 0, falselabel = 0;
		int count = 0;
		ArrayList<Integer> list = new ArrayList<>();
		ArrayList<Integer> res = new ArrayList<>();
		ArrayList<Integer> real = new ArrayList<>();
		GaussianMixtureModel.featuretest = new ArrayList<>();
		ConfusionMatrix cm = new ConfusionMatrix();
		ArrayList<ArrayList<Integer>> confusionmatrixdata = new ArrayList<>();
		for (int i = 0; i < 3; i++) {
			int nol = 0,satu = 0,dua=0;
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\btest\\" + i;
			//String path = "1Real\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = c.PredictBset(child.getAbsolutePath());
					//System.out.println(child.getName());
					list.add(result);
					res.add(result);
					real.add(i);
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
					else {
						dua++;
					}
					
					 //System.out.println(" file " + child.getName() +" real "+i+ " prediksi = " + result);

				}
				
			}
			
			ArrayList<Integer> cmdi = new  ArrayList<>();
			cmdi.add(nol);
			cmdi.add(satu);
			cmdi.add(dua);
			confusionmatrixdata.add(cmdi);
			if (i==0) {
				//System.out.println("nol");
				cm.increaseValue("Normal", "Normal",nol);
				cm.increaseValue("Normal", "Murmur",satu);
				cm.increaseValue("Normal", "Extrasystole",dua);
				nol=0;
				satu=0;
				dua=0;
			}
			if (i==1) {
				//System.out.println("satu");
				cm.increaseValue("Murmur", "Normal",nol);
				cm.increaseValue("Murmur", "Murmur",satu);
				cm.increaseValue("Murmur", "Extrasystole",dua);
				nol=0;
				satu=0;
				dua=0;
			}
			if (i==2) {
				//System.out.println("dua");
				cm.increaseValue("Extrasystole", "Normal",nol);
				cm.increaseValue("Extrasystole", "Murmur",satu);
				cm.increaseValue("Extrasystole", "Extrasystole",dua);
				nol=0;
				satu=0;
				dua=0;
			}
			
			
			
			
		}
		
		System.out.println("Predict");
		for (int i = 0; i < res.size(); i++) {
			System.out.println(res.get(i));
		}
		/*System.out.println("Real");
		for (int i = 0; i < real.size(); i++) {
			System.out.println(real.get(i));
		}
		
		System.out.println("number");
		for (int i = 1; i <= 195; i++) {
			System.out.println(i);
		}*/
		
		System.out.println(cm);
		System.out.println(cm.printLabelPrecRecFm());
		System.out.println(cm.getPrecisionForLabels());
		System.out.println(cm.getRecallForLabels());
		System.out.println(cm.printNiceResults());
		double[][] xtest = new double[GaussianMixtureModel.featuretest.size()][GaussianMixtureModel.featuretest.get(0).size()];
		for (int i = 0; i < GaussianMixtureModel.featuretest.size(); i++) {
			for (int j = 0; j < GaussianMixtureModel.featuretest.get(i).size(); j++) {
				xtest[i][j]=GaussianMixtureModel.featuretest.get(i).get(j);
			}
			
		}
		toJson(xtest,res);
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);

	}

	
	public static void predictbtest() {
		int truelabel = 0, falselabel = 0;
		int count = 0;
		for (int i = 0; i < 3; i++) {
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\btest\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = c.PredictBset(child.getAbsolutePath());
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

	}

	public static void predictbyin() {
		
		
		int truelabel = 0, falselabel = 0;
		int count = 0;
		ArrayList<Integer> res = new ArrayList<>();
		ArrayList<Integer> real = new ArrayList<>();
		ArrayList<Integer> list = new ArrayList<>();
		GaussianMixtureModel.featuretest = new ArrayList<>();
		ConfusionMatrix cm = new ConfusionMatrix();
		for (int i = 0; i < 3; i++) {
			int nol = 0,satu = 0,dua=0,tiga=0;
			//String path = "1Real\\" + i;
			String path ="C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\btest\\"+i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					count++;
					int result = c.PredictBsetYin(child.getAbsolutePath());
					//System.out.println(child.getName());
					list.add(i);
					res.add(result);
					real.add(i);
					if (i == result) {
						truelabel++;

					} else {
						falselabel++;
					}/*
					// System.out.println(" file " + child.getName() + " prediksi = " + result);
					System.out.println(child.getName());
				}
			}
		}
		
		System.out.println("Predict");
		for (int i = 0; i < res.size(); i++) {
			System.out.println(res.get(i));
		}
		System.out.println("Real");
		for (int i = 0; i < real.size(); i++) {
			System.out.println(real.get(i));
		}
		System.out.println("number");
		for (int i = 1; i <= 52; i++) {
			System.out.println(i);
		}
		*/
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
					 /*logger.info(" file " + child.getName() + " prediksi = " + result);
					 logger1.info(child.getName());
					 logger2.info(String.valueOf(result));
					 logger3.info(String.valueOf(i));*/
					//System.out.println(child.getName());
				}
			}
			if (i==0) {
				//System.out.println("nol");
				cm.increaseValue("Normal", "Normal",nol);
				cm.increaseValue("Normal", "Murmur",satu);
				cm.increaseValue("Normal", "Extrasystole",dua);
				
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
			if (i==1) {
				//System.out.println("satu");
				cm.increaseValue("Murmur", "Normal",nol);
				cm.increaseValue("Murmur", "Murmur",satu);
				cm.increaseValue("Murmur", "Extrasystole",dua);
				
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
			if (i==2) {
				//System.out.println("dua");
				cm.increaseValue("Extrasystole", "Normal",nol);
				cm.increaseValue("Extrasystole", "Murmur",satu);
				cm.increaseValue("Extrasystole", "Extrasystole",dua);
				
				nol=0;
				satu=0;
				dua=0;
				tiga=0;
			}
			
		}
	/*	System.out.println("Predict");
		for (int i = 0; i < res.size(); i++) {
			System.out.println(res.get(i));
		}
		System.out.println("Real");
		for (int i = 0; i < real.size(); i++) {
			System.out.println(real.get(i));
		}*/
		
		System.out.println("Predict");
		for (int i = 0; i < res.size(); i++) {
			System.out.println(res.get(i));
		}
		/*System.out.println("Real");
		for (int i = 0; i < real.size(); i++) {
			System.out.println(real.get(i));
		}
		
		System.out.println("number");
		for (int i = 1; i <= 195; i++) {
			System.out.println(i);
		}
		*/
		
		System.out.println(cm);
		System.out.println(cm.printLabelPrecRecFm());
		System.out.println(cm.getPrecisionForLabels());
		System.out.println(cm.getRecallForLabels());
		System.out.println(cm.printNiceResults());	
		
		double[][] xtest = new double[GaussianMixtureModel.featuretest.size()][GaussianMixtureModel.featuretest.get(0).size()];
		for (int i = 0; i < GaussianMixtureModel.featuretest.size(); i++) {
			for (int j = 0; j < GaussianMixtureModel.featuretest.get(i).size(); j++) {
				xtest[i][j]=GaussianMixtureModel.featuretest.get(i).get(j);
			}
			
		}
		toJson(xtest,res);
		
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);

		//Controller.toJsonPredict(list);
	/*	
		int nol = 0, satu = 0, dua = 0;
		for (int i = 0; i < 4; i++) {
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datatrain\\b\\" + i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			if (directoryListing != null) {
				for (File child : directoryListing) {
					int result = c.PredictBsetYin(child.getAbsolutePath());

					System.out.println("folder " + i + " file " + child.getName() + " prediksi = " + result);

				}
			}
		}*/
	}
	
	public static void insertb() {
		Database db = new Database();
		try {
			db.truncateDb(databaseb);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		for (int i = 0; i < 3; i++) {
			// TODO Auto-generated method stub
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\btest\\" + i;
			//String path = "1Real\\"+i;
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\b\\"+i;
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\b\\pisah\\"+i;
			//String path = "150\\"+i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			//double[][] fctore = new double[directoryListing.length][];
			int iter = 0;
			if (directoryListing != null) {
				for (File child : directoryListing) {
					/*double[] feature = null;
				try {
						Yin yin = new Yin(16000);
						yin.main(child.getAbsolutePath());
						ArrayList<Float> yins= yin.getPitchs();
						feature = new double[yins.size()];
						for (int i1 = 0; i1 < yins.size(); i1++) {
							feature[i1] = Double.valueOf(yins.get(i1).toString());
						}*/


					double[] feature = null;
					try {
						GTCC gtcc = new GTCC();
						double[] data = StdAudio.read(child.getAbsolutePath());
/*						for (int j = 0; j < data.length; j++) {
							data[j]+=10;
						}*/
						feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
					} catch (Exception e) {
						// TODO: handle exception
						System.out.println(child.getName()+"  hapus dari kelas "+i);
						child.deleteOnExit();
						e.printStackTrace();				
						}

					
					try {
						//Database db = new Database();
						db.insertGtccDataset(String.valueOf(i), feature, child.getName(), databaseb);
					} catch (Exception e) {
						// TODO: handle exception
						//return "Gagal insert feature ke database, dengan errorcode = " + e.toString();
					}

				}
			}
		}
	}
	
	public static void inserta() throws  IllegalArgumentException, IllegalStateException, InterruptedException {
		Database db = new Database();
		try {
			db.truncateDb(databasea);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		for (int i = 0; i < 4; i++) {
			// TODO Auto-generated method stub
			//String path = "0Real\\"+i;
			//String path = "0Norm\\"+i;
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\a\\"+i;
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\a\\pisah\\"+i;
			
			
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\atest\\"+i;
			
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			//double[][] fctore = new double[directoryListing.length][];
			
			
			
			
			int iter = 0;
			if (directoryListing != null) {
				for (File child : directoryListing) {
					/*double[] feature = null;
				try {
						Yin yin = new Yin(16000);
						yin.main(child.getAbsolutePath());
						ArrayList<Float> yins= yin.getPitchs();
						feature = new double[yins.size()];
						for (int i1 = 0; i1 < yins.size(); i1++) {
							feature[i1] = Double.valueOf(yins.get(i1).toString());
						}*/


					double[] feature = null;
					try {
						GTCC gtcc = new GTCC();
						double[] data = StdAudio.read(child.getAbsolutePath());
/*						for (int j = 0; j < data.length; j++) {
							data[j]+=10;
						}*/
						//feature = gtcc.GetFeatureVector(data, 0.9, 22000, 16000);
						feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
						
					} catch (Exception e) {
						// TODO: handle exception
						//child.deleteOnExit();
						e.printStackTrace();				}

					
					try {
						
						
						db.insertGtccDataset(String.valueOf(i), feature, child.getName(), databasea);
					} catch (Exception e) {
						// TODO: handle exception
						//return "Gagal insert feature ke database, dengan errorcode = " + e.toString();
					}

				}
			}
		}
	}
	
	
	/*public static void insertamat() throws EngineException, IllegalArgumentException, IllegalStateException, InterruptedException {
		for (int i = 0; i < 4; i++) {
			// TODO Auto-generated method stub
			String path = "0Real\\"+i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			//double[][] fctore = new double[directoryListing.length][];
			
			
			
			
			int iter = 0;
			if (directoryListing != null) {
				for (File child : directoryListing) {
					double[] feature = null;
				try {
						Yin yin = new Yin(16000);
						yin.main(child.getAbsolutePath());
						ArrayList<Float> yins= yin.getPitchs();
						feature = new double[yins.size()];
						for (int i1 = 0; i1 < yins.size(); i1++) {
							feature[i1] = Double.valueOf(yins.get(i1).toString());
						}


					double[] feature = null;
					try {
						//GTCC gtcc = new GTCC();
						
						//double[] data = StdAudio.read(child.getAbsolutePath());
						for (int j = 0; j < data.length; j++) {
							data[j]+=10;
						}
						//feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
						feature = gm.extractFeature(child.getAbsolutePath());
					} catch (Exception e) {
						// TODO: handle exception
						child.deleteOnExit();
						e.printStackTrace();				}

					
					try {
						Database db = new Database();
						db.insertGtccDataset(String.valueOf(i), feature, child.getName(), databaseamat);
					} catch (Exception e) {
						// TODO: handle exception
						//return "Gagal insert feature ke database, dengan errorcode = " + e.toString();
					}

				}
			}
		}
	}
	
	public static void insertbyin() {
		Database db = new Database();
		try {
			db.truncateDb(databasebyin);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		for (int i = 0; i < 3; i++) {
			// TODO Auto-generated method stub
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\btest\\" + i;
			//String path = "1Real\\"+i;
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\b\\"+i;
			//String path = "1\\"+i;
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\b\\pisah\\"+i;
			
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			//double[][] fctore = new double[directoryListing.length][];
			int iter = 0;
			if (directoryListing != null) {
				for (File child : directoryListing) {
					double[] feature = null;
				try {
						Yin yin = new Yin(16000);
						yin.main(child.getAbsolutePath());
						ArrayList<Float> yins= yin.getPitchs();
						feature = new double[yins.size()];
						for (int i1 = 0; i1 < yins.size(); i1++) {
							feature[i1] = Double.valueOf(yins.get(i1).toString());
						}
				}


					double[] feature = null;
					try {
						GTCC gtcc = new GTCC();
						double[] data = StdAudio.read(child.getAbsolutePath());
						for (int j = 0; j < data.length; j++) {
							data[j]+=10;
						}
						feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
					} catch (Exception e) {
						// TODO: handle exception
						//child.deleteOnExit();
						e.printStackTrace();				}

					
					try {
						//Database db = new Database();
						db.insertGtccDataset(String.valueOf(i), feature, child.getName(), databasebyin);
					} catch (Exception e) {
						// TODO: handle exception
						//return "Gagal insert feature ke database, dengan errorcode = " + e.toString();
					}

				}
			}
		}
	}*/
	
	public static void insertayin() {
		Database db = new Database();
		try {
			db.truncateDb(databaseayin);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		for (int i = 0; i < 4; i++) {
			// TODO Auto-generated method stub
			//String path = "0\\"+i;
			//String path = "0Real\\"+i;
			//String path = "0Norm\\"+i;
			String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\a\\pisah\\"+i;
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\thereal\\a\\"+i;
			
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			//double[][] fctore = new double[directoryListing.length][];
			int iter = 0;
			if (directoryListing != null) {
				for (File child : directoryListing) {
					double[] feature = null;
				try {
						Yin yin = new Yin(16000);
						yin.main(child.getAbsolutePath());
						ArrayList<Float> yins= yin.getPitchs();
						feature = new double[yins.size()];
						for (int i1 = 0; i1 < yins.size(); i1++) {
							feature[i1] = Double.valueOf(yins.get(i1).toString());
						}


					/*double[] feature = null;
					try {
						GTCC gtcc = new GTCC();
						double[] data = StdAudio.read(child.getAbsolutePath());
						for (int j = 0; j < data.length; j++) {
							data[j]+=10;
						}
						feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);*/
					} catch (Exception e) {
						// TODO: handle exception
						child.deleteOnExit();
						e.printStackTrace();				}

					
					try {
						//Database db = new Database();
						db.insertGtccDataset(String.valueOf(i), feature, child.getName(), databaseayin);
					} catch (Exception e) {
						// TODO: handle exception
						//return "Gagal insert feature ke database, dengan errorcode = " + e.toString();
					}

				}
			}
		}
	}

	public static void insertayinna() {
		for (int i = 0; i < 3; i++) {
			// TODO Auto-generated method stub
			String path = "0Real\\"+i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			//double[][] fctore = new double[directoryListing.length][];
			int iter = 0;
			if (directoryListing != null) {
				for (File child : directoryListing) {
					double[] feature = null;
				try {
						Yin yin = new Yin(16000);
						yin.main(child.getAbsolutePath());
						ArrayList<Float> yins= yin.getPitchs();
						feature = new double[yins.size()];
						for (int i1 = 0; i1 < yins.size(); i1++) {
							feature[i1] = Double.valueOf(yins.get(i1).toString());
						}


					/*double[] feature = null;
					try {
						GTCC gtcc = new GTCC();
						double[] data = StdAudio.read(child.getAbsolutePath());
						for (int j = 0; j < data.length; j++) {
							data[j]+=10;
						}
						feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);*/
					} catch (Exception e) {
						// TODO: handle exception
						child.deleteOnExit();
						e.printStackTrace();				}

					
					try {
						Database db = new Database();
						db.insertGtccDataset(String.valueOf(i), feature, child.getName(), databaseayin);
					} catch (Exception e) {
						// TODO: handle exception
						//return "Gagal insert feature ke database, dengan errorcode = " + e.toString();
					}

				}
			}
		}
	}
	
	
	
	public static void getBestFeature() {
		Database db = new Database();
		int truelabel = 0, falselabel = 0;
		int countindexdb = 1;
		int count = 0;
		try {
			count = db.getCount(databasea);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		while(countindexdb<=count) {
		for (int i = 0; i < 4; i++) {
			//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\atest\\" + i;
			String pathtemp ="temptest\\";
			String pathreal = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji1\\atest\\"+i;
			
			File dir = new File(pathreal);
			File[] directoryListing = dir.listFiles();
			String asal="";
			String tujuan="";
			if (directoryListing != null) {
				for (File child : directoryListing) {
					int result = c.PredictAset(child.getAbsolutePath());
					if (i == result) {
						truelabel++;

					} else {
						falselabel++;
						System.out.println(child.getAbsolutePath()+" salah");
						double[] feature = null;
						//String wrong = child.getName();
						GTCC gtcc = new GTCC();
						double[] data = StdAudio.read(child.getAbsolutePath());
						StdAudio.close();
						feature = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
						db.insertGtccDataset(String.valueOf(i), feature, child.getName(), databasea);
						try {
							
							System.out.println("pindah "+child.getAbsolutePath()+" ke "+pathtemp+child.getName());
							asal = child.getAbsolutePath();
							tujuan = pathtemp+child.getName();
							Files.copy(Paths.get(asal), Paths.get(tujuan), StandardCopyOption.REPLACE_EXISTING);
							//JOptionPane.showMessageDialog(null, "ayo apus "+asal);
							child.delete();
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						
						

							ObjectFeatureGtcc of = new ObjectFeatureGtcc();
							try {
								db.selectGtccDataset(countindexdb, databasea);
								FileInputStream fileIn = new FileInputStream("Out.ser");
								ObjectInputStream in = new ObjectInputStream(fileIn);
								of = new ObjectFeatureGtcc();
								of = (ObjectFeatureGtcc) in.readObject();
								in.close();
								fileIn.close();
							} catch (Exception e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
							String asals="";
							String tujuans="";

							System.out.println(of.classy + " "+of.name+" inivroh");
							String path = "0Real1\\" + of.classy;
							File dirs = new File(path);	
							File[] directoryListingTrain = dirs.listFiles();
							if (directoryListing != null) {
								for (File childs : directoryListingTrain) {
									System.out.println(childs.getName());
									if (childs.getName().equals(of.name)) {
										try {
											String pathtest = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji1\\atest\\"+of.classy+"\\";
											System.out.println(childs.getAbsolutePath()+" dari db masuk ke= "+pathtest+childs.getName());
											asals = childs.getAbsolutePath();
											tujuans = pathtest+childs.getName();
											Files.copy(Paths.get(asals), Paths.get(tujuans), StandardCopyOption.REPLACE_EXISTING);
											//childs.delete();
											//JOptionPane.showMessageDialog(null, "ayo pindahin "+asal+ " ke "+tujuan);
										} catch (IOException e) {
											// TODO Auto-generated catch block
											e.printStackTrace();
										}
									}
								}
							}
							System.out.println("index "+countindexdb+" hapus dari db");
							try {
								db.deleteData(countindexdb,databasea);
							} catch (Exception e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
							
							countindexdb++;
							System.out.println("Train Ulang");
							System.out.println(i+" sekarang");
							try {
								traina();
								i=0;
							} catch (Exception e) {
								// TODO: handle exception
								System.err.println("nah gagal kan trainnya , balikin lagi");
								try {
									Files.copy(Paths.get(tujuan), Paths.get(asal), StandardCopyOption.REPLACE_EXISTING);
									File file = new File(tujuan);
									file.delete();
									GTCC gtccs = new GTCC();
									double[] datas = StdAudio.read(asal);
									StdAudio.close();
									feature = gtccs.GetFeatureVector(datas, 0.9, 4000, 1600);
									db.insertGtccDataset(String.valueOf(i), feature, child.getName(), databasea);
									Files.copy(Paths.get(tujuans),Paths.get(asals), StandardCopyOption.REPLACE_EXISTING);
									File file1 = new File(tujuans);
									file1.delete();
									
									
								} catch (IOException e1) {
									// TODO Auto-generated catch block
									e1.printStackTrace();
								}
								
							}
							
							
							
					}

				}
			}
			if (i==4) {
				break;
			}
		}
		}
		double accuracy = ((double) truelabel / (double) count) * 100;
		System.out.println(
				"Jumlah Prediksi benar = " + truelabel + " ,salah= " + falselabel + " dengan akurasi = " + accuracy);
	}

	public static void main(String[] args) throws Exception {
		//GTCCMatlab gm = new GTCCMatlab();
		//gm.extractFeature("C:\\Users\\Lenovo\\Documents\\Johan\\sk1\\heartwebnew.zip_expanded\\gtccyingmm\\0\\0\\201101070538.wav");
       /* initMatlab();
		insertamat();
		gm.closeMatlab()*/;
		
		
		//getBestFeature();

		readProp();
		Controller con = new Controller();
		//con.TrainAsetWeb();
		//trainayin();
		con.TrainAsetWebYIN();
		con.TrainBsetWebYIN();

		
		//#GTCC A
		//inserta();
		//traina();
		//predicta();
		
		//#GTCC A NA
		//inserta();
		//traina();
		//predicta();
						
		//#GTCC B
		//insertb();
		//trainb();
		
	/*	Database db = new Database();
		db.insertTerpaksa();*/
		//predictb();
		
		//#YIN A
		//insertayin();
		//trainayin();
		//predictayin();
		
		//#GTCC A NA
		//inserta();
		//traina();
		//predicta();
		
		//#YIN B
		//insertbyin();
		//trainbyin();
		//predictbyin();
		
		
		
		/*GTCC gtcc = new GTCC();
		double[] data = StdAudio.read("C:\\Users\\Lenovo\\Documents\\Johan\\sk1\\heartwebnew.zip_expanded\\gtccyingmm\\0\\0\\201101070538.wav");
		gtcc.GetFeatureVector(data, 0.9, 4000, 1600);*/
		//String path = "C:\\Users\\Lenovo\\Documents\\Johan\\datauji\\atest\\"+i;
	///trainana();
		
		
	 
	 //predictbtest();
	 
	//predictana();
	
		
	//predictbyin();
	//trainayinna();
	 
	 //predictayinna();
	 
		 
	// System.out.println(c.PredictAset("C:\\Users\\Lenovo\\Documents\\Johan\\Bunlabelledtest\\Bunlabelledtest\\135_1306428972976_D.wav"));

/*		int[] nc = new int[4];
		int[] co = new int[4];
		//ArrayList<ArrayList<Double>> featurearrayold = new ArrayList<>();
		double [][] featureold = new double[117][];
		int iter=0;
		ArrayList<ArrayList<Double>> feature = new ArrayList<>();
		for (int i = 0; i < 3; i++) {
			String path = "0\\"+i;
			File dir = new File(path);
			File[] directoryListing = dir.listFiles();
			int notcor = 0;
			int cor = 0;
			
			if (directoryListing != null) {
				for (File child : directoryListing) {
					GTCC gtcc = new GTCC();
					ArrayList<Double> f = new ArrayList<>();
					double []data = StdAudio.read(child.getAbsolutePath());
					 featureold[iter] = gtcc.GetFeatureVector(data, 0.9, 4000, 1600);
					
					for (int j = 0; j < featureold[iter].length; j++) {
						f.add(featureold[iter][j]);
					}
					feature.add(f);
					iter++;
				}
			}
		}
			
			
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
			
			int [] result = c.PredictAsetArrays(featurearray);
				
					
			int notcor = 0;
			int cor = 0;
			int j =0;
			for (int i = 0; i < 4; i++) {
				String path = "0\\"+i;
				File dir = new File(path);
				File[] directoryListing = dir.listFiles();
				
				
				if (directoryListing != null) {
					for (File child : directoryListing) {
						if (i == result[j]) {
							cor++;

						} else {
							notcor++;
						}
						System.out.println("folder "+i+" file "+child.getName()+" prediksi = "+result[j]);
						j++;
					}
				}
			}
			
			double accuracy = ((double) cor / (double) j) * 100;
			System.out.println(
					"Jumlah Prediksi benar = " + cor + " ,salah= " + notcor + " dengan akurasi = " + accuracy);
			
			
	 
*/
}

}
