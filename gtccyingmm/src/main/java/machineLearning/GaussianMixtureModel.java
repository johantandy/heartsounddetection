package machineLearning;

import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.apache.commons.math3.distribution.MixtureMultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization;
import org.apache.commons.math3.fitting.leastsquares.GaussNewtonOptimizer.Decomposition;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

import Jama.LUDecomposition;
import Jama.Matrix;
import gtccyingmm.Main;
import plot.XyChart;
import util.ObjectParamGMM;

public class GaussianMixtureModel {
	private static KMeans km;
	private ArrayList<ArrayList<Integer>> ClusterIndex;
	public static ArrayList<ArrayList<Double>> featuretest = new ArrayList<>();
	public static double[][] means;
	public static double[] pi;
	public static double[][][] covariances;
	private int[] label;
	private int numCLuster;
	private int maximumiteration = 1000;
	private double[][] data;
	private double loglikelihoods = 0;

	public double getLoglikelihoods() {
		return loglikelihoods;
	}

	public double[][] getData() {
		return data;
	}

/*	public GaussianMixtureModel(double[][] data,double[][] means, double[][][] covariances, double[] pi) {
		// this.data = data;
		this.data=data;
		this.covariances = covariances;
		this.means = means;
		this.pi = pi;
		//initKmeans(this.data, 4);
		// EM(this.data,numClusters);
		EMApache(this.data);

	}*/
	
	//CONSTRUCTOR FOR EM APACHE & JOHAN ONLY
	public GaussianMixtureModel(double[][] means, double[][][] covariances, double[] pi,int numCluster) {
		this.numCLuster=numCluster;
		this.covariances = covariances;
		this.means = means;
		this.pi = pi;
	}
	
	//CONSTRUCTOR FOR EM KMEANS APACHE & JOHAN ONLY
		public GaussianMixtureModel(double[][] data, int numClusters, double[][] means) {
			this.numCLuster=numClusters;
			this.covariances = covariances;
			this.means = means;
			this.pi = pi;
			initKmeans(data, numClusters);
		}
	
	
	/*//ini yang dipake constructornya
	public GaussianMixtureModel(double[][] means, double[][][] covariances, double[] pi,int numCluster) {
		// this.data = data;
		this.covariances = covariances;
		this.means = means;
		this.pi = pi;
		this.numCLuster=numCluster;
		
		//initKmeans(this.data, numCluster,means);
		// EM(this.data,numClusters);
		// EMApache(this.data);

	}*/

/*	public GaussianMixtureModel(double[][] data, int numClusters, double[][] means) {
		this.data = data;
		this.means=means;
		this.numCLuster=numClusters;
		initKmeans(data, numClusters);
		EMApache(data);
		// System.out.println(Arrays.deepToString(covariances));
		//EMKmeans(this.data);
		//EMKmeansWork(this.data);
		// fit(this.data, 4);
		// predict(this.data);
	}*/

	public void initKmeans(double[][] data, int numClusters) {
		km = new KMeans(data);
		/*
		 * String fileName = "C:\\Users\\extre\\Downloads\\iris.csv"; String filename1 =
		 * "C:\\Users\\extre\\Desktop\\heart audio\\heartgtcc.csv"; km = new
		 * KMeans(fileName, null);
		 */

		// Sementara IRIS Dataset
		/*
		 * String fileName = "C:\\Users\\extre\\Downloads\\iris.csv"; km = new
		 * KMeanss(fileName, null);
		 */

	/*	double[] datas = new double[km.get_data().length * km.get_data()[0].length];
		int iter1 = 0;
		for (int i = 0; i < km.get_data().length; i++) {

			for (int j = 0; j < km.get_data()[i].length; j++) {

				datas[iter1] = km.get_data()[i][j];
				iter1++;
			}
		}*/

		/*
		 * DescriptiveStatistics ds = new DescriptiveStatistics(datas); double mean =
		 * ds.getMean(); double stdev = ds.getStandardDeviation(); double [][] newdata =
		 * new double[km.get_data().length][km.get_data()[0].length]; for (int i = 0; i
		 * < km.get_data().length; i++) { double[] temp = new
		 * double[km.get_data()[i].length]; for (int j = 0; j < km.get_data()[i].length;
		 * j++) {
		 * 
		 * temp[j]=(km.get_data()[i][j]-mean)/stdev; } newdata[i]=temp;
		 * 
		 * }
		 * 
		 * km.set_data(newdata);
		 */

		km.clustering(numClusters, this.maximumiteration,means);
		this.data = km.get_data();
		means = km.get_centroids();
		label = km.get_label();
		ClusterIndex = new ArrayList<>();
		km.printResults();

		for (int i = 0; i < numClusters; i++) {
			ArrayList<Integer> temp = new ArrayList<>();
			for (int j = 0; j < label.length; j++) {
				if (label[j] == i) {
					temp.add(j);

				}

			}
			ClusterIndex.add(temp);
		}

		covariances = computeCovariances(data);

		double[] numPoint = new double[ClusterIndex.size()];
		for (int i = 0; i < ClusterIndex.size(); i++) {
			numPoint[i] = (double) (ClusterIndex.get(i).size()) / (double) this.data.length;
		}

		pi = numPoint;
	}

	public void initSemiSupervised(double[][] data, int[] label) {

	}

	public double[][][] computeCovariances(double[][] data) {
		double[][][] covariances = new double[this.ClusterIndex.size()][this.ClusterIndex.get(0)
				.size()][this.ClusterIndex.size()];
		Matrix[] cova = new Matrix[this.ClusterIndex.size()];
		;

		for (int i = 0; i < this.ClusterIndex.size(); i++) {

			double[][] tempdata = new double[this.ClusterIndex.get(i).size()][this.ClusterIndex.size()];

			for (int j = 0; j < this.ClusterIndex.get(i).size(); j++) {
				tempdata[j] = data[this.ClusterIndex.get(i).get(j)];
			}

			System.out.println("data length " + tempdata.length);
			cova[i] = new Matrix(tempdata, ClusterIndex.get(i).size(), data[0].length);
			cova[i] = (cova[i].transpose().times(cova[i])).times(1.0d / ClusterIndex.get(i).size());

			if (cova[i].rank() < data[0].length) {
				// Arbitrary width used if variance collapses to zero: make it 'large' so
				// that centre is responsible for a reasonable number of points.
				double GMM_WIDTH = 1.0d;

				// add GMM_WIDTH*Identity to rank-deficient covariance matrices
				cova[i] = cova[i].plus(Matrix.identity(data[0].length, data[0].length).times(0.1));
				System.out.println("rank " + cova[i].rank());
			}

			/*
			 * RealMatrix mx = MatrixUtils.createRealMatrix(tempdata); RealMatrix cov = new
			 * Covariance(mx).getCovarianceMatrix();
			 */

			/*
			 * SimpleMatrix X = new SimpleMatrix(tempdata); int n = X.numRows();
			 * SimpleMatrix Xt = X.transpose(); int m = Xt.numRows();
			 * 
			 * SimpleMatrix x = new SimpleMatrix(m, 1); for(int r=0; r<m; r++ ){ x.set(r, 0,
			 * Xt.extractVector(true, r).elementSum() / n); }
			 * 
			 * SimpleMatrix S = new SimpleMatrix(m, m); for(int r=0; r<m; r++){ for(int c=0;
			 * c<m; c++){ if(r > c){ S.set(r, c, S.get(c, r)); } else { double cov1 =
			 * Xt.extractVector(true, r).minus( x.get((r), 0) ).dot(Xt.extractVector(true,
			 * c).minus( x.get((c), 0) ).transpose()); S.set(r, c, (cov1 / n)); } } }
			 * 
			 * for (int j = 0; j < tempdata.length; j++) { double [] ds = new
			 * double[tempdata.length]; for (int j2 = 0; j2 < tempdata.length; j2++) {
			 * 
			 * ds[j2]=S.get(j, j2); } tempdata[i]=ds;
			 * 
			 * }
			 */
			/*
			 * RealMatrix b = MatrixUtils.createRealMatrix(tempdata);
			 * SingularValueDecomposition svd = new SingularValueDecomposition(b);
			 * //CholeskyDecomposition cd = new CholeskyDecomposition(b); //RealMatrix a =
			 * new SingularValueDecomposition(b).getCovariance(0); DecompositionSolver l =
			 * new org.apache.commons.math3.linear.LUDecomposition(b).getSolver();
			 * RealMatrix cov = new Covariance(b).getCovarianceMatrix();
			 * tempdata=cov.getData();
			 */

			/*
			 * EigenDecomposition ed = new EigenDecomposition(cov);
			 * System.out.println("deter "+ed.getDeterminant());
			 */
			// covariances[i] = cov.getData();
			covariances[i] = cova[i].getArray();
			System.out.println("panjang covar" + covariances[i].length);
		}

		return covariances;
	}

	public double log2(double num) {
		if (num == 0)
			return 0;
		else
			return (Math.log(num) / Math.log(2));
	}
	
	
	public void EMKmeansWork(double[][] data) {
		Main.status="Perform Train process ... ";
		double TotalLoglikelihood = 0;
		double ChangeInLogLikelihood = 1;
		//for (int ites = 0; ites < 100; ites++) {
		
		while (ChangeInLogLikelihood > 0.0001) {

			double LogLikelihood = 0, ProbOfXn = 0;
			double[][] Rspb = new double[data.length][this.numCLuster];
			for (int i = 0; i < data.length; i++) {
				ProbOfXn = 0;
				System.err.println(i);
				for (int j = 0; j < this.numCLuster; j++) {

					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					// 
					Rspb[i][j] = pi[j] * mnd.density(data[i]);
					ProbOfXn = ProbOfXn + pi[j] * mnd.density(data[i]);

				}

				LogLikelihood = LogLikelihood + log2(ProbOfXn);
			}

			// E-step
			
			Main.status="Performing Expectation Maximization ... ";
			for (int i = 0; i < data.length; i++) {
				ProbOfXn = 0;
				for (int j = 0; j < this.numCLuster; j++) {
					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					
					System.out.println(i);
					ProbOfXn = ProbOfXn + Rspb[i][j];
				}

				for (int j = 0; j < this.numCLuster; j++) {
					Rspb[i][j] = Rspb[i][j] / ProbOfXn;
				}

			}

			// M-Step
			double[] N_k = new double[this.numCLuster];

			for (int k = 0; k < this.numCLuster; k++) { // Calculating N_k's
				for (int n = 0; n < data.length; n++) {
					N_k[k] = N_k[k] + Rspb[n][k];
				}
			}

			// Update Means
			for (int i = 0; i < means.length; i++) {
				for (int j = 0; j < means[0].length; j++) {
					means[i][j] = 0;
				}

				double sumposterior = 0;
				for (int i1 = 0; i1 < data.length; i1++) {
					double pw = Rspb[i1][i];
					sumposterior += pw;
					for (int j = 0; j < means[0].length; j++) {
						means[i][j] += (pw * data[i1][j]);
					}
				}

				double oosumposterior = 1. / sumposterior;
				for (int i1 = 0; i1 < means[0].length; i1++) {
					means[i][i1] *= oosumposterior;
				}
				
				System.out.println(Arrays.toString(means[i]));
			}

			// Update Covariances
			for (int i = 0; i < means.length; i++) {
				double[][] covv = covariances[i];
				double[] mean = means[i];

				for (int j = 0; j < means[0].length; j++) {
					for (int j2 = 0; j2 < means[0].length; j2++) {
						covv[j][j2] = 0.;
					}
				}

				double sumposterior = 0.;
				for (int l = 0; l < data.length; l++) {
					double[] d = data[l];
					double pw = Rspb[l][i];
					sumposterior += pw;
					for (int k = 0; k < means[0].length; k++) {
						for (int j = 0; j < means[0].length; j++) {
							covv[k][j] += (pw * (d[k] - mean[k]) * (d[j] - mean[j]));
						}

					}

				}

				double oosumposterior = 1. / sumposterior;
				for (int j = 0; j < means[0].length; j++) {
					for (int j2 = 0; j2 < means[0].length; j2++) {
						covv[j][j2] *= oosumposterior;
					}
				}
				
				Matrix cova = new Matrix(covv);
				cova = cova.plus(Matrix.identity(covv.length, covv.length).times(0.1));
				
				covariances[i] = cova.getArray();
				System.out.println(Arrays.deepToString(covariances[i]));

				//covariances[i] = covv;

			}

			// Update Pi
			for (int i = 0; i < means.length; i++) {
				double sum = 0;
				for (int j = 0; j < data.length; j++) {
					sum += Rspb[j][i];
				}
				pi[i] = sum / data.length;
			}

			double NewLogLikelihood = 0, ProbOfX = 0;

			for (int i = 0; i < data.length; i++) {
				ProbOfX = 0;
				for (int j = 0; j < this.numCLuster; j++) {
					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					System.out.println(i);
					ProbOfX = ProbOfX + pi[j] * mnd.density(data[i]);

				}
				NewLogLikelihood = NewLogLikelihood + log2(ProbOfX);
			}
			Main.status="Loglikelihood = "+ChangeInLogLikelihood;
			ChangeInLogLikelihood = NewLogLikelihood - LogLikelihood;
			TotalLoglikelihood = NewLogLikelihood;
			System.out.println("Total LogLikelihood: " + TotalLoglikelihood);
			System.out.println("Change in LogLikelihood: " + ChangeInLogLikelihood);

		}
	}

	public void fit(double[][] data, int numOfCluster) {
		this.data = data;
		double TotalLoglikelihood = 0;
		double ChangeInLogLikelihood = 1;
		double[][] Rspb = new double[data.length][numOfCluster];
		while (ChangeInLogLikelihood > 0.0001) {

			double LogLikelihood = 0, ProbOfXn = 0;
			for (int i = 0; i < data.length; i++) {
				ProbOfXn = 0;

				for (int j = 0; j < numOfCluster; j++) {

					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					System.out.println(j + "   " + mnd.density(data[i]));
					// System.out.println(mnd.density(data[i]));
					// Rspb[i][j] = pi[j] * mnd.density(data[i]);
					ProbOfXn = ProbOfXn + pi[j] * mnd.density(data[i]);
					// System.err.println(i);

				}

				LogLikelihood = LogLikelihood + log2(ProbOfXn);
			}

			System.out.println("E-STEP");
			// E-step

			for (int i = 0; i < data.length; i++) {
				ProbOfXn = 0;
				for (int j = 0; j < numOfCluster; j++) {
					// MultivariateNormalDistribution mnd = new
					// MultivariateNormalDistribution(means[j], covariances[j]);
					// System.err.println(i);

					ProbOfXn = ProbOfXn + Rspb[i][j];
				}

				for (int j = 0; j < numOfCluster; j++) {
					Rspb[i][j] = Rspb[i][j] / ProbOfXn;

				}

			}

			// M-Step
			double[] N_k = new double[numOfCluster];

			for (int k = 0; k < numOfCluster; k++) { // Calculating N_k's
				for (int n = 0; n < data.length; n++) {
					N_k[k] = N_k[k] + Rspb[n][k];
				}
			}

			// Update Means
			for (int i = 0; i < means.length; i++) {
				for (int j = 0; j < means[0].length; j++) {
					means[i][j] = 0;
				}

				double sumposterior = 0;
				for (int i1 = 0; i1 < data.length; i1++) {
					double pw = Rspb[i1][i];
					sumposterior += pw;
					for (int j = 0; j < means[0].length; j++) {
						means[i][j] += (pw * data[i1][j]);

					}
				}

				double oosumposterior = 1. / sumposterior;
				for (int i1 = 0; i1 < means[0].length; i1++) {
					means[i][i1] *= oosumposterior;
				}
				System.out.println(Arrays.toString(means[i]));
			}

			// Update Covariances
			for (int i = 0; i < means.length; i++) {
				double[][] covv = covariances[i];
				double[] mean = means[i];

				for (int j = 0; j < means[0].length; j++) {
					for (int j2 = 0; j2 < means[0].length; j2++) {
						covv[j][j2] = 0.;
					}
				}

				double sumposterior = 0.;
				for (int l = 0; l < data.length; l++) {
					double[] d = data[l];
					double pw = Rspb[l][i];
					sumposterior += pw;
					for (int k = 0; k < means[0].length; k++) {
						for (int j = 0; j < means[0].length; j++) {
							covv[k][j] += (pw * (d[k] - mean[k]) * (d[j] - mean[j]));
						}

					}

				}

				double oosumposterior = 1. / sumposterior;
				for (int j = 0; j < means[0].length; j++) {
					for (int j2 = 0; j2 < means[0].length; j2++) {
						covv[j][j2] *= oosumposterior;
					}
				}
				Matrix cova = new Matrix(covv);
				cova = cova.plus(Matrix.identity(covv.length, covv.length).times(0.1));

				covariances[i] = cova.getArray();
				System.out.println(Arrays.deepToString(covariances[i]));

			}

			// Update Pi
			for (int i = 0; i < means.length; i++) {
				double sum = 0;
				for (int j = 0; j < data.length; j++) {
					sum += Rspb[j][i];
				}
				pi[i] = sum / data.length;
			}

			double NewLogLikelihood = 0, ProbOfX = 0;

			for (int i = 0; i < data.length; i++) {
				ProbOfX = 0;
				for (int j = 0; j < numOfCluster; j++) {
					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					System.err.println(i);
					ProbOfX = ProbOfX + pi[j] * mnd.density(data[i]);

				}
				NewLogLikelihood = NewLogLikelihood + log2(ProbOfX);
			}

			ChangeInLogLikelihood = NewLogLikelihood - LogLikelihood;
			TotalLoglikelihood = NewLogLikelihood;
			System.out.println("Total LogLikelihood: " + TotalLoglikelihood);
			System.out.println("Change in LogLikelihood: " + ChangeInLogLikelihood);

		}
	}

	public void EMKmeans1(double [][] data) {
		int numCluster=this.numCLuster;
		this.data=data;
		//initKmeans(data, numCluster,null);
		double TotalLoglikelihood = 0;
		double ChangeInLogLikelihood = 1;
		while (ChangeInLogLikelihood > 0.0001) {
			double[][] Rspb = new double[data.length][numCluster];
			double LogLikelihood = 0, ProbOfXn = 0;
			for (int i = 0; i < data.length; i++) {
				ProbOfXn = 0;
				//MultivariateNormalDistribution mn = new MultivariateNormalDistribution(means[2], covariances[2]);
				/*if (mn.density(data[i]) != 0.) {
					JOptionPane.showMessageDialog(null, "ada nig");
				}*/
				// System.out.println("ngetes nih"+mn.density(data[i]));
				for (int j = 0; j < numCluster; j++) {

					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					//System.out.println(j + "   " + mnd.density(data[i]));
					 System.err.println(i);
					 Rspb[i][j] = pi[j] * mnd.density(data[i]);
					 ProbOfXn = ProbOfXn + pi[j] * mnd.density(data[i]);

				}

				LogLikelihood = LogLikelihood + log2(ProbOfXn);
			}

			// E-step
			

			for (int i = 0; i < data.length; i++) {
				ProbOfXn = 0;
				for (int j = 0; j < numCluster; j++) {
					//MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					
					// System.out.println("rspb"+Rspb[i][j]);
					ProbOfXn = ProbOfXn + Rspb[i][j];
				}

				for (int j = 0; j < numCluster; j++) {
					Rspb[i][j] = Rspb[i][j] / ProbOfXn;
				}

			}

			// M-Step
			double[] N_k = new double[numCluster];

			for (int k = 0; k < numCluster; k++) { // Calculating N_k's
				for (int n = 0; n < data.length; n++) {
					N_k[k] = N_k[k] + Rspb[n][k];
				}
			}

			// Update Means
			for (int i = 0; i < means.length; i++) {
				for (int j = 0; j < means[0].length; j++) {
					means[i][j] = 0;
				}

				double sumposterior = 0;
				for (int i1 = 0; i1 < data.length; i1++) {
					double pw = Rspb[i1][i];
					sumposterior += pw;
					// System.out.println("pw"+pw);
					for (int j = 0; j < means[0].length; j++) {
						means[i][j] += (pw * data[i1][j]);
						// System.out.println(means[i][j]);
					}
				}

				double oosumposterior = 1. / sumposterior;
				for (int i1 = 0; i1 < means[0].length; i1++) {
					means[i][i1] *= oosumposterior;
				}

				System.out.println(Arrays.toString(means[i]));
			}

			// Update Covariances
			for (int i = 0; i < means.length; i++) {
				double[][] covv = covariances[i];
				double[] mean = means[i];

				for (int j = 0; j < means[0].length; j++) {
					for (int j2 = 0; j2 < means[0].length; j2++) {
						covv[j][j2] = 0.;
					}
				}

				double sumposterior = 0.;
				for (int l = 0; l < data.length; l++) {
					double[] d = data[l];
					double pw = Rspb[l][i];
					sumposterior += pw;
					for (int k = 0; k < means[0].length; k++) {
						for (int j = 0; j < means[0].length; j++) {
							covv[k][j] += (pw * (d[k] - mean[k]) * (d[j] - mean[j]));
						}

					}

				}

				double oosumposterior = 1. / sumposterior;
				for (int j = 0; j < means[0].length; j++) {
					for (int j2 = 0; j2 < means[0].length; j2++) {
						covv[j][j2] *= oosumposterior;
					}
				}

				Matrix cova = new Matrix(covv);
				cova = cova.plus(Matrix.identity(covv.length, covv.length).times(0.1));

				covariances[i] = cova.getArray();
				System.out.println(Arrays.deepToString(covariances[i]));

				// covariances[i] = covv;

			}

			// Update Pi
			for (int i = 0; i < means.length; i++) {
				double sum = 0;
				for (int j = 0; j < data.length; j++) {
					sum += Rspb[j][i];
				}
				pi[i] = sum / data.length;
			}

			double NewLogLikelihood = 0, ProbOfX = 0;

			for (int i = 0; i < data.length; i++) {
				ProbOfX = 0;
				for (int j = 0; j < numCluster; j++) {
					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					System.err.println(i);
					ProbOfX = ProbOfX + pi[j] * mnd.density(data[i]);

				}
				NewLogLikelihood = NewLogLikelihood + log2(ProbOfX);
			}

			ChangeInLogLikelihood = NewLogLikelihood - LogLikelihood;
			TotalLoglikelihood = NewLogLikelihood;
			System.out.println("Total LogLikelihood: " + TotalLoglikelihood);
			System.out.println("Change in LogLikelihood: " + ChangeInLogLikelihood);

		}
	}

	public void EMKmeans(double[][] data) {

		double TotalLoglikelihood = 0;
		double ChangeInLogLikelihood = 1;
		while (ChangeInLogLikelihood > 0.0001) {

			double LogLikelihood = 0, ProbOfXn = 0;
			for (int i = 0; i < data.length; i++) {
				ProbOfXn = 0;

				for (int j = 0; j < this.ClusterIndex.size(); j++) {

					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					// System.err.println("sempet");
					System.out.println(j + "   " + mnd.density(data[i]));
					// System.out.println(mnd.density(data[i])+" "+i+" "+j);
					ProbOfXn = ProbOfXn + pi[j] * mnd.density(data[i]);

				}

				LogLikelihood = LogLikelihood + log2(ProbOfXn);
			}

			// E-step
			double[][] Rspb = new double[data.length][this.ClusterIndex.size()];

			for (int i = 0; i < data.length; i++) {
				ProbOfXn = 0;
				for (int j = 0; j < this.ClusterIndex.size(); j++) {
					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
					Rspb[i][j] = pi[j] * mnd.density(data[i]);
					ProbOfXn = ProbOfXn + Rspb[i][j];
				}

				for (int j = 0; j < this.ClusterIndex.size(); j++) {
					Rspb[i][j] = Rspb[i][j] / ProbOfXn;
				}

			}

			// M-Step
			double[] N_k = new double[this.ClusterIndex.size()];

			for (int k = 0; k < this.ClusterIndex.size(); k++) { // Calculating N_k's
				for (int n = 0; n < data.length; n++) {
					N_k[k] = N_k[k] + Rspb[n][k];
				}
			}

			// Update Means
			for (int i = 0; i < means.length; i++) {
				for (int j = 0; j < means[0].length; j++) {
					means[i][j] = 0;
				}

				double sumposterior = 0;
				for (int i1 = 0; i1 < data.length; i1++) {
					double pw = Rspb[i1][i];
					sumposterior += pw;
					for (int j = 0; j < means[0].length; j++) {
						means[i][j] += (pw * data[i1][j]);
					}
				}

				double oosumposterior = 1. / sumposterior;
				for (int i1 = 0; i1 < means[0].length; i1++) {
					means[i][i1] *= oosumposterior;
				}
			}

			// Update Covariances
			for (int i = 0; i < means.length; i++) {
				double[][] covv = covariances[i];
				double[] mean = means[i];

				for (int j = 0; j < means[0].length; j++) {
					for (int j2 = 0; j2 < means[0].length; j2++) {
						covv[j][j2] = 0.;
					}
				}

				double sumposterior = 0.;
				for (int l = 0; l < data.length; l++) {
					double[] d = data[l];
					double pw = Rspb[l][i];
					sumposterior += pw;
					for (int k = 0; k < means[0].length; k++) {
						for (int j = 0; j < means[0].length; j++) {
							covv[k][j] += (pw * (d[k] - mean[k]) * (d[j] - mean[j]));
						}

					}

				}

				double oosumposterior = 1. / sumposterior;
				for (int j = 0; j < means[0].length; j++) {
					for (int j2 = 0; j2 < means[0].length; j2++) {
						covv[j][j2] *= oosumposterior;
					}
				}

				Matrix cova = new Matrix(covv);
				cova = cova.plus(Matrix.identity(covv.length, covv.length).times(0.1));

				covariances[i] = cova.getArray();
				System.out.println(Arrays.deepToString(covariances[i]));

				// covariances[i] = covv;

			}

			// Update Pi
			for (int i = 0; i < means.length; i++) {
				double sum = 0;
				for (int j = 0; j < data.length; j++) {
					sum += Rspb[j][i];
				}
				pi[i] = sum / data.length;
			}

			double NewLogLikelihood = 0, ProbOfX = 0;

			for (int i = 0; i < data.length; i++) {
				ProbOfX = 0;
				for (int j = 0; j < this.ClusterIndex.size(); j++) {
					MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);

					ProbOfX = ProbOfX + pi[j] * mnd.density(data[i]);

				}
				NewLogLikelihood = NewLogLikelihood + log2(ProbOfX);
			}

			ChangeInLogLikelihood = NewLogLikelihood - LogLikelihood;
			TotalLoglikelihood = NewLogLikelihood;
			System.out.println("Total LogLikelihood: " + TotalLoglikelihood);
			System.out.println("Change in LogLikelihood: " + ChangeInLogLikelihood);

		}
	}

	public void EMApache(double[][] data) {
		MixtureMultivariateNormalDistribution mmnd = new MixtureMultivariateNormalDistribution(pi, means, covariances);
		MultivariateNormalMixtureExpectationMaximization mnmem = new MultivariateNormalMixtureExpectationMaximization(
				data);
		mnmem.fit(mmnd);
		mmnd = mnmem.getFittedModel();
		loglikelihoods = mnmem.getLogLikelihood();

		for (int i = 0; i < mmnd.getComponents().size(); i++) {

			covariances[i] = mmnd.getComponents().get(i).getValue().getCovariances().getData();
			means[i] = mmnd.getComponents().get(i).getValue().getMeans();
			pi[i] = mmnd.getComponents().get(i).getFirst();
		}
		
		System.out.println("log= "+loglikelihoods);

	}

	public int predict(double[] data,int length) {
		int index = 0;
		// int check = Integer.parseInt(JOptionPane.showInputDialog("Masukin"));
		/*RealVector rv = MatrixUtils.createRealVector(data);
		Matrix m = new Matrix(data, 1);
		Jama.SingularValueDecomposition svd = new Jama.SingularValueDecomposition(m);
		double [][] s = svd.getS().getArray();
		double [][]vt = svd.getV().transpose().getArray();
		double[] unew = null,vnew = null;
		for (int i = 0; i < 1; i++) {
		unew = s[i];
			vnew = vt[i];
		}

		Matrix mnew = svd.getU().times(new Matrix(unew, 1)).times(new Matrix(vnew,1));
		data=mnew.getArray()[0];*/
		double[] shortdata = new double[length];
		
			for (int i = 0; i < length; i++) {
			if (i<data.length) {
				shortdata[i] = data[i];
			}
			else {
				shortdata[i]=0.;
			}
			
		}
			ArrayList<Double> features = new ArrayList<>();
			for (int i = 0; i < shortdata.length; i++) {
				features.add(shortdata[i]);
			}
			
			featuretest.add(features);
			
			/*XYSeries xys = new XYSeries("Sampel 4");
			for (int i = 0; i < shortdata.length; i++) {
				xys.add(i, shortdata[i]);
			}
			XYSeriesCollection dataset = new XYSeriesCollection();
			dataset.addSeries(xys);
			XyChart.dataSet=dataset;
			XyChart chart = new XyChart("Sampel 4", "Sampel 4");
			chart.pack();
			RefineryUtilities.centerFrameOnScreen(chart);
			chart.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			chart.setVisible(true);*/
			
		double ProbOfXns = 0;
		double[] prob = new double[means.length];
		for (int j = 0; j < means.length; j++) {
			ProbOfXns = 0;
			MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means[j], covariances[j]);
			ProbOfXns = pi[j] * mnd.density(shortdata);
			prob[j] = ProbOfXns;
		}

			double set = prob[0];
			index = 0;
			for (int i1 = 1; i1 < prob.length; i1++) {

				if (prob[i1] > set) {
					set = prob[i1];
					index = i1;
				}

			}

		
		// System.out.println(index);
		return (index);
	}

	public void printGMMResult() {

		System.out.println("covariances baru = " + Arrays.deepToString(covariances));
		System.out.println("Means baru = " + Arrays.deepToString(means));
		System.out.println("Phi baru = " + Arrays.toString(pi));

	}

}
