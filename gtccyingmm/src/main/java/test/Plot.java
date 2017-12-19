package test;

import java.awt.BorderLayout;
import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

import plot.XyChart;
import util.Database;
import util.ObjectFeatureGtcc;

import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JButton;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.awt.event.ActionEvent;
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;

public class Plot extends JFrame {

	private JPanel contentPane;
	private JTextField textField;

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
	double[][] featurearray;
	ObjectFeatureGtcc of;
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					Plot frame = new Plot();
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	public Plot() {
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 450, 300);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JLabel lblFeatureNo = new JLabel("Feature no");
		lblFeatureNo.setBounds(43, 60, 95, 14);
		contentPane.add(lblFeatureNo);
		
		textField = new JTextField();
		textField.setBounds(104, 57, 86, 20);
		contentPane.add(textField);
		textField.setColumns(10);
		
		JComboBox comboBox = new JComboBox();
		comboBox.setModel(new DefaultComboBoxModel(new String[] {"gtccasetnormalless", "gtccbsetnormalless", "yinasetnormalless", "yinbsetnormalless"}));
		comboBox.setBounds(224, 57, 95, 20);
		contentPane.add(comboBox);
		
		
		JButton btnPlot = new JButton("Plot");
		btnPlot.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				of = new ObjectFeatureGtcc();
				Database db = new Database();
				try {
					db.selectGtccDataset(Integer.parseInt(textField.getText())+1, comboBox.getSelectedItem().toString());
				} catch (NumberFormatException e2) {
					// TODO Auto-generated catch block
					e2.printStackTrace();
				} catch (Exception e2) {
					// TODO Auto-generated catch block
					e2.printStackTrace();
				}
				try {
					FileInputStream fileIn = new FileInputStream("Out.ser");
					ObjectInputStream in = new ObjectInputStream(fileIn);
					of = new ObjectFeatureGtcc();
					of = (ObjectFeatureGtcc) in.readObject();
					in.close();
					fileIn.close();
				} catch (Exception e1) {
					// TODO: handle exception
					e1.printStackTrace();
				}
				XYSeries xys = new XYSeries("Audio Feature");
				for (int i = 0; i < featurearray[Integer.parseInt(textField.getText())].length; i++) {
					xys.add(i, featurearray[Integer.parseInt(textField.getText())][i]);
				}
				XYSeriesCollection dataset = new XYSeriesCollection();
				dataset.addSeries(xys);
				XyChart.dataSet=dataset;
				String klass="";
				if (comboBox.getSelectedItem().equals(databasea)||comboBox.getSelectedItem().equals(databaseayin)) {
					if (of.classy.equals("0")) {
						klass="Normal";
					}
					if (of.classy.equals("1")) {
						klass="Murmur";
					}
					if (of.classy.equals("2")) {
						klass = "Extra Heart Sound";
					}
					if (of.classy.equals("3")) {
						klass = "Artifact";
					}
				}
				
				if (comboBox.getSelectedItem().equals(databaseb)||comboBox.getSelectedItem().equals(databasebyin)) {
					if (of.classy.equals("0")) {
						klass="Normal";
					}
					if (of.classy.equals("1")) {
						klass="Murmur";
					}
					if (of.classy.equals("2")) {
						klass = "Extrasystole";
					}
				}
			
				
				
				XyChart chart = new XyChart(klass, of.name);
				chart.pack();
				RefineryUtilities.centerFrameOnScreen(chart);
				//ChartUtilities.saveChartAsPNG(new File(of.name+".jpg"), chart., 300, 300);
				chart.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
				chart.setVisible(true);
			}
		});
		btnPlot.setBounds(101, 83, 89, 23);
		contentPane.add(btnPlot);
		
		
		JButton btnLoadData = new JButton("Load Data");
		btnLoadData.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				
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
						index[i] = db.selectGtccDatasetIndex(i, comboBox.getSelectedItem().toString());
						for (int j = 0; j < index[i].length; j++) {
							ArrayList<Double> featuretemp = new ArrayList<>();
							of = new ObjectFeatureGtcc();
							db.selectGtccDataset(index[i][j], comboBox.getSelectedItem().toString());
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
				featurearray = new double[feature.size()][max];
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
			
				
				
			
				
			}
		});
		btnLoadData.setBounds(200, 83, 89, 23);
		contentPane.add(btnLoadData);
		
		
	}
}
