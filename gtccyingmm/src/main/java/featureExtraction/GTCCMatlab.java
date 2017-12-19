package featureExtraction;

import java.io.File;
import java.util.Arrays;
import java.util.concurrent.ExecutionException;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/*import com.mathworks.engine.EngineException;
import com.mathworks.engine.MatlabEngine;
import com.mathworks.engine.MatlabExecutionException;
import com.mathworks.engine.MatlabSyntaxException;*/

import util.Database;

public class GTCCMatlab {
	/* static MatlabEngine eng;
	public GTCCMatlab() throws EngineException, IllegalArgumentException, IllegalStateException, InterruptedException {
		System.out.println("init matlab");
		eng = MatlabEngine.startMatlab();
	}
	
	public static void closeMatlab() throws EngineException {
		eng.close();
	}
	
 public double[] extractFeature(String file) throws IllegalArgumentException, IllegalStateException, InterruptedException, MatlabExecutionException, MatlabSyntaxException, ExecutionException {
	 
		String path = "Matlab";
		File matlabfunction = new File(path);
		eng.eval("cd "+matlabfunction.getAbsolutePath());
		System.out.println("Extract");
  double[][] feature = eng.feval("gtfeatures",file);
  
  double [] features = new double[feature.length*feature[0].length];
  int iter =0;
  System.out.println("Storing");
  for (int i = 0; i < feature.length; i++) {
	for (int j = 0; j < feature.length; j++) {
		features[iter]=feature[i][j];
		iter++;
	}
}
  System.out.println("reduce dimension");
	
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
		multip = 22 * itera;
		itera++;
	}
	System.out.println("itera "+tempitera);
	itera=tempitera;

	SingularValueDecomposition svd = new SingularValueDecomposition(MatrixUtils.createRealMatrix(feature));
	RealMatrix u = svd.getU();
	RealMatrix s = svd.getS();
	RealMatrix vt = svd.getVT();
	double[][] sar = s.getData();
	double[][] vtar = vt.getData();
	double[][] news = new double[22][itera];
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
	double[][] featurereduced = featuresvd.getData();
	double[] featuresvdreduc = new double[featurereduced.length*featurereduced[0].length];
	int iterators = 0;
	for (int i = 0; i < featurereduced.length; i++) {
		for (int j = 0; j < featurereduced[0].length; j++) {
			featuresvdreduc[iterators] = featurereduced[i][j];
			iterators++;
		}
	}
  System.out.println("length feature svd "+featuresvdreduc.length);
  //System.out.println(features.length);
return featuresvdreduc;
 }*/
}
