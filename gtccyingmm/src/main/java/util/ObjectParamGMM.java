package util;

import java.sql.Date;

public class ObjectParamGMM implements java.io.Serializable {
	public double[][] means;
	public double[][][] covariances;
	public double[] pi;
	public String inserton;
	public double accuracy;
}
