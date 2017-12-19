package util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.sql.Connection;
import java.sql.Date;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;

import com.mysql.jdbc.Statement;

public class Database {

	public void insertGtccDataset(String classy, double[] feature, String name, String database) {

		try {
			ObjectFeatureGtcc ofg = new ObjectFeatureGtcc();
			ofg.name = name;
			ofg.data = feature;
			ofg.classy = classy;
			FileOutputStream fileout = new FileOutputStream("data.ser");
			ObjectOutputStream out = new ObjectOutputStream(fileout);
			out.writeObject(ofg);
			out.close();
			fileout.close();
			Connection conn = DB.getConnection();
			DB.insertobject("data.ser", conn, classy, name, database);
			conn.close();
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
	}

	public void insertYinDataset(String classy, double[] feature, String name, String database) {

		try {
			ObjectFeatureGtcc ofg = new ObjectFeatureGtcc();
			ofg.name = name;
			ofg.data = feature;
			ofg.classy = classy;
			FileOutputStream fileout = new FileOutputStream("data.ser");
			ObjectOutputStream out = new ObjectOutputStream(fileout);
			out.writeObject(ofg);
			out.close();
			fileout.close();
			Connection conn = DB.getConnection();
			DB.insertobject("data.ser", conn, classy, name, database);
			conn.close();
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
	}

	public void selectGtccDataset(int index, String database) throws Exception {
		Connection conn = DB.getConnection();
		String sql = "SELECT data from " + database + " where id=?";
		PreparedStatement pstmt = conn.prepareStatement(sql);
		pstmt.setInt(1, index);
		ResultSet rs = pstmt.executeQuery();

		File file = new File("Out.ser");
		FileOutputStream output = new FileOutputStream(file);

		while (rs.next()) {
			InputStream input = rs.getBinaryStream("data");
			byte[] buffer = new byte[100];
			while (input.read(buffer) > 0) {
				output.write(buffer);
			}
		}
		conn.close();

	}
	
	public void deleteData(int index, String database) throws Exception {
		Connection conn = DB.getConnection();
		String sql = "delete from " + database + " where id=?";
		PreparedStatement pstmt = conn.prepareStatement(sql);
		pstmt.setInt(1, index);
		pstmt.executeUpdate();
		pstmt.close();
	}
	
	public int getCount(String database) throws Exception {
		int count=0;
		Connection conn = DB.getConnection();
		String sql = "SELECT count(*) as count from " + database;
		PreparedStatement pstmt = conn.prepareStatement(sql);
		ResultSet rs = pstmt.executeQuery();
		while(rs.next()) {
			count = rs.getInt("count");
		}
		pstmt.close();
		return count;
	}


	public int[] selectGtccDatasetIndex(int classy, String database) throws Exception {
		Connection conn = DB.getConnection();
		String sql = "SELECT id from " + database + " where class=?";
		PreparedStatement pstmt = conn.prepareStatement(sql);
		pstmt.setInt(1, classy);
		ResultSet rs = pstmt.executeQuery();
		ArrayList<Integer> index = new ArrayList<>();
		while (rs.next()) {
			index.add(rs.getInt("id"));
		}
		int[] indexarray = new int[index.size()];

		for (int i = 0; i < index.size(); i++) {
			indexarray[i] = index.get(i);
		}
		conn.close();
		return indexarray;
		

	}

	public void insertGMMParam(double[][][] covariances, double[][] means, double[] pi, double accuracy,
			String database) {
		try {
			DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
			java.util.Date date = new java.util.Date();
			ObjectParamGMM opg = new ObjectParamGMM();
			opg.covariances = covariances;
			opg.means = means;
			opg.pi = pi;
			opg.inserton=dateFormat.format(date);
			opg.accuracy=accuracy;
			FileOutputStream fileout = new FileOutputStream("data.ser");
			ObjectOutputStream out = new ObjectOutputStream(fileout);
			out.writeObject(opg);
			out.close();
			fileout.close();
			Connection conn = DB.getConnection();
			DB.insertParamGmm("data.ser",opg.inserton,accuracy,conn, database);
			conn.close();
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
	}
	
	public void insertTerpaksa() {
		try {
			DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
			java.util.Date date = new java.util.Date();
			Connection conn = DB.getConnection();
			DB.insertParamGmm("data.ser",dateFormat.format(date),100,conn, "gtccgmmbsetnormalless");
			conn.close();
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
	}
	
	public int selectlength(String database) throws Exception {
		int length = 0;
		Connection conn = DB.getConnection();
		String sql = "SELECT length from " + database + " order by id desc limit 1";
		PreparedStatement pstmt = conn.prepareStatement(sql);
		ResultSet rs = pstmt.executeQuery();

		while (rs.next()) {
			length=rs.getInt("length");
			}
		return length;
		
	}

	public void selectGtccParam(String database) throws Exception {
		Connection conn = DB.getConnection();
		String sql = "SELECT param from " + database + " order by id desc limit 1";
		PreparedStatement pstmt = conn.prepareStatement(sql);
		ResultSet rs = pstmt.executeQuery();
		File file = new File("Out.ser");
		FileOutputStream output = new FileOutputStream(file);

		while (rs.next()) {
			InputStream input = rs.getBinaryStream("param");
			byte[] buffer = new byte[1024];
			while (input.read(buffer) > 0) {
				output.write(buffer);
			}
		}
		conn.close();

	}
	
	public void truncateDb(String database) throws Exception {
		Connection conn = DB.getConnection();
		String sql = "TRUNCATE TABLE " + database;
		PreparedStatement pstmt = conn.prepareStatement(sql);
		pstmt.executeQuery();
	}

}
