package util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import com.mysql.jdbc.Statement;

public class DB {
	
	static String sql = "INSERT INTO objectfeature (Object,String) VALUES(?,?)";
	
	public static Connection getConnection() throws Exception {
	    String driver = "org.gjt.mm.mysql.Driver";
	    String url = "jdbc:mysql://"+ObjectDB.url+"/GMMDB";
	    String username = ObjectDB.username;
	    String password = ObjectDB.password;
	    Class.forName(driver);
	    Connection conn = DriverManager.getConnection(url, username, password);
	    return conn;
	  }
	
	public static void writebyte(Connection conn, byte[] buffer,String name) throws Exception {
		PreparedStatement pstmt = conn.prepareStatement(sql);
	    pstmt.setBytes(1, buffer);
	    pstmt.setString(2, name);
	    pstmt.executeUpdate();
	    pstmt.close();
	  }
	
	public static void insertobject(String path ,Connection conn,String classy,String name,String Database) throws FileNotFoundException, SQLException{
		InputStream inputStream = new FileInputStream(new File(path));
		 
		String sql = "INSERT INTO "+Database+" (name,data,class) values (?,?,?)";
		PreparedStatement statement = conn.prepareStatement(sql);
		statement.setString(1, name);
		statement.setBlob(2, inputStream);
		statement.setString(3, classy);
		
		statement.executeUpdate();
		conn.close();
	}
	
	public static void insertParamGmm(String path,String date ,double accuracy,Connection conn,String Database) throws FileNotFoundException, SQLException{
		InputStream inputStream = new FileInputStream(new File(path));
		 
		String sql = "INSERT INTO "+Database+" (param,inserton,accuracy) values (?,?,?)";
		PreparedStatement statement = conn.prepareStatement(sql);
		statement.setBlob(1, inputStream);
		statement.setString(2, date);
		statement.setDouble(3, accuracy);
		statement.executeUpdate();
		conn.close();
	}
	
	public static void readobject(Connection conn,int i,String Database) throws SQLException, IOException{
		String sql = "SELECT data from "+Database+" where id=?";
		PreparedStatement pstmt = conn.prepareStatement(sql);
		pstmt.setInt(1, i);
		ResultSet rs=pstmt.executeQuery();
		
		File file = new File("Out.ser");
		FileOutputStream output = new FileOutputStream(file);
		
		while (rs.next()) {
			InputStream input = rs.getBinaryStream("data");
            byte[] buffer = new byte[1024];
            while (input.read(buffer) > 0) {
                output.write(buffer);
            }
		}
		conn.close();
		
	}
	
	public static void readallobject(Connection conn,int i,String Database) throws SQLException, IOException{
		String sql = "SELECT data from "+Database;
		java.sql.Statement stmt = conn.createStatement();
		ResultSet rs=stmt.executeQuery(sql);
		
		File file = new File("Out.ser");
		FileOutputStream output = new FileOutputStream(file);
		
		while (rs.next()) {
			InputStream input = rs.getBinaryStream("data");
            byte[] buffer = new byte[1024];
            while (input.read(buffer) > 0) {
                output.write(buffer);
            }
		}
		
	}
	
	public static int getcount(Connection con,String Database) throws SQLException{
		String sql = "select count(*) from "+Database;
		int count=0;
		java.sql.Statement stmt=con.createStatement();
		ResultSet rs = stmt.executeQuery(sql);
		while(rs.next()){
			count = rs.getInt(1);
		}
		con.close();
		rs.close();
		return count;
	}
	
	public static String getkata(Connection con,int i,String Database) throws SQLException{
		String sql = "select kata from "+Database+" where id=?";
		PreparedStatement pstmt = con.prepareStatement(sql);
		pstmt.setInt(1, i);
		ResultSet rs=pstmt.executeQuery();
		String kata=null;
		while(rs.next()){
			kata = rs.getString(1);
		}
		con.close();
		rs.close();
		return kata;
	}
	
	public static String getspeaker(Connection con,int i,String Database) throws SQLException{
		String sql = "select speaker from "+Database+" where id=?";
		PreparedStatement pstmt = con.prepareStatement(sql);
		pstmt.setInt(1, i);
		ResultSet rs=pstmt.executeQuery();
		String speaker=null;
		while(rs.next()){
			speaker = rs.getString(1);
		}
		con.close();
		rs.close();
		return speaker;
	}
	
	public static void insertobjectyin(String path ,Connection conn,String name,String speaker,String Database) throws FileNotFoundException, SQLException{
		InputStream inputStream = new FileInputStream(new File(path));
		 
		String sql = "INSERT INTO "+Database+" (data,kata,speaker) values (?,?,?)";
		PreparedStatement statement = conn.prepareStatement(sql);
		statement.setBlob(1, inputStream);
		statement.setString(2, name);
		statement.setString(3, speaker);
		statement.executeUpdate();
	}
	
	public static void readobjectyin(Connection conn,int i,String Database) throws SQLException, IOException{
		String sql = "SELECT data from "+Database+" where id=?";
		PreparedStatement pstmt = conn.prepareStatement(sql);
		pstmt.setInt(1, i);
		ResultSet rs=pstmt.executeQuery();
		
		File file = new File("Out.ser");
		FileOutputStream output = new FileOutputStream(file);
		
		while (rs.next()) {
			InputStream input = rs.getBinaryStream("data");
            byte[] buffer = new byte[1024];
            while (input.read(buffer) > 0) {
                output.write(buffer);
            }
		}
		conn.close();
		
	}
	
	public static int getcountyin(Connection con,String Database) throws SQLException{
		String sql = "select count(*) from "+Database;
		int count=0;
		java.sql.Statement stmt=con.createStatement();
		ResultSet rs = stmt.executeQuery(sql);
		while(rs.next()){
			count = rs.getInt(1);
		}
		con.close();
		rs.close();
		return count;
	}
	
	public static String getkatayin(Connection con,int i,String Database) throws SQLException{
		String sql = "select kata from "+Database+" where id=?";
		PreparedStatement pstmt = con.prepareStatement(sql);
		pstmt.setInt(1, i);
		ResultSet rs=pstmt.executeQuery();
		String kata=null;
		while(rs.next()){
			kata = rs.getString(1);
		}
		con.close();
		rs.close();
		return kata;
	}
	
	public static String getspeakeryin(Connection con,int i,String Database) throws SQLException{
		String sql = "select speaker from "+Database+" where id=?";
		PreparedStatement pstmt = con.prepareStatement(sql);
		pstmt.setInt(1, i);
		ResultSet rs=pstmt.executeQuery();
		String speaker=null;
		while(rs.next()){
			speaker = rs.getString(1);
		}
		con.close();
		rs.close();
		return speaker;
	}
	

	
}
