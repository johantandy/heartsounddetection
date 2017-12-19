package gtccyingmm;
import static spark.Spark.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.sql.Date;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

import javax.servlet.MultipartConfigElement;
import javax.servlet.http.Part;
import javax.servlet.Filter;

import org.apache.velocity.VelocityContext;
import org.eclipse.jetty.http.MetaData.Request;
import org.jfree.chart.plot.Plot;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import spark.ModelAndView;
import spark.template.velocity.VelocityTemplateEngine;
import spark.utils.IOUtils;
import util.ObjectDB;

public class Main {

	public static String status;
	
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
	
	public static void main(String[] args) {
		
	/*	File uploadDir = new File("upload");
        uploadDir.mkdir(); // create the upload directory if it doesn't exist*/

       readProp();
       // staticFileLocation("/public");
        externalStaticFileLocation("public");
        port(Integer.parseInt(ObjectDB.webport));
        //staticFileLocation(publik.toURI().getPath());
		
		// TODO Auto-generated method stub
        	
        get("/plot", (request, response) -> {
        	test.Plot plot = new test.Plot();
        	plot.main(null);
		      return new ModelAndView(new HashMap(), "templates/index.vtl");
		      //return new ModelAndView(new HashMap(), indexx.getAbsolutePath());
		    }, new VelocityTemplateEngine());
			
			  get("/", (request, response) -> {
			      return new ModelAndView(new HashMap(), "templates/index.vtl");
			      //return new ModelAndView(new HashMap(), indexx.getAbsolutePath());
			    }, new VelocityTemplateEngine());
			  
			  get("/homelogin", (request, response) -> {
			      return new ModelAndView(new HashMap(), "templates/homelogin.vtl");
			    }, new VelocityTemplateEngine());
			  get("/api/test", (request, response) -> {
			      return "helloworld";
			    });
			  
			  get("/about", (request, response) -> {
			      return new ModelAndView(new HashMap(), "templates/about.vtl");
			      //return new ModelAndView(new HashMap(), admin.getAbsolutePath());
			    }, new VelocityTemplateEngine());
			  	
			  
			  
			  get("/admin", (request, response) -> {
			      return new ModelAndView(new HashMap(), "templates/admin.vtl");
			      //return new ModelAndView(new HashMap(), admin.getAbsolutePath());
			    }, new VelocityTemplateEngine());
			  	
		
			
			//UPLOAD FILE DATASET A API
			post("/api/uploada", (req, res) -> {
				System.out.println("Upload A Dataset"+LocalDateTime.now());
				req.attribute("org.eclipse.jetty.multipartConfig", new MultipartConfigElement("upload\\0"));
	            String directory = "";
	            Part filePart=null;
	            directory = String.valueOf(req.queryParams("selectGtccYinA")); 
				filePart = req.raw().getPart("fileGtccYinA");
					String path = "upload\\0\\"+directory+"\\"+ filePart.getSubmittedFileName();
	            
	            try (InputStream inputStream = filePart.getInputStream()) {
	            		
	                OutputStream outputStream = new FileOutputStream(path);
	                IOUtils.copy(inputStream, outputStream);
	                outputStream.close();
	            }
	            
	            Controller con = new Controller();
	            String gtccreply = con.SaveFeatureDatasetA(path, directory, filePart.getSubmittedFileName());
	            String yinreply = con.SaveFeatureDatasetAYin(path, directory, filePart.getSubmittedFileName());
	            Map<String, String> model = new HashMap<>();
	            model.put("gtcc", gtccreply);
	            model.put("yin", yinreply);
	            Gson gson = new Gson();
	            return gson.toJson(model);
	        });
			
			post("/api/uploadtest", (req, res) -> {
				System.out.println("Upload A Test File "+LocalDateTime.now());
				req.attribute("org.eclipse.jetty.multipartConfig", new MultipartConfigElement("testing"));
	            Part filePart=null;
					
					filePart = req.raw().getPart("file");
					String path = "testing\\"+ filePart.getSubmittedFileName();
	            
	            try (InputStream inputStream = filePart.getInputStream()) {
	            	
	                OutputStream outputStream = new FileOutputStream(path);
	                IOUtils.copy(inputStream, outputStream);
	                outputStream.close();
	            }
	            
	            System.out.println("\"Upload\"");
	            System.out.println(path);
	            if (path.contains(".amr")) {
	            	File file = new File(path);
	            	File ffmpeg = new File("ffmpeg\\bin");
	            	ProcessBuilder builder = new ProcessBuilder(
	                        "cmd.exe", "/c", "ffmpeg.exe -i "+"\""+file.getAbsolutePath()+"\" "+"\""+file.getAbsolutePath().replace(".amr", ".wav")+"\"");
	                    builder.redirectErrorStream(true);
	                    Process p = builder.start();
	                    System.out.println("ffmpeg.exe -i "+"\""+file.getAbsolutePath()+"\" "+"\""+file.getAbsolutePath().replace(".amr", ".wav")+"\"");
	                    BufferedReader r = new BufferedReader(new InputStreamReader(p.getInputStream()));
	                    String line;
	                    while (true) {
	                        line = r.readLine();
	                        if (line == null) { break; }
	                        System.out.println(line);
	                    }
	                    path=file.getAbsolutePath().replace(".amr", ".wav");
	            
	            }
	            
	            if (path.contains(".m4a")) {
	            	File file = new File(path);
	            	File ffmpeg = new File("ffmpeg\\bin");
	            	ProcessBuilder builder = new ProcessBuilder(
	                        "cmd.exe", "/c", "ffmpeg.exe -i "+"\""+file.getAbsolutePath()+"\" "+"\""+file.getAbsolutePath().replace(".m4a", ".wav")+"\"");
	                    builder.redirectErrorStream(true);
	                    Process p = builder.start();
	                    System.out.println("ffmpeg.exe -i "+"\""+file.getAbsolutePath()+"\" "+"\""+file.getAbsolutePath().replace(".m4a", ".wav")+"\"");
	                    BufferedReader r = new BufferedReader(new InputStreamReader(p.getInputStream()));
	                    String line;
	                    while (true) {
	                        line = r.readLine();
	                        if (line == null) { break; }
	                        System.out.println(line);
	                    }
	                    path=file.getAbsolutePath().replace(".m4a", ".wav");
	            
	            }
	            
	            Controller con = new Controller();
	            
	            int gtcca = con.PredictAset(path);
	            int gtccb = con.PredictBset(path);
	            int yina = con.PredictAsetYin(path);
	            int yinb = con.PredictBsetYin(path);
	            String Agtcc = null,Bgtcc = null,Ayin = null,Byin = null;
	            if (gtcca==0) {
					Agtcc="Normal";
				}
	            if (gtcca==1) {
					Agtcc="Murmur";
				}
	            if (gtcca==2) {
					Agtcc="Extra Heart Sound";
				}
	            if (gtcca==3) {
					Agtcc="Artifact";
				}
	            
	            if (yina==0) {
					Ayin="Normal";
				}
	            if (yina==1) {
					Ayin="Murmur";
				}
	            if (yina==2) {
					Ayin="Extra Heart Sound";
				}
	            if (yina==3) {
					Ayin="Artifact";
				}
	            
	            if (gtccb==0) {
					Bgtcc="Normal";
				}
	            if (gtccb==1) {
					Bgtcc="Murmur";
				}
	            if (gtccb==2) {
					Bgtcc="Extrasystole";
				}
	            
	            if (yinb==0) {
					Byin="Normal";
				}
	            if (yinb==1) {
					Byin="Murmur";
				}
	            if (yinb==2) {
					Byin="Extrasystole";
				}
	            
	            
	            
	            Map<String, String> model = new HashMap<>();
	            model.put("gtcca", Agtcc);
	            model.put("gtccb", Bgtcc);
	            model.put("yina", Ayin);
	            model.put("yinb", Byin);
	            Gson gson = new Gson();
	           return gson.toJson(model) ;
	        });
			
			//UPLOAD FILE DATASET B API
			post("/api/uploadb", (req, res) -> {
				System.out.println("Upload B Dataset "+LocalDateTime.now());
				req.attribute("org.eclipse.jetty.multipartConfig", new MultipartConfigElement("upload\\1"));
	            String directory = "";
	            Part filePart=null;
					directory = String.valueOf(req.queryParams("selectGtccYinB")); 
					filePart = req.raw().getPart("fileGtccYinB");
					String path = "upload\\1\\"+directory+"\\"+ filePart.getSubmittedFileName();
	            
	            try (InputStream inputStream = filePart.getInputStream()) {
	                OutputStream outputStream = new FileOutputStream(path);
	                IOUtils.copy(inputStream, outputStream);
	                outputStream.close();
	            }
	            
	            Controller con = new Controller();
	            String gtccreply = con.SaveFeatureDatasetB(path, directory, filePart.getSubmittedFileName());
	            String yinreply = con.SaveFeatureDatasetBYin(path, directory, filePart.getSubmittedFileName());
	            Map<String, String> model = new HashMap<>();
	            model.put("gtcc", gtccreply);
	            model.put("yin", yinreply);
	            Gson gson = new Gson();
	            return gson.toJson(model);

	            
	        });
			
			//TRAIN API DATASET A & B
			post("/api/train", (req, res) -> {
				int choose = Integer.valueOf(req.queryParams("selectdataset"));
				Controller con = new Controller();
				if (choose==0) {
					return "Berhasil Training dengan loglikelihood = "+con.TrainAsetWork();
				}
				else {
					return "Berhasil Training dengan loglikelihood = "+con.TrainBsetWork();
				}
				
				
				});
			
			
			get("/api/trainGtccA", (request, response) -> {
				System.out.println("Train GTCC A");
				status="Training process called";
				Controller con = new Controller();
				Gson gson = new Gson();
			      return gson.toJson(con.TrainAsetWeb());
			    });
			
			get("/api/trainGtccB", (request, response) -> {
				System.out.println("Train GTCC B");
				Controller con = new Controller();
				Gson gson = new Gson();
			      return gson.toJson(con.TrainBsetWeb());
			    });
			
			get("/api/trainYinA", (request, response) -> {
				System.out.println("Train YIN A");
				Controller con = new Controller();
				Gson gson = new Gson();
			      return gson.toJson(con.TrainAsetWebYIN());
			    });
			
			get("/api/trainYinB", (request, response) -> {
				System.out.println("Train YIN B");
				Controller con = new Controller();
				Gson gson = new Gson();
			      return gson.toJson(con.TrainBsetWebYIN());
			    });
			
			
			post("/api/predict", (req, res) -> {
				int choose = Integer.valueOf(req.queryParams("selectdataset"));
				
				req.attribute("org.eclipse.jetty.multipartConfig", new MultipartConfigElement("test"));
	            String directory = "";
	            Part filePart=null;
					filePart = req.raw().getPart("myfiletest");
			String path = "test\\"+ filePart.getSubmittedFileName();
	            
	            try (InputStream inputStream = filePart.getInputStream()) {
	                OutputStream outputStream = new FileOutputStream(path);
	                IOUtils.copy(inputStream, outputStream);
	                outputStream.close();
	            }
				
				Controller con = new Controller();
				if (choose==0) {
					return "hasil prediksi = "+con.PredictAset(path);
				}
				else {
					return "hasil prediksi = "+con.PredictAset(path);
				}
				
				
				});
			




	}

}
