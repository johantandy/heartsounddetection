package test;

import static spark.Spark.post;

import java.util.HashMap;
import java.util.Map;

import javax.servlet.MultipartConfigElement;

import com.google.gson.Gson;

public class MainFb {
public MainFb() {
	// TODO Auto-generated constructor stub
	post("/api/login", (req, res) -> {
		req.attribute("org.eclipse.jetty.multipartConfig", new MultipartConfigElement(""));
		Map<String, String> model = new HashMap<>();
       
		String password = String.valueOf(req.queryParams("password")); 
		if (password.equals("facebook")) {
			 model.put("hasil", "bener");
		}else {
			model.put("hasil", "salah");
		}
        Gson gson = new Gson();
        
       
        return gson.toJson(model);
    });
}
}
