package ca.mcmaster.magarveylab.prism.web;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Timer;

import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.http.HttpEntity;
import org.apache.http.HttpResponse;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.HttpClientBuilder;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.database.JsonInput;
import ca.mcmaster.magarveylab.prism.util.PrismFileWriter;
import ca.mcmaster.magarveylab.prism.web.html.PrismReport;
import ca.mcmaster.magarveylab.wasp.WebApplicationSubmit;
import ca.mcmaster.magarveylab.wasp.exception.ExceptionHandler;
import ca.mcmaster.magarveylab.wasp.session.BasicSession;
import ca.mcmaster.magarveylab.wasp.session.Session;
import ca.mcmaster.magarveylab.wasp.session.SessionListener;
import ca.mcmaster.magarveylab.wasp.session.SessionManager;
import ca.mcmaster.magarveylab.wasp.util.TimeUtil;

/**
 * Servlet to load saved results of a PRISM search, in JSON form, from the
 * database without form submission.
 * 
 * @author skinnider
 */
@WebServlet("/load")
public class LoadJson extends WebApplicationSubmit {

	private static final long serialVersionUID = -1846321311064350374L;
	private static final String authToken = "db_prod S76RRcN6QzFUaSBdoXPs"; //"dev nyCNmeLPpyyGgcepmxc6";
	private static final String host = "http://magarveylab-ws.mcmaster.ca/v1";
	
	@Override
	protected void doGet(HttpServletRequest request, HttpServletResponse response) 
			throws ServletException, IOException {
		// register session
		String sessionID = createSessionID();
		registerSession(sessionID, request);
		Session session = SessionManager.getSessionManager().getSession(sessionID);
				
		// get JSON 
		String id = request.getParameter("id");
		System.out.println("Getting JSON with ID " + id);
		executeHttpRequest(id, session);

		try {
			// parse in PRISM
			String filepath = session.dir() + id + ".json";
			Prism prism = JsonInput.read(filepath, session);
			session.setWebapp(prism);
			
			System.out.println("webapp: " + session.webapp());
			
			// set up report
			PrismReport report = (PrismReport) session.report();
			Timer timer = new Timer();
			timer.schedule(report, 0, 1_000);

			// write reports 
			PrismFileWriter.writeAllFiles(prism.genome(), session);
			prism.terminate();			
		} catch (ParseException e) {
			throw new IOException("Could not read JSON file " + id + "!");
		}
		
		// redirect user 
		String htmlPath = "/prism/tasks/" + sessionID + "/index.html";
		response.sendRedirect(htmlPath);
	}
	
	public String createSessionID() {
		// create a new session with random ID
		Date date = new Date();
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMMdd-kkmm");
		String prefix = dateFormat.format(date);
		
		String sessionID = prefix + "-" + (int) Math.floor(Math.random() * 1000000000); // convert int to string
		return sessionID;
	}
	
	public void executeHttpRequest(String id, Session session) 
			throws ClientProtocolException, IOException {
		HttpClient client = HttpClientBuilder.create().build();
        HttpGet httpget = new HttpGet(host + "/model/prism_result/id/" + id + "/document/raw_data");
        httpget.addHeader("Authorization", authToken);
        httpget.addHeader("accept", "application/json");
        HttpResponse response = client.execute(httpget);
        response.setHeader("Content-Type", "application/json");
        HttpEntity entity = response.getEntity();
        entity.writeTo(new FileOutputStream(new File(session.dir() + id + ".json")));        
	}
	
	public void registerSession(String sessionID, HttpServletRequest request) {
		// if an old session with the same id already exists, unregister it
		SessionManager sessionManager = SessionManager.getSessionManager();
		Session old = sessionManager.getSession(sessionID);
		if (old != null)
			sessionManager.removeSession(sessionID);
		
		// create new session and set ID, root, exception handler & heartbeat
		Session session = createNewSession(sessionID);

		SessionListener listener = new SessionListener();
		session.setListener(listener);
		
		String root = request.getServletContext().getRealPath(File.separator);
		session.setRoot(root);
		
		String dir = root + "tasks" + File.separator + sessionID + File.separator;
		session.setDir(dir);
		new File(dir).mkdir();
		
		ExceptionHandler exceptionHandler = new ExceptionHandler(listener);
		session.setExceptionHandler(exceptionHandler);

		setHeartBeat(session);
		
		// register the session
		sessionManager.addSession(sessionID, session);
		
		// make sure session directory exists
		File sessionDir = new File(session.dir());
		if (!sessionDir.isDirectory()) 
			sessionDir.mkdir();
		
		// set report
		PrismReport report = new PrismReport(session);
		session.setReport(report); 
		
		// set context name
		String context = request.getContextPath().replace("/", "");
		session.setContext(context);
	}

	public static Session createNewSession(String sessionID) {
		Session session = new BasicSession();
		session.setID(sessionID);
		return session;
	}
	
	public void setHeartBeat(Session session) {
		String heartbeat = TimeUtil.getTimeTag();
		session.setLastHeartBeat(heartbeat);
	}
	
}
