package ca.mcmaster.magarveylab.prism.web;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

import javax.servlet.ServletInputStream;
import javax.servlet.http.HttpServletRequest;

import org.apache.tomcat.util.http.fileupload.FileItem;
import org.apache.tomcat.util.http.fileupload.servlet.ServletFileUpload;
import org.apache.tomcat.util.http.fileupload.servlet.ServletRequestContext;

import ca.mcmaster.magarveylab.prism.orfs.GenePredictionModes;
import ca.mcmaster.magarveylab.prism.web.html.PrismReport;
import ca.mcmaster.magarveylab.wasp.session.FileUploadListener;
import ca.mcmaster.magarveylab.wasp.session.Session;
import ca.mcmaster.magarveylab.wasp.session.SessionManager;

public class PrismSubmitParser {
	
	/**
	 * Parse a HttpServletRequest to instantiate a new PrismConfig configuration object for this session.
	 * @param request	the request to parse
	 * @param uploadHandler	the request's upload handler
	 * @param session	the current session
	 * @return	the configuration object
	 * @throws Exception
	 */
	public static PrismConfig parseConfig(HttpServletRequest request, ServletFileUpload uploadHandler, Session session) 
			throws Exception {
		PrismConfig config = new PrismConfig();
		config.web = true;
		config.saveSequences = true;
		
		Date date = new Date();
		config.date = date;
		
		List<FileItem> items = uploadHandler.parseRequest(new ServletRequestContext(request));
		Iterator<FileItem> itr = items.iterator();
		// iterate over form items
		while (itr.hasNext()) {
			FileItem item = (FileItem) itr.next();
			if (item.isFormField()) {
				// handle form items
				switch (item.getFieldName()) {
					case "allOrfs":
						System.out.println("[PrismSubmitParser] Finding all potential coding sequences");
						config.genePredictionModes.add(GenePredictionModes.ALL_ORFS);
						break;
					case "prodigal": 
						System.out.println("[PrismSubmitParser] Using Prodigal to predict orfs");
						config.genePredictionModes.add(GenePredictionModes.PRODIGAL);
						break;
					case "genemark":
						config.genePredictionModes.add(GenePredictionModes.GENEMARK);
						break;
					case "maxSubstrates":
						int maxSubstrates = Integer.parseInt(item.getString());
						if (maxSubstrates > 0)
							config.display = maxSubstrates;
						System.out.println("[PrismSubmitParser] Set maximum number of "
								+ "substrates to display to " + maxSubstrates);
						break;
					case "window":
						int window = Integer.parseInt(item.getString());
						if (window > 0)
							config.window = window;
						System.out.println("[PrismSubmitParser] Set window to " + window);
						break;
					case "score":
						config.score = true;
						break;
					case "thiotemplated":
						config.thiotemplated = true;
						break;
					case "sugar":
						config.sugar = true;
						break;
					case "resistance":
						config.resistance = true;
						break;
					case "regulation":
						config.regulation = true;
						break;
					case "ribosomal":
						config.ribosomal = true;
						break;
					case "tanimotoCutoff":
						double tanimotoCutoff = Double.parseDouble(item.getString());
						if (tanimotoCutoff > 0)
							config.tanimotoCutoff = tanimotoCutoff;
						System.out.println("[PrismSubmitParser] Set tanimoto cutoff to " + tanimotoCutoff);
						break;
					case "homologyCutoff":
						double homologyCutoff = Double.parseDouble(item.getString());
						if (homologyCutoff > 0)
							config.homologyCutoff = homologyCutoff;
						System.out.println("[PrismSubmitParser] Set homology cutoff to " + homologyCutoff);
						break;
					case "scaffoldLimit":
						int scaffoldLimit = Integer.parseInt(item.getString());
						config.scaffoldLimit = scaffoldLimit;
						System.out.println("[PrismSubmitParser] Set scaffold limit to " + scaffoldLimit);
						break;
					default:
						break;
				}
			} else { 
				// handle genome
				File localFile = new File(item.getName());
				String webFileName = localFile.getName();
				String webSafeFileName = webFileName.replaceAll(" ", "_");
				
				// write genome file to session directory
				File genomeFile = new File(session.dir() + webSafeFileName);
				if (!genomeFile.exists())
					item.write(genomeFile);
				System.out.println("[PrismSubmitParser] Wrote incoming file to " + genomeFile.getAbsolutePath());

				if (item.getFieldName().equals("genomeFile")) {
					config.input = genomeFile.getAbsolutePath();
				}
			}
		}
		
		return config;
	}
	
	/**
	 * Parse a HttpServletRequest to create a new session and set its upload handler. 
	 * @param request	the request to parse
	 * @param uploadHandler	the upload handler to associate with the session
	 * @return	the new session
	 * @throws IOException
	 */
	public static Session parseSession(HttpServletRequest request, ServletFileUpload uploadHandler) throws IOException {
		SessionManager sessionManager = SessionManager.getSessionManager();
		Session session = null;

		ServletInputStream in = request.getInputStream();
		byte b[] = new byte[256];
		for (int i = 0; i < 4; i++) {
			Arrays.fill(b, (byte) 0);
			in.readLine(b, 0, 256);
		}
		String sessionID = new String(b, 0, 256);
		sessionID = sessionID.trim();
		System.out.println("[PrismSubmitParser] Session: " + sessionID);	
		session = sessionManager.getSession(sessionID);

		// make sure session directory exists
		File sessionDir = new File(session.dir());
		if (!sessionDir.isDirectory()) 
			sessionDir.mkdir();
		
		// set upload progress listener
		FileUploadListener uploadListener = new FileUploadListener();
		uploadHandler.setProgressListener(uploadListener);
		session.listener().uploadListener = uploadListener;
		
		// set report
		PrismReport report = new PrismReport(session);
		session.setReport(report); 
		report.run();
		
		// set context name
		String context = request.getContextPath().replace("/", "");
		session.setContext(context);
		
		return session;
	}

}
