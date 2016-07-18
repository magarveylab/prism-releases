package ca.mcmaster.magarveylab.prism.web;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Timer;

import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.tomcat.util.http.fileupload.FileItem;
import org.apache.tomcat.util.http.fileupload.disk.DiskFileItemFactory;
import org.apache.tomcat.util.http.fileupload.servlet.ServletFileUpload;
import org.apache.tomcat.util.http.fileupload.servlet.ServletRequestContext;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.database.JsonInput;
import ca.mcmaster.magarveylab.prism.util.PrismFileWriter;
import ca.mcmaster.magarveylab.prism.web.html.PrismReport;
import ca.mcmaster.magarveylab.wasp.WebApplicationSubmit;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Servlet to load saved results of a PRISM search, in JSON form. 
 * 
 * @author skinnider
 */
@WebServlet("/PrismJsonSubmit")
public class PrismJsonSubmit extends WebApplicationSubmit {

	private static final long serialVersionUID = -1846321311064350374L;

	@Override
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		DiskFileItemFactory fileItemFactory = new DiskFileItemFactory();
		ServletFileUpload uploadHandler = new ServletFileUpload(fileItemFactory);
		
		Session session = PrismSubmitParser.parseSession(request, uploadHandler);
		PrismReport report = (PrismReport) session.report();
		Timer timer = new Timer();
		timer.schedule(report, 0, 1_000);

		try {
			String filepath = parseJsonFile(request, uploadHandler, session);
			Prism prism = JsonInput.read(filepath, session);
			session.setWebapp(prism);
			PrismFileWriter.writeAllFiles(prism.genome(), session);
			prism.terminate();
		} catch (Exception e) {
			session.listener().throwException(e);
		} finally {
			timer.cancel();
		}
	}

	public String parseJsonFile(HttpServletRequest request, ServletFileUpload uploadHandler, Session session) 
			throws Exception {
		List<FileItem> items = uploadHandler.parseRequest(new ServletRequestContext(request));
		Iterator<FileItem> itr = items.iterator();
		
		String jsonFilepath = null; 
		
		// iterate over form items
		while (itr.hasNext()) {
			FileItem item = (FileItem) itr.next();
			if (!item.isFormField()) {
				// handle file
				File localFile = new File(item.getName());
				String webFileName = localFile.getName();
				String webSafeFileName = webFileName.replaceAll(" ", "_");
				
				// write file to session directory
				File jsonFile = new File(session.dir() + webSafeFileName);
				if (!jsonFile.exists())
					item.write(jsonFile);
				System.out.println("[PrismJsonSubmit] Wrote incoming file to " + jsonFile.getAbsolutePath());
				jsonFilepath = jsonFile.getAbsolutePath();
			}
		}
		
		return jsonFilepath;
	}
	
}
