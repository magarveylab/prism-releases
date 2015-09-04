package ca.mcmaster.magarveylab.prism.web;

import java.io.IOException;
import java.util.Timer;

import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.tomcat.util.http.fileupload.disk.DiskFileItemFactory;
import org.apache.tomcat.util.http.fileupload.servlet.ServletFileUpload;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.web.html.PrismReport;
import ca.mcmaster.magarveylab.wasp.WebApplicationSubmit;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Servlet for a genome NRPS/PKS search.
 * 
 * @author skinnider
 */
@WebServlet("/PrismSubmit")
public class PrismSubmit extends WebApplicationSubmit {

	private static final long serialVersionUID = -1846321311064350374L;

	@Override
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		DiskFileItemFactory fileItemFactory = new DiskFileItemFactory();
		ServletFileUpload uploadHandler = new ServletFileUpload(fileItemFactory);
		
		Session session = PrismSubmitParser.parseSession(request, uploadHandler);
		PrismReport report = (PrismReport) session.report();
		Timer timer = new Timer();
		timer.schedule(report, 0, 1000);

		try {
			PrismConfig config = PrismSubmitParser.parseConfig(request, uploadHandler, session);
			Prism prism = new Prism(config, session);
			session.setWebapp(prism);
			prism.run();
		} catch (Exception e) {
			session.listener().throwException(e);
		} finally {
			timer.cancel();
		}
	}

}
