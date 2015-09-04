package ca.mcmaster.magarveylab.prism.web.session;

import javax.servlet.annotation.WebServlet;
import ca.mcmaster.magarveylab.wasp.session.SessionSubmit;

/**
 * Servlet to register a PRISM session. 
 * @author skinnider
 * 
 */
@WebServlet("/PrismSessionSubmit")
public class PrismSessionSubmit extends SessionSubmit {
	
	private static final long serialVersionUID = -5671961856226916741L;
	
}
