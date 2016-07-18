package ca.mcmaster.magarveylab.prism.util.exception;

public class ScaffoldGenerationException extends Exception {
	
	private static final long serialVersionUID = 709275824264401205L;

	public ScaffoldGenerationException() {
		super();
	}
	
	public ScaffoldGenerationException(String message) {
		super(message);
	}
	
	public ScaffoldGenerationException(String m, Exception e) {
		super(m,e);
	}
	
	public ScaffoldGenerationException(Exception e) {
		super(e);
	}

}
