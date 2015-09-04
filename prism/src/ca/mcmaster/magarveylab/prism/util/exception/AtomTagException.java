package ca.mcmaster.magarveylab.prism.util.exception;

public class AtomTagException extends ScaffoldGenerationException {
	
	private static final long serialVersionUID = 709275824264401205L;

	public AtomTagException() {
		super();
	}
	
	public AtomTagException(String message) {
		super(message);
	}
	
	public AtomTagException(String m, Exception e) {
		super(m,e);
	}

}
