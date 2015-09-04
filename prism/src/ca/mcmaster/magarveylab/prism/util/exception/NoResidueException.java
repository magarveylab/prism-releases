package ca.mcmaster.magarveylab.prism.util.exception;

public class NoResidueException extends ScaffoldGenerationException {
	
	private static final long serialVersionUID = 709275824264401205L;

	public NoResidueException() {
		super();
	}
	
	public NoResidueException(String message) {
		super(message);
	}
	
	public NoResidueException(String m, Exception e) {
		super(m,e);
	}

}
