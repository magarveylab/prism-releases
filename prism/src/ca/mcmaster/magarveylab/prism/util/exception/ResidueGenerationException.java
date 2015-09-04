package ca.mcmaster.magarveylab.prism.util.exception;

public class ResidueGenerationException extends ScaffoldGenerationException {
	
	private static final long serialVersionUID = 709275824264401205L;

	public ResidueGenerationException() {
		super();
	}
	
	public ResidueGenerationException(String message) {
		super(message);
	}
	
	public ResidueGenerationException(String m, Exception e) {
		super(m,e);
	}

}
