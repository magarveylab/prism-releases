package ca.mcmaster.magarveylab.prism.util.exception;

public class BondFormationException extends ScaffoldGenerationException {
	
	private static final long serialVersionUID = 709275824264401205L;

	public BondFormationException() {
		super();
	}
	
	public BondFormationException(String message) {
		super(message);
	}
	
	public BondFormationException(String m, Exception e) {
		super(m,e);
	}

}
