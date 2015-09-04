package ca.mcmaster.magarveylab.prism.util.exception;

public class SugarGeneException extends Exception {
	
	private static final long serialVersionUID = 709275824264401205L;

	public SugarGeneException() {
		super();
	}
	
	public SugarGeneException(String message) {
		super(message);
	}
	
	public SugarGeneException(String m, Exception e) {
		super(m,e);
	}

}
