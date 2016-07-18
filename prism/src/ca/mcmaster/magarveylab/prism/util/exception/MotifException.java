package ca.mcmaster.magarveylab.prism.util.exception;

public class MotifException extends Exception {
	
	private static final long serialVersionUID = 709275824264401205L;

	public MotifException() {
		super();
	}
	
	public MotifException(String message) {
		super(message);
	}
	
	public MotifException(String m, Exception e) {
		super(m,e);
	}

}
