package ca.mcmaster.magarveylab.prism.util.exception;

public class NoBitsetException extends Exception {
	
	private static final long serialVersionUID = 709275824264401205L;

	public NoBitsetException() {
		super();
	}
	
	public NoBitsetException(String message) {
		super(message);
	}
	
	public NoBitsetException(String m, Exception e) {
		super(m,e);
	}

}