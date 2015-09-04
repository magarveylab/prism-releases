package ca.mcmaster.magarveylab.prism.util.exception;

public class TailoringSubstrateException extends Exception {
	
	private static final long serialVersionUID = 709275824264401205L;

	public TailoringSubstrateException() {
		super();
	}
	
	public TailoringSubstrateException(String message) {
		super(message);
	}
	
	public TailoringSubstrateException(String m, Exception e) {
		super(m,e);
	}

}
