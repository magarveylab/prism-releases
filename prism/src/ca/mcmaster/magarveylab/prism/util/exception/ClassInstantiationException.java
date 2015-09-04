package ca.mcmaster.magarveylab.prism.util.exception;

public class ClassInstantiationException extends Exception {
	
	private static final long serialVersionUID = 709275824264401205L;

	public ClassInstantiationException() {
		super();
	}
	
	public ClassInstantiationException(String message) {
		super(message);
	}
	
	public ClassInstantiationException(String m, Exception e) {
		super(m,e);
	}

}
