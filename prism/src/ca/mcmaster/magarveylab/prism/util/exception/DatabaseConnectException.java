package ca.mcmaster.magarveylab.prism.util.exception;

public class DatabaseConnectException extends Exception {
	
	private static final long serialVersionUID = -4162530027934490284L;

	public DatabaseConnectException() {
		super();
	}
	
	public DatabaseConnectException(String message) {
		super(message);
	}
	
	public DatabaseConnectException(String m, Exception e) {
		super(m,e);
	}

}
