package ca.mcmaster.magarveylab.prism.util.exception;

import java.io.IOException;

public class NoSequenceException extends IOException {
	
	private static final long serialVersionUID = 3136531587233621335L;

	public NoSequenceException() {
		super();
	}
	
	public NoSequenceException(String message) {
		super(message);
	}
	
	public NoSequenceException(String m, Exception e) {
		super(m,e);
	}

}
