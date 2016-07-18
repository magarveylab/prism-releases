package ca.mcmaster.magarveylab.prism.util.exception;

/**
 * Indicates an error running a Fimo search. 
 * @author Robyn Edgar
 *
 */
public class FimoSearchException extends Exception {

	private static final long serialVersionUID = 1L;

	public FimoSearchException() {
		super();
	}
	
	public FimoSearchException(String message) {
		super(message);
	}

}