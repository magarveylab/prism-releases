package ca.mcmaster.magarveylab.prism.util.exception;

/**
 * An error predicting genes within an input sequence with Prodigal.
 * 
 * @author skinnider
 *
 */
public class ProdigalSearchException extends Exception {

	private static final long serialVersionUID = 6636418602615180761L;

	public ProdigalSearchException() {
		super();
	}

	public ProdigalSearchException(String message) {
		super(message);
	}

}