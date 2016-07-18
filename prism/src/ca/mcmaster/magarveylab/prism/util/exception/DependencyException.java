package ca.mcmaster.magarveylab.prism.util.exception;

/**
 * Thrown when PRISM cannot locate, or is not compatible with, one of its
 * dependencies.
 * 
 * @author skinnider
 *
 */
public class DependencyException extends Exception {

	private static final long serialVersionUID = 7168713672880499037L;

	public DependencyException() {
		super();
	}

	public DependencyException(String message) {
		super(message);
	}

}
