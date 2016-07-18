package ca.mcmaster.magarveylab.prism.util.exception;

/**
 * Thrown when an binary dependency is not installed.
 * 
 * @author skinnider
 *
 */
public class DependencyBinaryException extends DependencyException {

	private static final long serialVersionUID = -1025978320935769218L;

	public DependencyBinaryException() {
		super();
	}

	public DependencyBinaryException(String message) {
		super(message);
	}

}
