package ca.mcmaster.magarveylab.prism.util.exception;

/**
 * Thrown when PRISM is not compatible with the version number of one of its
 * dependencies at runtime.
 * 
 * @author skinnider
 *
 */
public class VersionCompatibilityException extends DependencyException {

	private static final long serialVersionUID = 1187300452043305404L;

	public VersionCompatibilityException() {
		super();
	}

	public VersionCompatibilityException(String message) {
		super(message);
	}

}
