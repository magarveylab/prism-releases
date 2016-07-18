package ca.mcmaster.magarveylab.prism.util.exception;

/**
 * Thrown when a version is not in the format major.minor.patch (e.g. 1.3.1).
 * 
 * @author skinnider
 *
 */
public class VersionFormatException extends DependencyException {

	private static final long serialVersionUID = 3246955240801837146L;

	public VersionFormatException() {
		super();
	}

	public VersionFormatException(String message) {
		super(message);
	}

}
