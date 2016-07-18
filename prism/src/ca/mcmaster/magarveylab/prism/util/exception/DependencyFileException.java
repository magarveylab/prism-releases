package ca.mcmaster.magarveylab.prism.util.exception;

/**
 * Thrown when a file required to execute a dependency (i.e., a bash script
 * required to run an installed binary from Java) cannot be found.
 * 
 * @author skinnider
 *
 */
public class DependencyFileException extends DependencyException {

	private static final long serialVersionUID = -5058206864932881409L;

	public DependencyFileException() {
		super();
	}

	public DependencyFileException(String message) {
		super(message);
	}

}
