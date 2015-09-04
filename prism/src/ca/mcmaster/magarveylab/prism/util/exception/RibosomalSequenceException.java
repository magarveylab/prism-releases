package ca.mcmaster.magarveylab.prism.util.exception;

/**
 * Indicates an error detecting ribosomal 16S sequence within user-input genome. 
 * @author skinnider
 *
 */
public class RibosomalSequenceException extends Exception {

	private static final long serialVersionUID = 6636418602615180761L;
	
	public RibosomalSequenceException() {
		super();
	}
	
	public RibosomalSequenceException(String message) {
		super(message);
	}

}
