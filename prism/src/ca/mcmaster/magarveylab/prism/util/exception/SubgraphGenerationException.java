package ca.mcmaster.magarveylab.prism.util.exception;

public class SubgraphGenerationException extends ScaffoldGenerationException {
	
	private static final long serialVersionUID = 709275824264401205L;

	public SubgraphGenerationException() {
		super();
	}
	
	public SubgraphGenerationException(String message) {
		super(message);
	}
	
	public SubgraphGenerationException(String m, Exception e) {
		super(m,e);
	}

}
