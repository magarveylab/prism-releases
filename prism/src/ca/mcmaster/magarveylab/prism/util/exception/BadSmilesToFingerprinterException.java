package ca.mcmaster.magarveylab.prism.util.exception;

public class BadSmilesToFingerprinterException extends Exception {

	private static final long serialVersionUID = 1L;

	public BadSmilesToFingerprinterException(String smiles, int fingerprintCount) {
		System.out.println("Couldn't make fingerprints: only " + fingerprintCount + " were completed for: " + smiles);
	}
}
