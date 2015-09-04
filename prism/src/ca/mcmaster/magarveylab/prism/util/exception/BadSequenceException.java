package ca.mcmaster.magarveylab.prism.util.exception;

import java.io.IOException;

public class BadSequenceException extends IOException {

	private static final long serialVersionUID = 6636418602615180761L;
	
	public BadSequenceException() {
		super();
	}
	
	public BadSequenceException(String message) {
		super(message);
	}

}
