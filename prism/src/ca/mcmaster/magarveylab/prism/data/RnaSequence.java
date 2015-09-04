package ca.mcmaster.magarveylab.prism.data;

import java.util.HashMap;

public class RnaSequence {
	
	Integer start;
	Integer end;
	char frame;
	String name;
	String product;
	HashMap<String, String> data;
	boolean sixteenSmallRna;
	
	public RnaSequence(Integer s, Integer e, String f, String n, String p){
		this.start = s;
		this.end = e;
		this.frame = f.charAt(0);
		this.name = n;
		this.product = p;
		
		if (name.equals("16S_rRNA")) {
			this.sixteenSmallRna = true;
		} else {
			this.sixteenSmallRna = false;
		}
	}
	
	public RnaSequence(Integer s, Integer e, HashMap<String, String> data, String strand){
		this.start = s;
		this.end = e;
		this.data = data;
		if (strand.equals("1")) {
			this.frame = '+';
		}else if (strand.equals("-1")) {
			this.frame = '-';
		}else {
			this.frame = '0';
		}
	}
	
	public Integer getStart() {
		return this.start;
	}
	
	public char getFrame() {
		return this.frame;
	}
	
	public Integer getEnd() {
		return this.end;
	}
	
	public String getName(){
		return this.name;
	}
	
	public String getProduct() {
		return this.product;
	}
	
	public HashMap<String, String> getData() {
		return this.data;
	}
	
	public void setData(HashMap<String, String> data) {
		this.data = data;
	}
	
	public boolean isSixteenSmallRna() {
		return this.sixteenSmallRna;
	}
	
}