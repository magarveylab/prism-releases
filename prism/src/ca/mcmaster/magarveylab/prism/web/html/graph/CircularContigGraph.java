package ca.mcmaster.magarveylab.prism.web.html.graph;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.List;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import ca.mcmaster.magarveylab.enums.ClusterFamilies;
import ca.mcmaster.magarveylab.enums.Colors;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;

/**
 * 
 * A circular graph of all the biosynthetic clusters within a given contig. 
 * This is generated regardless of how many contigs are present
 * 
 * @author nishanth
 * 
 */
public class CircularContigGraph implements Serializable {
	
	private static final long serialVersionUID = 4239398439206827099L;
	
	private static final Color darkGrey = Color.decode(Colors.DARK_GREY.hex());
	private static final Color lightGrey = Color.decode(Colors.LIGHT_GREY.hex());
	private static final Color black = Color.decode(Colors.BLACK.hex());
	private static final Color white = Color.decode(Colors.WHITE.hex());
	
	private Contig contig;
	private Cluster cluster;
	private int size;
	private String link;

	/**
	 * Create a new circular genome graph for the entire genome (i.e., no single
	 * cluster will be graphed as the "currently selected cluster."
	 * 
	 * @param genome
	 *            genome to graph
	 * @param size
	 *            size of the graph, in pixels
	 */
	public CircularContigGraph(Contig contig, int size) {
		this.contig = contig;
		this.size = size;
	}

	/**
	 * Create a new circular genome graph with a single cluster selected to be
	 * displayed differently.
	 * 
	 * @param genome
	 *            genome to graph
	 * @param cluster
	 *            selected cluster
	 * @param size
	 *            size of the graph, in pixels
	 */
	public CircularContigGraph(Contig contig, Cluster cluster, int size) {
		this.contig = contig;
		this.cluster = cluster;
		this.size = size;
	}

	/**
	 * Get the HTML code to embed this genome graph on a PRISM web page.
	 */
	public String html(PrismConfig config) {
		StringBuffer sb = new StringBuffer();
		String graphClass = "circularGenomeGraph";
		
		sb.append("<object type='image/svg+xml' data='" + link + "' class='" + graphClass + "'>");
		sb.append("Your browser does not support SVG"); // fallback text
		sb.append("</object>");
	 
		String legend = Legend.getContigGraphLegend(contig, config);
		
		sb.append("<br>");
		
		sb.append(legend);
		
		return sb.toString();
	}
	
	/**
	 * Write this genome graph to an SVG file.
	 * 
	 * @param path
	 *            the location to write the SVG file to
	 * @throws Exception
	 */
	public void write(String path) throws Exception {
		int canvasSize = size + 10; // size + 5 px stroke x2
		
		SVGGraphics2D generator = getGenerator();
		Dimension svgCanvasSize = new Dimension(canvasSize, canvasSize);
		generator.setSVGCanvasSize(svgCanvasSize);
		
		drawChromosome(generator);
		drawArcs(generator);
		drawCenter(generator);
		
		writeToFile(generator, path);
		
		link = File.separator + "prism" + File.separator + path.split("prism" + File.separator)[1]; 
	}

	/**
	 * Draw the "chromosome" or circular canvas on which arcs are superimposed.
	 * 
	 * @param g2d
	 *            current Graphics2D object
	 */
	public void drawChromosome(Graphics2D g2d) {
		Rectangle2D bounds = new Rectangle2D.Double(5, 5, size, size);
		Arc2D chromosome = new Arc2D.Double(bounds, 100, 350, Arc2D.PIE);
		
		// stroke
		g2d.setPaint(darkGrey);
		g2d.setStroke(new BasicStroke(5.0f));
		g2d.draw(chromosome);
		
		// fill
		g2d.setPaint(lightGrey);
		g2d.fill(chromosome);
	}
	
	/**
	 * Draw arcs within the circular chromosome corresponding to identified
	 * clusters.
	 * 
	 * @param g2d
	 *            current Graphics2D object
	 */
	private void drawArcs(Graphics2D g2d) {
		
		int length = contig.length();
		
		for (Cluster cluster : contig.clusters()) {
			int clusterSize = cluster.size();
			int start = cluster.start();
			double offset = ((double) start / length * 350) + 100;
			double extent = (double) clusterSize / length * 350;
			
			List<ClusterFamilies> families = cluster.families();
			int numFamilies = families.size();
			System.out.println("[CircularGenomeGraph] Graphing cluster with " + numFamilies + " families");
			
			if (this.cluster != null && cluster == this.cluster) {
				Rectangle2D bounds = new Rectangle2D.Double(5, 5, size, size);
				Arc2D arc = new Arc2D.Double(bounds, offset, extent, Arc2D.PIE);
				
				g2d.setPaint(black);
				g2d.fill(arc);
				
			} else {
				int o = 5;
				int s = size;
				for (int i = 0; i < numFamilies; i++) {
					Rectangle2D bounds = new Rectangle2D.Double(o, o, s, s);
					Arc2D arc = new Arc2D.Double(bounds, offset, extent, Arc2D.PIE);
					
					ClusterFamilies family = families.get(i);
					g2d.setPaint(family.color());
					g2d.fill(arc);
					
					o += (50/numFamilies);
					s -= (100/numFamilies);
				}
			}
		}
	}
	
	/**
	 * Draw the white center of the genome graph above the arcs to finish the
	 * "chromosome."
	 * 
	 * @param g2d
	 *            current Graphics2D objet
	 */
	private void drawCenter(Graphics2D g2d) {
		int width = size * 3/4;
		int offset = size / 8 + 5;
		
		Rectangle2D bounds = new Rectangle2D.Double(offset, offset, width, width);
		Arc2D center = new Arc2D.Double(bounds, 100, 350, Arc2D.PIE);

		// stroke
		g2d.setPaint(darkGrey);
		g2d.setStroke(new BasicStroke(5.0f));
		g2d.draw(center);
		
		// fill 
		g2d.setPaint(white);
		g2d.fill(center);
		
		// fix for arc stroke
		Ellipse2D circle = new Ellipse2D.Double(offset, offset, width, width);
		g2d.setPaint(white);
		g2d.fill(circle);
	}

	/**
	 * Get the SVG generator.
	 * 
	 * @return the SVG generator
	 * @throws Exception
	 */
	private SVGGraphics2D getGenerator() throws Exception {
		DOMImplementation dom = GenericDOMImplementation.getDOMImplementation();
		Document doc = dom.createDocument(null, "svg", null);
		SVGGraphics2D generator = new SVGGraphics2D(doc);
		return generator;
	}
	
	/**
	 * Write an SVG graphic to a file.
	 * 
	 * @param generator
	 *            the SVG generator
	 * @param path
	 *            path to write the SVG file to
	 * @throws IOException
	 */
	private void writeToFile(SVGGraphics2D generator, String path) throws IOException {
		File file = new File(path);
		if (!file.exists())
			file.createNewFile();
		FileWriter fw = new FileWriter(file);
		PrintWriter writer = new PrintWriter(fw);
		generator.stream(writer);
		writer.close();
	} 
	
}
