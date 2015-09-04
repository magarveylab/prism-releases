package ca.mcmaster.magarveylab.prism.util;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.generators.IGenerator;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;

/**
 * Generates images of SMILES structures.
 * @author skinnider
 *
 */
public class SmilesImageGenerator {
	
	/**
	 * Generate a PNG image of a SMILES structure.
	 * @param smiles		the SMILES to graph 
	 * @param out			location to write the file to
	 * @param width			width of the image
	 * @param height		height of the image
	 * @throws CDKException
	 * @throws IOException
	 */
	public static void draw(String smiles, String out, int width, int height) throws CDKException, IOException {
		// set draw area and image size
		Rectangle drawArea = new Rectangle(width, height);
		Image image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		
		// set molecule
		IAtomContainer molecule = SmilesIO.molecule(smiles);
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule((IMolecule) molecule);
	//	sdg.setBondLength(1.0);
		sdg.generateCoordinates();
		molecule = sdg.getMolecule();
		
		// generators create image elements
		List<IGenerator<IAtomContainer>> generators = new ArrayList<IGenerator<IAtomContainer>>();
		generators.add(new BasicSceneGenerator());
		generators.add(new BasicBondGenerator());
		generators.add(new BasicAtomGenerator());
		
		// set up renderer
		AtomContainerRenderer renderer = new AtomContainerRenderer(generators, new AWTFontManager());
		renderer.setup(molecule, drawArea);
		
		// black & white
		RendererModel model = renderer.getRenderer2DModel();
		model.set(BasicAtomGenerator.ColorByType.class, false);
				
		// paint background
		Graphics2D g2d = (Graphics2D) image.getGraphics();
		g2d.setColor(Color.WHITE);
		g2d.fillRect(0, 0, width, height);
		
		// paint molecule w/ toolkit-specific renderer
		renderer.paint(molecule, new AWTDrawVisitor(g2d), drawArea, true);
		
		// write file
		ImageIO.write((RenderedImage) image, "PNG", new File(out));
	}

}
