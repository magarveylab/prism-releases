package ca.mcmaster.magarveylab.prism.database;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import ca.mcmaster.magarveylab.prism.database.data.SmallMolecule;
import ca.mcmaster.magarveylab.prism.tanimoto.FingerprintUtil;
import ca.mcmaster.magarveylab.prism.util.Strings;
import ca.mcmaster.magarveylab.prism.util.exception.DatabaseConnectException;

public class AccessDatabase {
		
	/**
	 * Get the SQL connection object to Phil's BIRD database.
	 * @param url	URL of the BIRD database
	 * @return		database connection object
	 * @throws SQLException 
	 * @throws DatabaseConnectException 
	 */
	public static Connection getConnection(String url) throws Exception {
		Connection conn = null;
		try {
			Class.forName("org.neo4j.jdbc.Driver");
			conn = DriverManager.getConnection(url); //130.113.155.220
		} catch (ClassNotFoundException e) {
			throw new DatabaseConnectException("Couldn't find database drivers: " + e.getMessage(), e);
		} catch (Exception e) {
			throw e;
		}
		return conn;
	}

	/**
	 * Get the ResultSet corresponding to a cypher query.
	 * @param query	cipher query
	 * @param con	database connection
	 * @return		database ResultSet for given query 
	 */
	public static ResultSet getQuery(String query, Connection con) {
		ResultSet results = null;
		try {
			Statement all = con.createStatement();
			results = all.executeQuery(query);

		} catch (SQLException e) {
			System.out.println("Error querying database:" + e.getMessage());
		}
		return results;
	}
	
	public static List<SmallMolecule> getSmallMolecules(Connection conn) throws Exception {
		List<SmallMolecule> molecules = new ArrayList<SmallMolecule>();

		ResultSet results = AccessDatabase.getQuery("MATCH (s:Small_Molecule) "
				+ "return s.name, s.smile, s.finger_fcfp6, s.finger_ecfp6;", conn);

		try {
			while (results.next()) {
				SmallMolecule molecule = new SmallMolecule();

				// parse smiles
				String smiles = results.getString("s.smile");
				molecule.setSmiles(smiles);
				
				// parse names
				List<String> names = null;
				String name = null;
				Object object = results.getObject("s.name");
				if (object != null) {
					@SuppressWarnings("unchecked")
					List<String> array = (ArrayList<String>) object;
					names = array;
				}
				if (names != null) {
					name = Strings.arrayAsCommaDelinatedList(names);
					molecule.setName(name);
				} else {
					continue;
				}

				// parse fingerprints 
				String fcfp6Raw = results.getString("s.finger_fcfp6");
				String ecfp6Raw = results.getString("s.finger_ecfp6");
				if (fcfp6Raw != null && ecfp6Raw != null) {
					BitSet fcfp6 = FingerprintUtil.bitset(fcfp6Raw, 1024);
					BitSet ecfp6 = FingerprintUtil.bitset(ecfp6Raw, 1024);
					molecule.setFcfp6Fingerprint(fcfp6);
					molecule.setEcfp6Fingerprint(ecfp6);
				} else {
					continue;
				}
				
				molecules.add(molecule);
			}
		} catch (SQLException e) {
			System.out.println("Error querying database:" + e.getMessage());
		}
		
		return molecules;
	}
	
}