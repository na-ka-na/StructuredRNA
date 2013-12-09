package iiit;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.List;

import org.apache.commons.lang.StringUtils;

/**
 * 
 * This class generates Hbond info from the Hydrogen added PDB file and certain extra reference files.
 * 
 * @author sanchay.h AT gmail DOT com
 *
 */
public class HBonds {

	/* Reference data structures */
	static Map<String,List<String>> dictD_atoms = new HashMap<String,List<String>>();
	static Map<String,List<String>> dictA_atoms = new HashMap<String,List<String>>();
	static Map<String,List<String>> dictDH = new HashMap<String,List<String>>();
	static List<String> baseAtoms = new ArrayList<String>();
	static List<String> sugarAtoms = new ArrayList<String>();
	static List<String> backBoneAtoms = new ArrayList<String>();
	static List<String> ions = new ArrayList<String>();
	static Map<String,List<String>> WC = new HashMap<String,List<String>>();
	static Map<String,List<String>> H = new HashMap<String,List<String>>();
	static Map<String,List<String>> S = new HashMap<String,List<String>>();
	static List<Float> dDA = new ArrayList<Float>();
	static List<Float> dHA = new ArrayList<Float>();
	static List<Float> DHA = new ArrayList<Float>();
		
	// Respective dictionaries of donor, acceptor, hydrogen atoms.
	static List<Atom> donorList = new ArrayList<Atom>();
	static List<Atom> acceptorList = new ArrayList<Atom>();
	static Map<String,Atom> dictH = new HashMap<String, Atom>();
	
	/**
	 * Path of the cuminp.txt file
	 */
	public static String cuminpFileName = "C:\\Docs\\ramsagar\\cuminp.txt";
	
	/**
	 * Path of the cummap.txt file
	 */
	public static String cummapFileName = "C:\\Docs\\ramsagar\\cummap.txt";
	
	/**
	 * Path of the PDB file
	 */
	public static String pdbFileName    = "C:\\Docs\\ramsagar\\1VQO-min.coor";
	
	/**
	 * Path where the hbond file is to be persisted
	 */
	public static String HBondFileName  = "C:\\Docs\\ramsagar\\1VQO-min.coor.out";
	
	/**
	 * Main method
	 * 
	 * Does the following in order
	 * 
	 * <ol>
	 * <li> Read cuminp.txt and store reference info. </li>
	 * <li> Read cummap.txt and store reference info. </li>
	 * <li> Read the PDB file and form donor atom list, acceptor atom list, hydrogen atom list. </li>
	 * <li> Extracts Hbonds and writes them to a file. </li>
	 * </ol>
	 * 
	 * @param args
	 * @throws NumberFormatException
	 * @throws IOException
	 */
	public static void main(String[] args) throws NumberFormatException, IOException  {
		long t = -System.nanoTime();
		
		// Considers O of HOH , O5' and O3' in the code itself.
		System.out.println("Reading cuminp .. ");
		readCuminp(cuminpFileName);
		
		System.out.println("Reading cummap .. ");
		readCummap(cummapFileName);
		
		System.out.println("Reading PDB .. ");
		readPDB(pdbFileName);
		
		System.out.println("Writing Hbonds .. ");
		writeHBonds(HBondFileName);
		
		t += System.nanoTime();
		System.out.println("Done. time taken = " + t/1000000);
	}
	
	/**
	 * 
	 * This method does the actual calculation of Hbonds and prints each Hbond to a file.
	 * 
	 * Algorithm is quite simple - <br>
	 * For each donor we find acceptor atoms which satisfy <br> 
	 * <ol>
	 * <li> (distDA >= min_dDA) && (distDA <= max_dDA) </li>
	 * <li> (distHA >= min_dHA) && (distHA <= max_dHA) </li>
	 * <li> (DHAbondAngle >= min_DHA) && (DHAbondAngle <= max_DHA) </li>
	 * </ol>
	 * 
	 * @param hBondFileName
	 * @throws IOException
	 */
	public static void writeHBonds(String hBondFileName) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(hBondFileName));
		bw.write("Bond\t\tDonor Information\t\t\tAcceptor Information\t\tDistance\tAngle\tClassification");
		bw.write("\n\n");
		bw.write("H-Bond\tNo.D\tat\trt\trn\tchNo\tNo.A\tat\trt\trn\tchNo\tdDA\tdHA\taDHA\tty\tsub-ty");
		bw.write("\n\n");
		
		// Sort the acceptorList based on x-coordinate
		// See the compareTo function of the Atom class
		Collections.sort(acceptorList);
		
		final float min_dDA    = dDA.get(0);
		final float max_dDA    = dDA.get(1);
		
		final float min_DHA    = DHA.get(0);
		final float max_DHA    = DHA.get(1);
		
		final float min_dHA    = dHA.get(0);
		final float max_dHA    = dHA.get(1);
		
		Atom atom1 = new Atom();
		Atom atom2 = new Atom();
		
		/*
		 * Each donor and acceptor must satisfy the following criterion for forming a hbond
		 * 1. (distDA >= min_dDA) && (distDA <= max_dDA)
		 * 2. (distHA >= min_dHA) && (distHA <= max_dHA)
		 * 3. (DHAbondAngle >= min_DHA) && (DHAbondAngle <= max_DHA)
		 */
		for (Atom atomD : donorList){
			int flag = 0;
			atom1.x = atomD.x - max_dDA;
			atom2.x = atomD.x + max_dDA;
			
			int k1 = -Collections.binarySearch(acceptorList, atom1) - 1;
        	int k2 = -Collections.binarySearch(acceptorList, atom2) - 1;
			
			for (int k=k1; (k<=k2) && (k<acceptorList.size()); k++){
				Atom atomA = acceptorList.get(k);

				if (!atomA.residue_num.equals(atomD.residue_num)){
					
					float distDA = ((atomD.x - atomA.x) * (atomD.x - atomA.x))
								 + ((atomD.y - atomA.y) * (atomD.y - atomA.y))
								 + ((atomD.z - atomA.z) * (atomD.z - atomA.z));
					distDA = (float) Math.sqrt(distDA);
				
					if ((distDA >= min_dDA) && (distDA <= max_dDA)){
						if ("MG".equals(atomD.atom_name) || "POT".equals(atomD.atom_name)){
							String s = getHbondString(atomA, atomD, distDA, 0, 0)+"\n";
							bw.write(s);
							flag = 1;
						}
						else{
							for (String atomname : dictDH.get(atomD.atom_name)){
								String getKey = atomname + "$" + atomD.residue_num;
								if (dictH.containsKey(getKey)){
									Atom atomH = dictH.get(getKey);
									float distHA = ((atomH.x - atomA.x) * (atomH.x - atomA.x))
					         		 			 + ((atomH.y - atomA.y) * (atomH.y - atomA.y))
					         		 			 + ((atomH.z - atomA.z) * (atomH.z - atomA.z));
									distHA = (float) Math.sqrt(distHA);
									
									if ((distHA >= min_dHA) && (distHA <= max_dHA)){
									
										float distDH = ((atomD.x - atomH.x) * (atomD.x - atomH.x))
													 + ((atomD.y - atomH.y) * (atomD.y - atomH.y))
													 + ((atomD.z - atomH.z) * (atomD.z - atomH.z));
										distDH = (float) Math.sqrt(distDH);

										float bondAngle = (float) Math.acos( ((distDH*distDH)-(distDA*distDA)+(distHA*distHA))/(2*distDH*distHA) );
										bondAngle = (float) (bondAngle*(180*7.0/22.0));
									
										if ((bondAngle >= min_DHA) && (bondAngle <= max_DHA)){
											// Finally we know these donor-acceptor form hbond
											String s = getHbondString(atomA, atomD, distDA, distHA, bondAngle)+"\n";
											bw.write(s);
											flag = 1;
										}
									}
								}
							}
						}
					}
				}
			}
			if (flag == 1)
				bw.write("\n\n");
		}
		
		bw.close();
	}
	
	/**
	 * 
	 * Classifies a hbond and forms a String representation to print it. 
	 * 
	 * @param atomA
	 * @param atomD
	 * @param distDA
	 * @param distHA
	 * @param bondAngle
	 * @return String representing a hbond
	 */
	public static String getHbondString(Atom atomA, Atom atomD, float distDA,
			float distHA, float bondAngle) {
		String first = "";
		// To find out the type of Hydrogen bond checks for the presence of
		// carbon,nitrogen,oxygen in the donor and acceptor atoms
		if (atomD.atom_name.contains("C")) {
			first = "C-H";
		} else if (atomD.atom_name.contains("N")) {
			first = "N-H";
		} else if (atomD.atom_name.contains("O")) {
			first = "O-H";
		}
		if (atomA.atom_name.contains("N")) {
			first += "--N";
		}
		if (atomA.atom_name.contains("O")) {
			first += "--O";
		}

		String outstr;
		outstr = first + "\t" + atomD.atom_num + "\t" + atomD.atom_name + "\t";
		outstr += atomD.residue_name + "\t"
				+ Integer.parseInt(atomD.residue_num.substring(1).trim())
				+ "\t ";
		outstr += atomD.residue_num.substring(0, 1) + "\t";
		outstr += atomA.atom_num + "\t" + atomA.atom_name + "\t";
		outstr += atomA.residue_name + "\t"
				+ Integer.parseInt(atomA.residue_num.substring(1).trim())
				+ "\t ";
		outstr += atomA.residue_num.substring(0, 1) + "\t";
		outstr += round(distDA, 3) + "\t";
		outstr += round(distHA, 3) + "\t";
		outstr += round(bondAngle, 3) + "\t";

		// To find the type of atom that donor and acceptor belong to. whether a
		// base/backbone/sugar atom.
		if (baseAtoms.contains(atomD.atom_name)) {
			outstr += "b";
		} else if (backBoneAtoms.contains(atomD.atom_name)) {
			outstr += "B";
		} else if (sugarAtoms.contains(atomD.atom_name)) {
			outstr += "s";
		} else if (ions.contains(atomD.atom_name)) {
			outstr += "I";
		}

		if (baseAtoms.contains(atomA.atom_name)) {
			outstr += "-b" + "\t";
		} else if (backBoneAtoms.contains(atomA.atom_name)) {
			outstr += "-B" + "\t";
		} else if (sugarAtoms.contains(atomA.atom_name)) {
			outstr += "-s" + "\t";
		}

		// To find the sub-type of the hydrogen bond. whether it is watson crick
		// or etc...
		if (baseAtoms.contains(atomD.atom_name)
				&& baseAtoms.contains(atomA.atom_name)) {
			if (WC.get(atomD.residue_name).contains(atomD.atom_name)) {
				outstr += "WC" + atomD.residue_name;
			} else if (H.get(atomD.residue_name).contains(atomD.atom_name)) {
				outstr += "H" + atomD.residue_name;
			} else if (S.get(atomD.residue_name).contains(atomD.atom_name)) {
				outstr += "S" + atomD.residue_name;
			}
			if (WC.get(atomA.residue_name).contains(atomA.atom_name)) {
				outstr += "-WC" + atomA.residue_name;
			} else if (H.get(atomA.residue_name).contains(atomA.atom_name)) {
				outstr += "-H" + atomA.residue_name;
			} else if (S.get(atomA.residue_name).contains(atomA.atom_name)) {
				outstr += "-S" + atomA.residue_name;
			}
		}
		return outstr;
	}
	
	/**
	 * Reads the Hydrogen added PDB file.
	 * 
	 * Extracts the following - <br>
	 * <ol> 
	 * <li> List of acceptor atoms </li>
	 * <li> List of donor atoms </li>
	 * <li> Hydrogens corresponding to donors </li>
	 * </ol>
	 * 
	 * @param pdbFileName
	 * @throws IOException
	 */
	public static void readPDB(String pdbFileName) throws IOException{
		int count=0; // Respective counts of total atoms, donor, acceptor, hydrogen atoms.
		String keyH;

		float x,y,z;
		String atom_name;
		int atom_num;
		String residue_name;
		String residue_num;
		Atom a;
		
		BufferedReader br = new BufferedReader(new FileReader(pdbFileName));
		
		String line = br.readLine();
		while(line != null){
			
			// These are according to the pdb file format
			atom_num     = Integer.parseInt(line.substring(6,11).trim());
			atom_name    = line.substring(12,16).trim();
			residue_name = line.substring(17,20).trim();
			residue_num  = line.substring(21,26).trim();
			x            = Float.parseFloat(line.substring(30,38).trim());
			y            = Float.parseFloat(line.substring(38,46).trim());
			z            = Float.parseFloat(line.substring(46,54).trim());
			a            = new Atom(x,y,z,atom_name,atom_num,residue_name,residue_num);
			
			line = br.readLine();
			
			if (   dictD_atoms.containsKey(residue_name) 
				&& (   dictD_atoms.get(residue_name).contains(atom_name)
					|| dictD_atoms.get("O").contains(atom_name))){
				donorList.add(a);
			}
			if (   dictA_atoms.containsKey(residue_name) 
				&& (   dictA_atoms.get(residue_name).contains(atom_name)
					|| dictA_atoms.get("O").contains(atom_name))){
				acceptorList.add(a);
			}
			else if ("OH2".equals(atom_name) && "HOH".equals(residue_name)){
				donorList.add(a);
				acceptorList.add(a);
			}
			else if (atom_name.contains("H")){
				keyH = atom_name + "$" + residue_num;
				dictH.put(keyH, a);
			}
			else if ("MG".equals(atom_name) && "MG".equals(residue_name)){
				donorList.add(a);
			}
			else if ("POT".equals(atom_name) && "POT".equals(residue_name)){
				donorList.add(a);
			}
			
			count++;
		}
		br.close();
	}
	
	/**
	 * Reads cummap.txt
	 * 
	 * @param cummapFileName
	 * @throws IOException
	 */
	public static void readCummap(String cummapFileName) throws IOException{
		// Note : since , C & U bases , both have C5 one of them is not included in the mapping files as both result in the same.
		//It contains the mapping of which a respective donor interacts with a hydrogen atom. For e.g.: C2-H2 , C5'-H5' etc..

		BufferedReader br = new BufferedReader(new FileReader(cummapFileName));
		String line;
		String[] tmp;
		while((line = br.readLine()) != null){
			tmp = line.split("-");
			if (!dictDH.containsKey(tmp[0])){
				dictDH.put(tmp[0], new ArrayList<String>());
			}
			dictDH.get(tmp[0]).add(tmp[1]);
		}
		br.close();
	}
	
	/**
	 * Reads cuminp.txt 
	 * 
	 * @param cuminpFileName
	 * @throws NumberFormatException
	 * @throws IOException
	 */
	public static void readCuminp(String cuminpFileName) throws NumberFormatException, IOException{
		// The input criteria file consists of the criteria for finding the bond. Like, the Donor-Acceptor Bond distance (dDA) and Bond angle (DHA) etc..
		// It also contains residue specific base atoms for finding the sub-type of bond formed between bases. like DA/AA => Donor Adenine / Acceptor Adenine.
		// Similarly the list of base atoms, backbone atoms etc.. and also atoms specific to them SA => Sugar Adenine atoms

		BufferedReader br = new BufferedReader(new FileReader(cuminpFileName));
		String line;
		String[] s;
		while((line = br.readLine()) != null){
			if (line.contains("dDA")){
				s = (line.split("="))[1].split("-");
				for (int i=0; i<s.length; i++)
					dDA.add(Float.parseFloat(s[i]));
			}
			else if (line.contains("dHA")){
				s = (line.split("="))[1].split("-");
				for (int i=0; i<s.length; i++)
					dHA.add(Float.parseFloat(s[i]));
			}
			else if (line.contains("DHA")){
				s = (line.split("="))[1].split("-");
				for (int i=0; i<s.length; i++)
					DHA.add(Float.parseFloat(s[i]));
			}
			else if (line.contains("DA")){
				dictD_atoms.put("A",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("DG")){
				dictD_atoms.put("G",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("DC")){
				dictD_atoms.put("C",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("DU")){
				dictD_atoms.put("U",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("DO")){
				dictD_atoms.put("O",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("AA")){
				dictA_atoms.put("A",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("AG")){
				dictA_atoms.put("G",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("AC")){
				dictA_atoms.put("C",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("AU")){
				dictA_atoms.put("U",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("AO")){
				dictA_atoms.put("O",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("baseAtoms")){
				baseAtoms.addAll(Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("sugarAtoms")){
				sugarAtoms.addAll(Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("backBoneAtoms")){
				backBoneAtoms.addAll(Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("ions")){
				ions.addAll(Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("WCA")){
				WC.put("A",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("WCG")){
				WC.put("G",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("WCC")){
				WC.put("C",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("WCU")){
				WC.put("U",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("HA")){
				H.put("A",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("HG")){
				H.put("G",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("HC")){
				H.put("C",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("HU")){
				H.put("U",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("SA")){
				S.put("A",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("SG")){
				S.put("G",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("SC")){
				S.put("C",Arrays.asList((line.split("="))[1].split(",")));
			}
			else if (line.contains("SU")){
				S.put("U",Arrays.asList((line.split("="))[1].split(",")));
			}
		}
		br.close();
	}
	
	static float round(float value, int decimalPlace) {
		double power_of_ten = 1;
		// floating point arithmetic can be very tricky.
		// that's why I introduce a "fudge factor"
		double fudge_factor = 0.05;
		while (decimalPlace-- > 0) {
			power_of_ten *= 10.0d;
			fudge_factor /= 10.0d;
		}
		return (float) (Math.round((value + fudge_factor) * power_of_ten) / power_of_ten);
	}
	
	/**
	 * Atom class to hold each atom.
	 * 
	 * Has the following fields <br>
	 * <ul>
	 * <li>x</li>
	 * <li>y</li>
	 * <li>z</li>
	 * <li>atom_num</li>
	 * <li>atom_name</li>
	 * <li>residue_num</li>
	 * <li>residue_name</li>
	 * <ul> 
	 * 
	 * Implements Comparable with the compareTo method considering only the x coordinate.
	 */
	private static class Atom implements Comparable<Atom>{
		float x,y,z;
		String atom_name;
		int atom_num;
		String residue_name;
		String residue_num;
		
		Atom(float x, float y, float z, String atomName, int atomNum,
				String residueName, String residueNum) {
			super();
			this.x = x;
			this.y = y;
			this.z = z;
			atom_name = atomName;
			atom_num = atomNum;
			residue_name = residueName;
			residue_num = residueNum;
		}
		
		Atom() {}

		@Override
		public int compareTo(Atom o) {
			return (x >= o.x) ? 1 : -1;
		}
		
		@Override
		public String toString() {
			return StringUtils.join( new String[]{""+atom_num,atom_name,residue_num,residue_name,
					String.format("%.3f", x),String.format("%.3f", y),String.format("%.3f", z)}, 
					"|");
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + atom_num;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Atom other = (Atom) obj;
			if (atom_num != other.atom_num)
				return false;
			return true;
		}
		
	}
	
	
}



























