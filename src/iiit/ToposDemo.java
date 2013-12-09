package iiit;

import iiit.RNAAsGraph;

import java.util.Iterator;

import clojure.lang.*;

/**
 * Demo Class for {@link iiit.Topos iiit.Topos} 
 * 
 * @author Sanchay Harneja
 *
 */
public class ToposDemo {

	/**
	 * See source code for demos of
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
				
		final String atomsPDBFileName1 = "C:\\Docs\\ramsagar\\Input.pdb";
		final String hBondsFileName1   = "C:\\Docs\\ramsagar\\Input.pdb.out.grepped";
		
		final String atomsPDBFileName2 = "C:\\Docs\\ramsagar\\1VQO-min.coor";
		final String hBondsFileName2   = "C:\\Docs\\ramsagar\\1VQO-min.coor.out.grepped";
		
		RNAAsGraph RAG; 
		RAG = demo_formRNAGraph(atomsPDBFileName1, hBondsFileName1);
		//*
		demo_getRNAGraphStats(RAG);
		demo_getAllResidues(RAG);
		demo_getAllCBonds(RAG);
		demo_getAllHBonds(RAG);
		demo_printRNAGraphToFile(RAG, "RAG1.txt");
		demo_allToposCounts(RAG);
		demo_getNLets(RAG, 24);
		demo_getStems(RAG, 2);
		demo_getStemLoops(RAG, 4, 2);
		demo_getBulges(RAG,5,5,3);
		
		RAG = demo_formRNAGraph(atomsPDBFileName2, hBondsFileName2);
		demo_getRNAGraphStats(RAG);
		demo_getAllResidues(RAG);
		demo_getAllCBonds(RAG);
		demo_getAllHBonds(RAG);
		demo_printRNAGraphToFile(RAG, "RAG2.txt");
		demo_allToposCounts(RAG);
		demo_getNLets(RAG, 33);
		demo_getStems(RAG, 16);
		demo_getStemLoops(RAG, 5, 3);
		demo_getBulges(RAG,4,5,3);
		demo_getKissingStemLoops(RAG, 6, 6, 6, 7);
		//*/
	}
	
	/**
	 * Prints down counts of various kinds of topologies in a RAG all at once.
	 * 
	 * @param RAG
	 */
	public static void demo_allToposCounts(RNAAsGraph RAG){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_allToposCounts");
		
		// N-lets
		for (int n=2; n<=40; n++){
			int count = Topos.getNLets(RAG, n).count();
			if (count != 0)
				System.out.println("Number of "+n+"-lets = " + count);
		}
		System.out.println();
		
		// Stems
		for (int s=2; s<=30; s++){
			int count = Topos.getStems(RAG, s).count();
			if (count != 0)
				System.out.println("Number of Stems("+s+")= " + count);
		}
		System.out.println();
		
		// Stem loops
		for (int l=1; l<20; l++){
			for (int s=1; s<15; s++){
				int count = Topos.getStemLoops(RAG, s, l).count();
				if (count != 0)
					System.out.println("Number of StemLoops("+s+","+l+")= " + count);
			}
		}
		System.out.println();
		
		// Bulges
		for (int s1=2; s1<16; s1++){
			for (int s2=s1; s2<16; s2++){
				for (int b=1; b<20; b++){
					int count = Topos.getBulges(RAG, s1, s2, b).count();
					if (count != 0)
						System.out.println("Number of Bulges("+s1+","+s2+","+b+")= " + count);
				}
			}
		}
		System.out.println();
		
		// Kissing Stem loops
		for (int s1=2; s1<20; s1++){
			for (int l1=1; l1<20; l1++){
				for (int s2=s1; s2<20; s2++){
					for (int l2=l1; l2<20; l2++){
						int count = Topos.getKissingStemLoops(RAG, s1, l1, s2, l2).count();
						if (count != 0)
							System.out.println("Number of KissingStemLoops("+s1+","+l1+","+s2+","+l2+")= " + count);
					}
				}
			}
		}
		System.out.println();
	}
	
	/**
	 * Form a RAG
	 * 
	 * @param atomsPDBFileName
	 * @param hBondsFileName
	 * @return
	 */
	public static RNAAsGraph demo_formRNAGraph(String atomsPDBFileName, String hBondsFileName){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_formRNAGraph");
		Keyword atomsPDBFileNameKey = Keyword.intern(Symbol.create("atoms-pdb-file-name"));
		Keyword hBondsFileNameKey   = Keyword.intern(Symbol.create("hbonds-file-name"));
				
		PersistentArrayMap p = new PersistentArrayMap(
				new Object[]{atomsPDBFileNameKey, atomsPDBFileName, 
						     hBondsFileNameKey,   hBondsFileName}
				);
		
		RNAAsGraph RAG = Topos.formRNAGraph(p);
		return RAG;
	}

	/**
	 * Get some RAG stats
	 * 
	 * @param RAG
	 */
	public static void demo_getRNAGraphStats(RNAAsGraph RAG){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_getRNAGraphStats");
		PersistentArrayMap stats = Topos.getRNAGraphStats(RAG);
		for (Iterator itr = stats.keySet().iterator(); itr.hasNext();){
			Object key = itr.next();
			System.out.println(key + " => " + stats.get(key));
		}
	}
	
	/**
	 * Get all residue nodes comprising RAG
	 * 
	 * @param RAG
	 */
	public static void demo_getAllResidues(RNAAsGraph RAG){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_getAllResidues");
		PersistentHashMap residues = Topos.getAllResidues(RAG);
		System.out.println("Number of residue nodes: "+residues.count());
		for (Iterator itr = residues.keySet().iterator(); itr.hasNext();){
			String resName = (String) itr.next();
			Residue res = (Residue) residues.get(resName);
			System.out.println(resName + " => [id=" + res.id + "; name=" + res.name + "; chain=" + res.chain + "; type=" + res.type + "]");
		}
	}
	
	/**
	 * Get all cbond edges comprising RAG
	 * 
	 * @param RAG
	 */
	public static void demo_getAllCBonds(RNAAsGraph RAG){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_getAllCBonds");
		PersistentHashMap cbonds = Topos.getAllCBonds(RAG);
		for (Iterator itr = cbonds.keySet().iterator(); itr.hasNext();){
			String resName = (String) itr.next();
			PersistentVector v = (PersistentVector) cbonds.get(resName);
			String prevCbondResName = (String) v.get(0);
			String nextCbondResName = (String) v.get(1);
			System.out.println(resName + " => [prevCbond=" + prevCbondResName + "; nextCbond=" + nextCbondResName+"]");
		}
	}
	
	/**
	 * Get all hbond edges comprising RAG
	 * 
	 * @param RAG
	 */
	public static void demo_getAllHBonds(RNAAsGraph RAG){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_getAllHBonds");
		PersistentHashMap hbonds = Topos.getAllHBonds(RAG);
		for (Iterator itr = hbonds.keySet().iterator(); itr.hasNext();){
			String resName = (String) itr.next();
			PersistentHashSet s = (PersistentHashSet) hbonds.get(resName);
			System.out.print(resName + " => #{ ");
			for (Iterator itr1 = s.iterator(); itr1.hasNext();){
				String resName1 = (String) itr1.next();
				System.out.print(resName1+ ", ");
			}
			System.out.println("}");
		}
		
	}
	
	/**
	 * Prints the RAG graph to a file
	 * 
	 * @param RAG
	 * @param fileName
	 */
	public static void demo_printRNAGraphToFile(RNAAsGraph RAG, String fileName){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_printRNAGraphToFile");
		Topos.printRNAGraphToFile(RAG, fileName);
	}
	
	/**
	 * Demo for n-lets
	 * 
	 * @param RAG
	 * @param n
	 */
	public static void demo_getNLets(RNAAsGraph RAG, int n){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_getNLets");
		PersistentHashSet nlets = Topos.getNLets(RAG, n);
		int i=1;
		for (Iterator itr = nlets.iterator(); itr.hasNext(); ){
			ArraySeq nlet = (ArraySeq) itr.next();
			System.out.print(n+"-let["+i+"] = [ ");
			for (Iterator itr1 = nlet.iterator(); itr1.hasNext();){
				String res = (String) itr1.next();
				System.out.print(res+", ");
			}
			System.out.println("]");
			i++;
		}
	}
	
	private static void printResidues(Iterator itr1){
		System.out.print("[ ");
		for (;itr1.hasNext();) {
			String res = (String) itr1.next();
			System.out.print(res + ", ");
		}
		System.out.println("]");
	}
	
	private static void printStem(PersistentVector stem){
		PersistentVector stemLeft = (PersistentVector) stem.get(0);
		PersistentVector stemRight = (PersistentVector) stem.get(1);
		System.out.print("LeftPart = ");
		printResidues(stemLeft.iterator());
		System.out.print("RightPart = ");
		printResidues(stemRight.iterator());
	}
	
	/**
	 * Demo for stems
	 * 
	 * @param RAG
	 * @param stemLen
	 */
	public static void demo_getStems(RNAAsGraph RAG, int stemLen){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_getStems");
		LazySeq stems = Topos.getStems(RAG,stemLen);
		int i=1;
		for (Iterator itr = stems.iterator(); itr.hasNext();){
			System.out.println("stem-"+stemLen+"["+i+"] = [ ");
			PersistentVector stem = (PersistentVector) itr.next();
			printStem(stem);
			System.out.println("]");
			i++;
		}
	}
	
	private static void printStemLoop(PersistentVector stemLoop){
		PersistentVector stem = (PersistentVector)stemLoop.get(0);
		PersistentVector loop = (PersistentVector)stemLoop.get(1);
		printStem(stem);
		System.out.print("Loop = ");
		printResidues(loop.iterator());
	}
	
	/**
	 * Demo for stem loops
	 * 
	 * @param RAG
	 * @param stemLen
	 * @param loopLen
	 */
	public static void demo_getStemLoops(RNAAsGraph RAG, int stemLen, int loopLen){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_getStemLoops");
		LazySeq stemLoops = Topos.getStemLoops(RAG, stemLen, loopLen);
		int i=1;
		for (Iterator itr = stemLoops.iterator(); itr.hasNext();){
			System.out.println("stemLoop-("+stemLen+","+loopLen+")["+i+"] = [ ");
			PersistentVector stemLoop = (PersistentVector) itr.next();
			printStemLoop(stemLoop);
			System.out.println("]");
			i++;
		}
	}
	
	/**
	 * Demo for bulges
	 * 
	 * @param RAG
	 * @param stem1Len
	 * @param stem2Len
	 * @param bulgeLen
	 */
	public static void demo_getBulges(RNAAsGraph RAG, int stem1Len, int stem2Len, int bulgeLen){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_getBulges");
		LazySeq stemLoops = Topos.getBulges(RAG, stem1Len, stem2Len, bulgeLen);
		int i=1;
		for (Iterator itr = stemLoops.iterator(); itr.hasNext();){
			System.out.println("bulge-("+stem1Len+","+stem2Len+","+bulgeLen+")["+i+"] = [ ");
			PersistentVector bulge = (PersistentVector) itr.next();
			PersistentVector stem1 = (PersistentVector)bulge.get(0);
			PersistentVector stem2 = (PersistentVector)bulge.get(1);
			PersistentVector bulgeLoop = (PersistentVector)bulge.get(2);
			printStem(stem1);
			printStem(stem2);
			System.out.print("BulgeLoop = ");
			printResidues(bulgeLoop.iterator());
			System.out.println("]");
			i++;
		}
	}
	
	/**
	 * Demo for Kissing loops
	 * 
	 * @param RAG
	 * @param stem1Len
	 * @param loop1Len
	 * @param stem2Len
	 * @param loop2Len
	 */
	public static void demo_getKissingStemLoops(RNAAsGraph RAG, int stem1Len, int loop1Len, int stem2Len, int loop2Len){
		System.out.println("--------------------------------------------------------------------");
		System.out.println("Inside demo_getKissingStemLoops");
		LazySeq kissingStemLoops = Topos.getKissingStemLoops(RAG, stem1Len, loop1Len, stem2Len, loop2Len);
		int i=1;
		for (Iterator itr = kissingStemLoops.iterator(); itr.hasNext();){
			System.out.println("kissingStemLoop-("+stem1Len+","+loop1Len+","+stem2Len+","+loop2Len+")["+i+"] = [ ");
			PersistentVector kissingStemLoop = (PersistentVector) itr.next();
			PersistentVector stemLoop1 = (PersistentVector) kissingStemLoop.get(0);
			PersistentVector stemLoop2 = (PersistentVector) kissingStemLoop.get(1);
			PersistentHashSet hBondIntersectionRes = (PersistentHashSet) kissingStemLoop.get(2);
			printStemLoop(stemLoop1);
			printStemLoop(stemLoop2);
			System.out.print("hBondIntersectionRes in StemLoop2 = ");
			printResidues(hBondIntersectionRes.iterator());
			System.out.println("]");
			i++;
		}
	}
	
}
