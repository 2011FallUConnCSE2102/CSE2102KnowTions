/*
 * Created on Jan 24, 2010
 */
package org.knowceans.corpus;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.knowceans.util.ArrayUtils;
import org.knowceans.util.Print;

/**
 * CorpusResolver resolves indices into names.
 * 
 * @author gregor
 */
public class CorpusResolver {

	public static void main(String[] args) {
		CorpusResolver cr = new CorpusResolver("corpus-example/nips");
		System.out.println(cr.getLabel(20));
		System.out.println(cr.getDoc(501));
		System.out.println(cr.getTerm(1));
		System.out.println(cr.getTermId(cr.getTerm(1)));
	}

	/**
	 * these indices correspond to ILabelCorpus.L* + 2 (docs and terms are -2
	 * and -1 there)
	 */
	public final String[] EXTENSIONS = { "docs", "vocab", "authors.key",
			"labels.key", "vols.key", "docnames", "docs.key" };

	HashMap<String, Integer> termids;
	String[][] data = new String[EXTENSIONS.length][];
	String filebase;

	private boolean parmode;

	public CorpusResolver(String filebase) {
		this(filebase, false);
	}

	/**
	 * control paragraph mode (possibly different vocabulary)
	 * 
	 * @param filebase
	 * @param parmode
	 */
	public CorpusResolver(String filebase, boolean parmode) {
		this.parmode = parmode;
		this.filebase = filebase;
		for (int i = 0; i < EXTENSIONS.length; i++) {
			String base = filebase;
			// read alternative vocabulary for paragraph mode
			if (parmode && EXTENSIONS[i].equals("vocab")) {
				base += ".par";
			}
			File f = new File(base + "." + EXTENSIONS[i]);
			if (f.exists()) {
				data[i] = load(f);
			}
		}
	}

	/**
	 * check whether labels exist
	 * 
	 * @param i
	 * @return 0 for no label values, 1 for yes, 2 for loaded, -1 for illegal
	 *         index
	 */
	public int hasValues(int i) {
		if (i >= EXTENSIONS.length || i < 0) {
			return -1;
		}
		// in the current impl, labels are pre-fetched
		return (data[i] != null ? 2 : 0);
	}

	/**
	 * load from file removing every information after a = sign in each line
	 * 
	 * @param f
	 * @return array of label strings
	 */
	private String[] load(File f) {
		String[] strings = null;
		try {
			ArrayList<String> a = new ArrayList<String>();
			BufferedReader br = new BufferedReader(new FileReader(f));
			String line = null;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				int ii = line.indexOf('=');
				if (ii > -1) {
					a.add(line.substring(0, ii).trim());
				} else {
					a.add(line.trim());
				}
			}
			br.close();
			strings = a.toArray(new String[] {});
			return strings;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return strings;
	}

	/**
	 * filter corpus with term subset with new indices.
	 * 
	 * @param old2new element (old index) contains new index
	 */
	public void filterTerms(int[] old2new) {
		HashMap<String, Integer> newids = new HashMap<String, Integer>();
		List<String> terms = new ArrayList<String>();
		// replace term ids.
		for (int i = 0; i < old2new.length; i++) {
			if (old2new[i] >= 0) {
				newids.put(getTerm(i), old2new[i]);
				terms.add(old2new[i], getTerm(i));
			}
		}
		data[LabelNumCorpus.LTERMS + 2] = (String[]) terms
				.toArray(new String[0]);
		termids = newids;
	}

	/**
	 * write the term set to the file
	 * 
	 * @param file (full file name, no .vocab appended)
	 * @throws IOException
	 */
	public void writeTerms(String file) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		for (String term : data[LabelNumCorpus.LTERMS + 2]) {
			bw.append(term).append('\n');
		}
		bw.close();
	}

	/**
	 * resolve the numeric term id
	 * 
	 * @param t
	 * @return
	 */
	public String getTerm(int t) {
		if (data[1] != null) {
			return data[1][t];
		} else {
			return null;
		}
	}

	/**
	 * find id for string term or return -1 if unknown
	 * 
	 * @param term
	 * @return
	 */
	public int getTermId(String term) {
		if (termids == null) {
			termids = new HashMap<String, Integer>();
			for (int i = 0; i < data[1].length; i++) {
				termids.put(data[1][i], i);
			}
		}
		Integer id = termids.get(term);
		if (id == null) {
			id = -1;
		}
		return id;
	}

	/**
	 * resolve the numeric label id
	 * 
	 * @param i
	 * @return
	 */
	public String getLabel(int i) {
		if (data[3] != null) {
			return data[3][i];
		} else {
			return null;
		}
	}

	/**
	 * resolve the numeric author id
	 * 
	 * @param i
	 * @return
	 */
	public String getAuthor(int i) {
		if (data[2] != null) {
			return data[2][i];
		} else {
			return null;
		}
	}

	/**
	 * resolve the numeric term id
	 * 
	 * @param i
	 * @return
	 */
	public String getDoc(int i) {
		if (data[0] != null) {
			return data[0][i];
		} else {
			return null;
		}
	}

	/**
	 * resolve the numeric term id
	 * 
	 * @param i
	 * @return
	 */
	public String getDocName(int i) {
		if (data[5] != null) {
			return data[5][i];
		} else {
			return null;
		}
	}

	/**
	 * resolve the numeric term id
	 * 
	 * @param i
	 * @return
	 */
	public String getDocKey(int i) {
		if (data[6] != null) {
			return data[6][i];
		} else {
			return null;
		}
	}

	/**
	 * filters the documents according to the new index
	 * 
	 * @param index
	 */
	public void filterDocuments(int[] old2new) {
		throw new Error();
	}

	/**
	 * resolve the numeric volume id
	 * 
	 * @param i
	 * @return
	 */
	public String getVol(int i) {
		if (data[4] != null) {
			return data[4][i];
		} else {
			return null;
		}
	}

	public String getLabel(int type, int id) {
		if (type == LabelNumCorpus.LTERMS) {
			return getTerm(id);
		} else if (type == LabelNumCorpus.LAUTHORS) {
			return getAuthor(id);
		} else if (type == LabelNumCorpus.LCATEGORIES) {
			return getLabel(id);
		} else if (type == LabelNumCorpus.LVOLS) {
			return getVol(id);
		} else if (type == LabelNumCorpus.LREFERENCES) {
			return getDocKey(id);
		} else if (type == LabelNumCorpus.LDOCS) {
			return getDoc(id);
		}
		return null;
	}

	public int getId(int type, String label) {
		if (type == LabelNumCorpus.LTERMS) {
			return getTermId(label);
		} else if (type == LabelNumCorpus.LAUTHORS) {
			return Arrays.asList(data[type]).indexOf(label);
		}
		return -1;
	}
}
