/*
 * (C) Copyright 2005, Gregor Heinrich (gregor :: arbylon : net) (This file is
 * part of the lda-j (org.knowceans.lda.*) experimental software package.)
 */
/*
 * lda-j is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 */
/*
 * lda-j is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 */
/*
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*
 * Created on Dec 3, 2004
 */
package org.knowceans.corpus;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import org.knowceans.util.Vectors;

/**
 * Represents a corpus of documents, using numerical data only.
 * <p>
 * 
 * @author heinrich
 */
public class LabelNumCorpus extends NumCorpus implements ILabelCorpus {

	/**
	 * the extensions for the label type constants in ILabelCorpus.L*
	 */
	public static final String[] EXTENSIONS = { ".authors", ".labels", ".vols",
	/* ".refs", */".cite", ".tags", ".years" };

	/**
	 * array of labels. Elements are filled as soon as readlabels is called.
	 */
	int[][][] labels;
	/**
	 * total count of labels
	 */
	int[] labelsW;
	/**
	 * total range of labels
	 */
	int[] labelsV;

	String dataFilebase = null;

	private CorpusResolver resolver;

	/**
     * 
     */
	public LabelNumCorpus() {
		super();
		init();
	}

	/**
	 * @param dataFilebase (filename without extension)
	 */
	public LabelNumCorpus(String dataFilebase) {
		super(dataFilebase + ".corpus");
		this.dataFilebase = dataFilebase;
		init();
	}

	/**
	 * @param dataFilebase (filename without extension)
	 * @param parmode if true read paragraph corpus
	 */
	public LabelNumCorpus(String dataFilebase, boolean parmode) {
		super(dataFilebase + (parmode ? ".par" : "") + ".corpus");
		this.dataFilebase = dataFilebase;
		init();
	}

	/**
	 * @param dataFilebase (filename without extension)
	 * @param readlimit number of docs to reduce corpus when reading (-1 =
	 *        unlimited)
	 * @param parmode if true read paragraph corpus
	 */
	public LabelNumCorpus(String dataFilebase, int readlimit, boolean parmode) {
		super(dataFilebase + (parmode ? ".par" : "") + ".corpus", readlimit);
		this.dataFilebase = dataFilebase;
		init();
	}

	/**
	 * create label corpus from standard one
	 * 
	 * @param corp
	 */
	public LabelNumCorpus(NumCorpus corp) {
		this.docs = corp.docs;
		this.numDocs = corp.numDocs;
		this.numTerms = corp.numTerms;
		this.numWords = corp.numWords;
		init();
	}

	protected void init() {
		labels = new int[EXTENSIONS.length][][];
		labelsW = new int[EXTENSIONS.length];
		labelsV = new int[EXTENSIONS.length];
	}

	/**
	 * checks whether the corpus has labels
	 * 
	 * @param kind according to label constants ILabelCorpus.L*
	 * @return 0 for no label values, 1 for yes, 2 for loaded, -1 for illegal
	 */
	public int hasLabels(int kind) {
		if (kind >= labels.length || kind < -2) {
			return -1;
		}
		if (kind < 0) {
			// we have docs and terms loaded
			if (docs != null)
				return 2;
		} else {
			// any metadata
			if (labels[kind] != null) {
				return 2;
			}
			File f = new File(this.dataFilebase + EXTENSIONS[kind]);
			if (f.exists()) {
				return 1;
			}
		}
		// not loaded
		return 0;
	}

	/**
	 * loads and returns the document labels of given kind
	 */
	// @Override
	public int[][] getDocLabels(int kind) {
		if (labels[kind] == null)
			readLabels(kind);
		return labels[kind];
	}

	/**
	 * get the labels for one document
	 * 
	 * @param m
	 * @param kind
	 * @return
	 */
	public int[] getDocLabels(int kind, int m) {
		if (labels[kind] == null)
			readLabels(kind);
		return labels[kind][m];
	}

	/**
	 * return the maximum number of labels in any document
	 * 
	 * @param kind
	 * @return
	 */
	public int getLabelsMaxN(int kind) {
		int max = 0;
		for (int m = 0; m < labels[kind].length; m++) {
			max = max < labels[kind][m].length ? labels[kind][m].length : max;
		}
		return max;
	}

	// @Override
	public int getLabelsW(int kind) {
		return labelsW[kind];
	}

	// @Override
	public int getLabelsV(int kind) {
		return labelsV[kind];
	}

	/**
	 * read a label file with one line per document and associated labels
	 * 
	 * @param kind
	 * @return
	 */
	private void readLabels(int kind) {
		ArrayList<int[]> data = new ArrayList<int[]>();
		int W = 0;
		int V = 0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(dataFilebase
					+ EXTENSIONS[kind]));
			String line;
			int j = 0;
			while ((line = br.readLine()) != null) {
				// remove additional info
				int c = line.indexOf(" : ");
				if (c > -1) {
					line = line.substring(0, c);
				}
				line = line.trim();
				if (line.length() == 0) {
					data.add(new int[0]);
					continue;
				}
				String[] parts = line.split(" ");
				int[] a = new int[parts.length];
				for (int i = 0; i < parts.length; i++) {
					a[i] = Integer.parseInt(parts[i].trim());
					if (a[i] >= V) {
						V = a[i] + 1;
					}
				}
				W += a.length;
				data.add(a);
				j++;
			}
			br.close();
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		labels[kind] = data.toArray(new int[0][0]);
		labelsW[kind] = W;
		labelsV[kind] = V;
	}

	@Override
	public void split(int order, int split, Random rand) {
		// get plain num corpora
		super.split(order, split, rand);

		// now also split labels
		int Mtest = splitstarts[split + 1] - splitstarts[split];
		labelsW = new int[EXTENSIONS.length];

		int[][][] trainLabels = new int[EXTENSIONS.length][numDocs - Mtest][];
		int[][][] testLabels = new int[EXTENSIONS.length][Mtest][];
		int[] trainLabelsW = new int[EXTENSIONS.length];
		int[] testLabelsW = new int[EXTENSIONS.length];

		int mstart = splitstarts[split];
		// for each label type
		for (int type = 0; type < EXTENSIONS.length; type++) {
			if (labels[type] == null) {
				continue;
			}
			int mtrain = 0;
			// before test split
			for (int m = 0; m < splitstarts[split]; m++) {
				trainLabels[type][mtrain] = labels[type][splitperm[m]];
				mtrain++;
			}
			// after test split
			for (int m = splitstarts[split + 1]; m < numDocs; m++) {
				trainLabels[type][mtrain] = labels[type][splitperm[m]];
				mtrain++;
			}
			// test split
			for (int m = 0; m < Mtest; m++) {
				testLabels[type][m] = labels[type][splitperm[m + mstart]];
				testLabelsW[type] += testLabels[type][m].length;
			}
			trainLabelsW[type] = labelsW[type] - testLabelsW[type];
		}
		// construct subcorpora
		trainCorpus = new LabelNumCorpus((NumCorpus) getTrainCorpus());
		testCorpus = new LabelNumCorpus((NumCorpus) getTestCorpus());
		LabelNumCorpus train = (LabelNumCorpus) trainCorpus;
		train.labels = trainLabels;
		train.labelsV = labelsV;
		train.labelsW = trainLabelsW;
		LabelNumCorpus test = (LabelNumCorpus) testCorpus;
		test.labels = testLabels;
		test.labelsV = labelsV;
		test.labelsW = testLabelsW;
	}

	@Override
	public void write(String pathbase) throws IOException {
		super.write(pathbase);
		// TODO: write authors and labels
	}

	/**
	 * test corpus reading and splitting
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		LabelNumCorpus nc = new LabelNumCorpus("berry95/berry95");
		nc.getDocLabels(LAUTHORS);
		nc.split(10, 0, new Random());

		System.out.println("train");
		LabelNumCorpus ncc = (LabelNumCorpus) nc.getTrainCorpus();
		System.out.println(ncc);
		int[][] x = ncc.getDocWords(new Random());
		System.out.println(Vectors.print(x));
		System.out.println("labels");
		int[][] a = ncc.getDocLabels(LAUTHORS);
		System.out.println(Vectors.print(a));

		System.out.println("test");
		ncc = (LabelNumCorpus) nc.getTestCorpus();
		System.out.println(ncc);
		x = ncc.getDocWords(new Random());
		System.out.println(Vectors.print(x));
		System.out.println("labels");
		a = ncc.getDocLabels(LAUTHORS);
		System.out.println(Vectors.print(a));

		System.out.println("document mapping");
		System.out.println(Vectors.print(nc.getOrigDocIds()));
	}

	/**
	 * get a resolver that acts on this corpus
	 * 
	 * @return
	 */
	public CorpusResolver getResolver() {
		if (resolver == null) {
			resolver = new CorpusResolver(dataFilebase);
		}
		return resolver;
	}
}
