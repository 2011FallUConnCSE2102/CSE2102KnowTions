This is the NIPS0-12 test collection

@author gregor :: arbylon . net
@date 2010-08-23 (release)
  
Contents and line format: 
 - nips.authors        -- document author ids, 1 line per document (line 1 = document 0), 
                          space-separated
 - nips.authors.key    -- author names, 1 line per author (line 1 = author 0)
 - nips.corpus         -- term vectors, 1 line per document, svm-light line format: 
                          nterms (termid:termfreq)+
 - nips.docs           -- document titles, 1 line per document (line 1 = document 0)
 - nips.labels         -- class labels, 1 line per document (line 1 = document 0)
 - nips.labels.extract -- class merging information (original to final labels)
 - nips.label.key      -- class label names, 1 line per label (line 1 = class 0)
 - nips.split          -- permutation of the document set, can be used for random splitting
 - nips.vocab          -- vocabulary index, 1 line per term (line 1 = term 0)
 - nips.vols           -- volume information, 1 line per document (line 1 = volume 0 
                          = year 1987, one volume per year continuously)
 - nips.refs           -- citation information, 1 line per document


Sources:
OCR'ed files: http://nips.djvuzone.org/
Original files: http://books.nips.cc/
 --> year/volume information
 --> labels (similar labels have been merged, nips.labels.extract)
Preprocessed MAT files: http://cs.nyu.edu/~roweis/data.html
 --> term vectors
 --> document titles
 --> vocabulary  
Text files: http://cs.nyu.edu/~roweis/data.html, http://ai.stanford.edu/~gal/
 --> citation entries (via searches for nips title and abbreviations and manual
     association with author/title/docnames)
