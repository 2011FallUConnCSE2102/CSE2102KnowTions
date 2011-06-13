knowceans-tools

@date 2005-03-28 update 2010-08-23
@author Gregor Heinrich questions / remarks email gregor :: arbylon . net

License: These files are licensed under the GPL version 3. See license.txt for details.

This package provides a couple of helper classes for diverse programming problems:

 - org.knowceans.util.* 
   miscellaneous classes, such as 
   
   o ExternalProcess : an interface to external processes that captures output in a 
     buffered streamreader that runs in its own thread
     
   o Conf : a properties file reader that extends the properties file format by a simple 
     way to define variables and simplify complex configurations that would require tedious
     editing.
   
   o WebLogAnalyzer : a "hack" to analyse Tomcat logs whose target log file format can be 
     easily be adjusted via regular expressions.
     
   o TableList : a list implementation that supports parallel sorting and lookup like a 
     table in a relational database. KeyTableList extends this functionality by HashMap
     keys to lookup list entries.
     
   o Densities, Samplers, Cokus, CokusRandom : random number generation and likelihood esti-
     mation for various probability distributions.
     
   o PatternString : simplifies working with regex patterns, with a syntax inspired by Perl.
   
   o SetVector, SetArrayList : implement both the set and list interface (i.e., sorted sets)
   
   o Vectors, ArrayUtils : many helper functions for arrays, specifically of primitive 
     elements.
     
   o StopWatch, Which : runtime introspection for memory, running time and stack traces.
     
   o ...
   
   
 - org.knowceans.map.*
   Java 1.5 implementation of 1:1, 1:n, and m:n relations in top of type-safe HashMaps, 
   allowing reverse lookup and partly access to keys via wildcards and regular expressions
   
   o BijectiveHashMap : 1:1 mapping with fast search.
   
   o BijectiveTreeMap : 1:1 mapping with sorting in both domain (keys) and co-domain (values).
   
   o HashMultiMap : m:n mapping, i.e., a HashMap that allows a key map to different values,
     however without reverse lookup (from values to keys). Allows wildcard and regex lookup.
   
   o InvertibleHashMap : 1:n mapping with reverse lookup (cf. HashMap, which is a 1:n mapping
     without reverse lookup.)
     
   o InvertibleHashMultiMap : m:n mapping with reverse lookup. Allows wildcard and regex lookup.
   
   + The implementation fully uses Java-1.5's generics features. (But beware: testing has not 
     been extensive after the generics have been introduced to the *.map package)

   + Different interfaces allow to develop other implementations, e.g., TreeMap-based or such that
     leave the implementation of Map generic
     
   org.knowceans.corpus.*
   Java-based handling of text corpora in a format similar to svmlight, which was originally used 
   by lda-c by Dave Blei and consequently my Java port, lda-j. added this package 2010-08-23 as part
   of codebase maintenance, code was originally written between 2005 and 2009. An example of a corpus
   is provided in the directory corpus-example. This is the widely used NIPS test collection that I
   have manually added year/volume information as well as labels as far as they could be possibly 
   be inferred from the NIPS conference page. 
   
   o Document: represents a text document

   o DisjointDocTerms: partitions the corpus in sets of disjoint documents and terms (to sync
     parallel threads)
   
   o NumCorpus: a corpus reader, accessing a *.corpus file with term vectors in svmlight format, 
     a *.docs file with document names and a *.vocab file with term indices (one line per term). 
     Allows splitting for cross validation.
     
   o LabelNumCorpus: a corpus reader, accessing the files of NumCorpus, plus optional labels,
     for which different kinds can be defined: authors, classes, years, references. Allows multimodal
     models to be fed with data. There are two types of label files: data (one line per entry, matching
     the document indices) and key (one line per entry, matching label indices). This is similar to
     the approach to the vocabulary
     
   o I*Corpus: interfaces to define the behaviour of the corpus implementations, including support 
     for labels, cross-validation, term lookup with the CorpusResolver.
     
   o CorpusResolver: class used to actually resolve indices of items into names.
   
 - org.knowceans.sandbox.* 
   some classes for quick playing around
   
The classes have been developed with the Eclipse platform, for which the .project and .classpath 
files are provided. Further, an ant build.xml exists.

