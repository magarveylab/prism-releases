## PRISM  

PRISM (PRediction Informatics for Secondary Metabolomes) is a software package for natural product genome mining. PRISM identifies biosynthetic gene clusters and generates libraries of hypothetical structures based on the identified genetic information. 

## Usage

PRISM is a Java 8 package, which can be deployed to a Tomcat 7 web server as a web application, or run from the command line. 

A GUI is available at http://magarveylab.ca/prism. 

A sample command line run might look like: 

```
$ java -jar prism.jar -a -p -f genome.fa -tt -sug -res -rib -w 10000
```

This runs PRISM with open reading frames predicted both by Prodigal and by identifying all possible coding sequences (Prodigal orfs are prioritized) on the file genome.fa, enabling searches for thiotemplated (nonribosomal peptide and polyketide), deoxysugar, resistance, and ribosomal natural product biosynthetic enzymes, with a cluster window of 10,000 bp. 

To see a detailed description of all command line options, do:

```
$ java -jar prism.jar -h 
``` 

## Dependencies

Core dependencies include HMMER (3.1) and BLAST (2.2.25+). Optional but highly recommended dependencies include BioPerl (1.006924) for GenBank file input, Prodigal (2.6.3) for prokaryotic open reading frame prediction, and FIMO (for ribosomally produced and post-translationally modified peptide precursor cleavage). 

## Citing PRISM

Please cite: 

> Skinnider MA, Dejong CA, Rees PN, Johnston CW, Li H, Webster ALH, Wyatt MA, Magarvey NA (2015). Genomes to natural products PRediction Informatics for Secondary Metabolomes (PRISM). _Nucleic Acids Research_, **16**, 9645-62. doi: [10.1093/nar/gkv1012](http://dx.doi.org/10.1093/nar/gkv1012)