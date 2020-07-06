## pdb_mapper
For a given input protein, to map its motif onto to a PDB complex(s) with with known motif-domain interactions through 3DID & ELM

#### If the source is [3DID](https://3did.irbbarcelona.org/)
1. For a given pairt of Motif and its interacting-domain, retrieve all the PDB complexes from 3DID database along with their chain and start/end annotations of the pair
2. Superimpose in PyMOL the RegEx match (peptide) of the motif on the input protein to that of the motif on the motif-chain (coordinates of the latter known from 3DID as mentioned in step 1)
3. Color and annotate the chains on the complexes and save their PyMOL sessions as well as the new PDB coordinates (which include the peptide)

#### If the source is [ELM](http://elm.eu.org/)
1. For a given pair of motif and its interacting-domain, extract all the PDB complexes known through the Motif class' instances and have the interacting domain on any of the  (latter using [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/index.html) PDB_CHAIN_TO_PFAM_ID)
2. Using SIFTS PDB_CHAIN_TO_UNIPROT_ACC, further filter and select only those complexes that have the instance's motif present on chain other than the one having the interacting domain
3. Map the start/end positions of the RegEx match on the instance to the motif chain using BLASTp (retrieve the FASTA sequence of the instance using its UniProt accession and that of the chain using PDB-AA)
4. Find the start/end positions of the domain on the interacting-domain chain using HMMSCAN
5. Superimpose in PyMOL the RegEx match (peptide) of the motif on the input protein to that of the motif on the motf-chain
6. Color and annotate the chains on the complexes and save their PyMOL sessions as well as the new PDB coordinates (which include the peptide)

#### References
1. Roberto Mosca, Arnaud Ceol, Amelie Stein, Roger Olivella & Patrick Aloy, 3did: a catalogue of domain-based interactions of known three-dimensional structure Nucl. Acids Res. 2014, 42(D1):D374-D379. [link](https://sbnb.irbbarcelona.org/3did_2013)
2. Manjeet Kumar, Marc Gouw, Sushama Michael, Hugo Sámano-Sánchez, Rita Pancsa, Juliana Glavina, Athina Diakogianni, Jesús Alvarado Valverde, Dayana Bukirova, Jelena Čalyševa, Nicolas Palopoli, Norman E Davey, Lucía B Chemes, Toby J Gibson, ELM—the eukaryotic linear motif resource in 2020, Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D296–D306. [link](https://doi.org/10.1093/nar/gkz1030)
3. Jose M Dana, Aleksandras Gutmanas, Nidhi Tyagi, Guoying Qi, Claire O’Donovan, Maria Martin, Sameer Velankar, SIFTS: updated Structure Integration with Function, Taxonomy and Sequences resource allows 40-fold increase in coverage of structure-based annotations for proteins, Nucleic Acids Research, Volume 47, Issue D1, 08 January 2019, Pages D482–D489. [link](https://doi.org/10.1093/nar/gky1114)
