## pdb_mapper
To map a given pair of motif-domain or domain-domain onto to the best PDB complex

#### If the source is [3DID](https://3did.irbbarcelona.org/)
1. Make a BLAST database of all the PDB complexes for the given motif-domain
2. BLASTp the query_fasta (which has the motif) against the database
3. Select the query_chain that has highest identity and includes the given  motif
4. Retrieve the interacting_chain of the best query_chain from 3DID
5. Map them onto their PDB complex, highlighting the motif and domain

#### If the source is [ELM](http://elm.eu.org/)
1. Extract all the PDB complexes known for the given motif through its instances
2. Using SIFTS pdb to pfam, filter out only those complexes that had the given interacting domain
3. In the filtered PDB complexes' list, make a BLAST database of the all the chains other than the ones having the domain, and call them as query_chains
4. Besides, run a regex of the given motif on the query_chains
5. BLASTp the query_fasta with query_chains
6. Finally, select only those complexes that have at least one chain in the query_chains' list that has the motif present, and at least one another (not the same as the query_chain) that has the given domain
7. Use BLASTp to remove the PDB complexes that do not have the regex of its chain overlapping with that of the query's start and end of the motif
8. Use HMMSCAN to find the start and end of the domain on a chain
9. Map the motif and domain onto to their PDB complex
