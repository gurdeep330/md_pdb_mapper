## md_pdb_mapper
To map a given pair of motif-domain or domain-domain onto to the best PDB complex

#### If the source is [3DID](https://3did.irbbarcelona.org/)
1. Make a BLAST database of all the PDB complexes for the given motif-domain
2. BLASTp the query_fasta (which has the motif) against the database
3. Select the query_chain that has highest identity and includes the given  motif
4. Retrieve the interacting_chain of the best query_chain from 3DID
5. Map them onto their PDB complex, highlighting the motif and domain
