This repository contains beta version of the graphlet kernel code and 
supplementary information for:

> Vacic V, Iakoucheva LM, Lonardi S, Radivojac P. (2010) "Graphlet 
kernels for prediction of functional residues in protein structures." 
*Journal of Computational Biology*. 17(1):55-72. 
[PMID:20078397](https://pubmed.ncbi.nlm.nih.gov/20078397/).


We introduced a novel graph-based kernel method for annotating
functional residues in protein structures. A structure is first modeled
as a protein contact graph, where nodes correspond to residues and
edges connect spatially neighboring residues. Each vertex in the graph
is then represented as a vector of counts of labeled non-isomorphic
subgraphs (graphlets), centered on the vertex of interest. A similarity
measure between two vertices is expressed as the inner product of their
respective count vectors and is used in a supervised learning framework
to classify protein residues. The key benefit of the graphlet
representation is its ability to capture neighborhood similarities in
protein structures via enumerating the patterns of local connectivity
in the corresponding labeled graphs.

A detailed description of the method and discussion of parameters can be 
found in (Vacic *et al.*, 2010).  


## Contents:

* supplement - Supplementary information for (Vacic *et al.*, 2010)
* src - graphlet kernel source code


## Supplementary information

[Supplementary Table S1](supplement/Table_S1_Phos_PDB.xls): A
nonredundant subset of phosphorylated sites mined from PDB and results
of the search for the structures of proteins that can be found both in
the phosphorylated and the unphosphorylated states.

[Supplementary Table S2](supplement/Table_S2.pdf): Summary of the set of
sequences with experimentally annotated phosphorylation sites. NR:
non-redundant, here defined as having less than 40% sequence identity
in the 25 residue-long fragment centered at S, T or Y.

[Supplementary Table S3](supplement/Table_S3.pdf): Performance of
method/parameter combinations on the CSA dataset. Within each group of
predictors, the one with the highest AUC was awarded a point. In the
case of ties, the point would be evenly split between the methods which
performed equally well.

[Supplementary Table S4](supplement/Table_S4.pdf): Performance of
method/parameter combinations on the PHOS dataset. Within each group of
predictors, the one with the highest AUC was awarded a point.

[Supplementary Table S5](supplement/Table_S5.pdf): Jaccard coefficients
between sets of edges in residue interaction networks obtained using
different methods. Jaccard similarity coefficient for two sets A and B
is defined as J(A,B) = |A \ B|/|A [ B|. Values over 70% are in bold
face and marked with an asterisk.

[Supplementary Figure S1](supplement/Figure_S1.pdf): Schematic
representation of the BLOSUM50 matrix-based amino acid alphabet
reduction.

[Supplementary Figures S2-S24](supplement/Figures_S2_S24.pdf): 
Performance of all predictors on the CSA and PHOS datasets. 


## Code:

### Compilation

To compile the graphlet kernel code, go to the `src` directory and type
`make` on the command prompt. This will generate program files
`parse_pdb` and `make_kernel`. 

NOTE: Global variable PIVOT in `gkernel.h` controls inclusion of the pivot
vertex in the counts of graphlets. If PIVOT is set to 0 then only
non-pivot label mismatches are counted. If PIVOT is set to 1 then label
mismatches will also include the pivot node.


### Program options

```
Usage: parse_pdb [options] PDB_FILE
Options:

  -h           Displays this message.

  -c CHAIN     PDB chain.
               Defaults to the first chain found in the PDB file.
  -m METHOD    Method used for determining interacting residues.
               Can be C_ALPHA, C_BETA, ALL_ATOMS, ALL_VDW_RADIUS.
               Defaults to C_ALPHA.
  -d DISTANCE  Threshold distance in Angstroms.
               Defaults to 6A.

  -r RESIDUE   Only the subgraph centered at this residue.
               Residues should be specified by their residue sequence
               number and an insertion code, if one exists (resSeq and
               iCode, columns 23-26 and 27 in the ATOM/HETATM lines).
  -s RADIUS    Includes only the residues with C alphas within the 
               sphere of the given RADIUS, centered at the C alpha
               of the -r residue.

  -E or -C     Output file format, edge list or compact.
               Defaults to edge list.
  -o OUTFILE   Output file for the results.
```

```
Usage: make_kernel -p FILE -n FILE -a PATH -[k|m|s] OUTPUT [...]
Options:

  -h         Displays this message.

  -p FILE    List of positives.
  -n FILE    List of negatives.
  -a PATH    Path to subgraph files.

  -r REDUCT  Alphabet reduction scheme. Can be NO_REDUCTION,
             UNLABELED, BLOSUM_2, BLOSUM_3, BLOSUM_4, BLOSUM_5,
             BLOSUM_6, BLOSUM_8, BLOSUM_10 or BLOSUM_15.
             Defaults to NO_REDUCTION.

  -N         Normalize the kernel matrix.
             Defaults to false.

  -k KERNEL  Output file for the binary kernel matrix.
   or
  -m MATLAB  Output file for the sparse attribute matrix (Matlab).
   or
  -s SVML    Output file for the sparse attribute matrix (SVM^Light).
   or
  -t TEXT    Output file for the plain text kernel matrix.
             Defaults to binary kernel matrix.

  -l LABELS  Output file for the example labels.

  -V         Verbose (prints progress messages).
             Defaults to false.

  -D         Debug (prints debug messages).
             Defaults to false.
```


### Usage example 

The data subdirectory contains a single PDB file (`3LCK.pdb`), and two
files with positives (functional) and negative (not-functional) residues,
`LYS_3LCK.positives` and `LYS_3LCK.negatives`. As an illustration, all 
lysines inside 3LCK were arbitrarily split into positives and negatives.

PDB files can be downloaded from the Protein Data Bank
(https://www.rcsb.org).

The "positives" and "negatives" files are tab-separated, with three
fields, PDB_CODE, CHAIN and RESIDUE.  

Preparing a kernel matrix from PDB files is a two step process:

(1) `parse_pdb` is used to parse PDB files into protein contact graphs.
The user can specify the connection method (C_ALPHA, B_BETA, ALL_ATOMS,
etc.), distance thresholds, etc. See (Vacic *et al.*, 2010) for
details.

(2) `make_kernel` is used to generate the graphlet count representation
based on the graph files. There are several output options, out of
which SVM^Light format is probably the easiest to use, because it can
be readily read by SVM^Light (see http://svmlight.joachims.org). In a
nutshell, this is a space-separated file with the first field equal to
1(for positives) or -1 (for negatives), and all other non-zero entries
are feature_key:feature_value pairs.  

A sample `example.sh` shell script is provided as an illustrative
example of the process. 


### Comment regarding large datasets

Using the SVM^Light format directly on very large datasets is not
efficient, because the dot-product between feature vectors has to be
unnecessarily recomputed over and over again. There is an option to
precompute the kernel matrix and save it as a binary file, however,
using it with SVM^Light is slightly more complicated, because the
reader for this file format has to be custom coded and according to the
SVM^Light copyright agreement we are not allowed to distribute modified
SVM^Light code. 

If efficiency is an issue and you decide to tweak SVM^Light, in
`svm_learn_main.c`, function `main()` you may add something in the
likes of:

```
CFLOAT *distmat;
FILE *f;
unsigned int numgraphs, numcells;

if (kernel_parm.kernel_type == 4) {
    f = fopen(kernel_parm.custom, "rb");

    fread(&numgraphs, sizeof(unsigned int), 1, f);
    numcells = (numgraphs+1) * numgraphs / 2;

    distmat = (CFLOAT *) my_malloc(numcells * sizeof(CFLOAT));
    fread(distmat, sizeof(CFLOAT), numcells, f);

    fclose(f);

    kernel_parm.distmat = distmat;
}
```

This will load the pre-computed dot-product matrix.

And in `svm_common.c`, function 
`CFLOAT single_kernel(KERNEL_PARM *kernel_parm,  SVECTOR *a, SVECTOR *b)`:

```
case 4:
    i = atoi(a->userdefined);
    j = atoi(b->userdefined);

    if (i >= j)
        return kernel_parm->distmat[(i+1)*i/2 + j];
    else
        return kernel_parm->distmat[(j+1)*j/2 + i];
```

For every pair of vectors, this will read their pre-computed dot-product from
the matrix.


### Summary of changes

1.1 (Jose Lugo-Martinez, School of Informatics and Computing, Indiana University, Bloomington)

* Added pivot as a paramater of make_key() function, which optionally includes pivot in the graphlets labeling representation. 
* Added additional graphlets counts for case 0123 to include all paths that lead to this case.
* Fixed bugs in counting orbits 12 and 13 and a small typo in case 0112.
* Fixed a bug in label generation in make_key() function for cases 8 and 15.

1.01

* Added basic index checking into the graph reading function.
