/* GraphletKernel class. Contains methods for reading graphs,
 * efficiently computing a kernel matrix and writing the output in
 * several formats.
 *
 * Vladimir Vacic
 * Algorithms and Computational Biology Lab
 * Department of Computer Science and Engineering
 * University of California, Riverside
 * Riverside, CA 92521, USA
 *
 * May-7-2008
 */

#ifndef __GKERNEL_H__
#define __GKERNEL_H__

#define GRAPHLET_SIZE 4
#define LOG_ALPHA_SIZE 5
#define ZERO_CHAR '@'

//JLM: If 0 exclude pivot from label generation/representation of
//graphlets and only consider non-pivot label matches(less stringent),
//otherwise include pivot in graphlet label representation
//(more stringent).
#define PIVOT 0

#include "aminoacid.h"
#include "protgraph.h"
#include <fstream>
#include <map>
#include <vector>
using namespace std;

typedef unsigned Key;


class GraphletKernel  {
public:
    GraphletKernel()  {}
    ~GraphletKernel()  {}

    inline static void set_verbose()  {
        GraphletKernel::VERBOSE = true;
    }

    inline static void set_debug()  {
        GraphletKernel::DEBUG = true;
    }

    inline static void set_normalize()  {
        GraphletKernel::NORMALIZE = true;
    } 

    inline static void set_reduction(AminoAcid::ReductionScheme scheme)  {
        GraphletKernel::REDUCTION = scheme;
    }
     
    /** Reads a list of protein contact graph files. */ 
    void read_graphs(string, const vector<string> &);

    /** Computes a graphlet kernel matrix over graphs. */
    void compute_matrix();

    /** Writes binary kernel matrix. */
    void write_binary_kernel(const char*);
    
    /** Writes kernel matrix in Matlab format. */
    void write_sparse_matlab(const char*);

    /** Writes kernel matrix in SVM^Light format. */
    void write_sparse_svml(const char *);

    /** Writes lower triangular kernel matrix. */
    void write_text_matrix(const char*);

    /** Set binary classification labels (positives vs. negatives). */
    inline void set_labels(const vector<int> &lab)  {
        labels = lab;
    }

    /** Writes binary classification labels (positives vs. negatives). */
    void write_labels(const char *);


private:
    static string print_key(Key);

    static Key make_key(char, char, char, unsigned, char);

    static void increment_hash(map<Key,float>&, const Key&);

    /** Returns the counts of labeled graphlets. */
    static map<Key,float> count(ProteinContactGraph&);

    /** Normalizes the kernel matrix using the method for normalizing the 
     *  spectral kernel matrix. */
    static void normalize_spectral(map<Key,float>&);

    /** Computes graphlet distance between two graphs.*/
    static float distance_hash_join(map<Key,float>, map<Key,float>);


    static bool VERBOSE;
    static bool DEBUG;
    static bool NORMALIZE;
    static AminoAcid::ReductionScheme REDUCTION;

    vector<int> labels;
    vector<ProteinContactGraph> graphs;
    vector<map<Key,float> > hashes;
    vector<vector<float> >  kernel;
};

#endif
