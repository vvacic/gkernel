/* ProteinContactGraph class
 *
 * Vladimir Vacic
 * Algorithms and Computational Biology Lab
 * Department of Computer Science and Engineering
 * University of California, Riverside
 * Riverside, CA 92521, USA
 *
 * May-14-2008
 */

#ifndef __PROTEIN_CONTACT_GRAPH_H__
#define __PROTEIN_CONTACT_GRAPH_H__

#include "pdb.h"
#include <climits>
#include <iostream>
#include <map>
#include <vector>
using namespace std;


class ProteinContactGraph  {
public:
    ProteinContactGraph()  {}
    ProteinContactGraph(Chain *);
    ProteinContactGraph(Chain *, Chain::ConnectionMethod, float);
    ~ProteinContactGraph()  {}


    /** Read a protein contract graph file and create a graph. */
    static ProteinContactGraph read_compact(const char*);

    /** Breadth-first assignment of distances from the center. */ 
    vector<unsigned> breadth_first_sort() const;

    /** Generates a protein structure graph based on the method and 
     *  distance threshold. Method can be one of the several listed
     *  in the ConnectionMethod enumeration. */
    void connect(Chain::ConnectionMethod, float);

    /** Includes only the residues with C-alphas within the sphere of
     *  the given radius, centered at the C-alpha of the residue of
     *  interest.
     *
     *  In the resulting ProteinContactGraph, center of the sphere
     *  becomes the first residue. The remaining residues which fall
     *  within the sphere are added in the order of increasing shortest
     *  path from the center, i.e. all residues 1 hop waya are added
     *  first, then all residues 2 hops away, and so forth. */
    ProteinContactGraph filter_sphere(const string&, float);

    /** Describe all residues located within a sphere of given radius from
     *  an input parameter residue. */
    vector<int> describe_shell(const string&, float, float);

    /** Generates a GraphViz file. */
    void write_dot_file(const char *);

    /** Generates a GraphViz file. */
    void write_dot_file(const char *, const vector<unsigned>&, unsigned);

    /** Print the edges of the protein contact graph. Typical lines
     *  of the resulting file look something like:
     *
     *      ARG105	PHE84	3.80736
     *      ARG105	GLU85	3.80736
     *      ARG105	VAL104	6.80172	backbone
     *
     * The columns are tab-separated, and the edges are from the lower
     * triangular adjacency matrix, i.e. the second residue has a lower
     * index than the first. 
     * 
     * NOTE: Backbone connections are not explicitly saved during
     * filtering, because filtering orders residues according to the
     * shortest path distance from the center, and not from N to C
     * terminus. 
     * 
     * Likewise, a filtered graph does not have a real Chain data
     * member, that is, the Chain does not contain any residues.
     */
    void write_edge_list(const char *);

    /** Print the protein contact graph in compact form. */
    void write_compact(const char *);
    
    Chain *ch;
    vector<string> nodes;
    vector<vector<unsigned> > edges;
    vector<vector<float> > distances;

private:
    /** This method adds both i-j and j-i edges in the adjacency list, so
     *  it is important not to invoke it both for i and j in order not to
     *  have duplicate edges in the graph. */
    void add_edge(unsigned, unsigned, float); 
};

#endif
