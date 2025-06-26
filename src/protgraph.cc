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

#include "aminoacid.h"
#include "protgraph.h"
#include <fstream>
#include <set>
#include <queue>


ProteinContactGraph::ProteinContactGraph(Chain *c)  {
    ch = c; 

    for (unsigned i=0; i<ch->residues.size(); i++)
        nodes.push_back(ch->residues[i].residue_id());

    edges.resize(nodes.size());
    distances.resize(nodes.size());
}


ProteinContactGraph::ProteinContactGraph(Chain *c, Chain::ConnectionMethod method, float threshold)  {
    ch = c; 

    for (unsigned i=0; i<ch->residues.size(); i++)
        nodes.push_back(ch->residues[i].residue_id());

    edges.resize(nodes.size());
    distances.resize(nodes.size());

    connect(method, threshold);
}


ProteinContactGraph ProteinContactGraph::read_compact(const char *file)  {
    ifstream rin(file, ios::in);

    if (rin.fail())  {
        cerr << "ERROR: Protein contact graph file " << file << " could not be opened." << endl;
        exit(1);
    }

    ProteinContactGraph g;

    string line;
    if (getline(rin, line))
        g.nodes = split(line, "");

    g.edges.resize(g.nodes.size());

    // Parse PCG file 
    while(getline(rin, line))  {
        vector<string> tokens = split(line, "\t");
        unsigned i = to_i(tokens[0]);

        for (unsigned j=1; j < tokens.size(); j++)
            g.edges[i].push_back( to_i(tokens[j]) ); 
    }

    return g;
}


vector<unsigned> ProteinContactGraph::breadth_first_sort() const  {
    vector<unsigned> dist(nodes.size(), UINT_MAX);    
    dist[0] = 0;
    
    queue<unsigned> Q;
    Q.push(0);
    
    while (!Q.empty())  {
        unsigned i = Q.front();
        
        for (unsigned j=0; j<edges[i].size(); j++)  {
            unsigned k = edges[i][j];
            
            if (UINT_MAX==dist[k])  {
                dist[k] = dist[i] + 1;
                //if (dist[k] < bound-1)
                Q.push(k);
            }
        }
        Q.pop();
    }

    return dist;
}


void ProteinContactGraph::connect(Chain::ConnectionMethod method, float threshold)  {
    float threshold_sq = threshold * threshold;
    Atom  *atom_t, *atom_u;

    for (unsigned i=1; i < ch->residues.size(); i++)  {
        atom_t = NULL;

        // For the C_ALPHA and C_BETA methods, find the pivot atoms only once. 
        if (Chain::C_ALPHA == method)  {
            atom_t = ch->residues[i].get_atom("CA");

            if (NULL == atom_t)  continue;  // It happens. 
        }
        else if (Chain::C_BETA == method)  {
            if ("GLY" == ch->residues[i].res_name)
                atom_t = ch->residues[i].get_atom("CA");
            else
                atom_t = ch->residues[i].get_atom("CB");

            if (NULL == atom_t)  continue;
        }

        float min_sq, min_dist;
		    
        for (unsigned j=0; j<i; j++)  {
            switch (method)  {
            case Chain::C_ALPHA:
                atom_u = ch->residues[j].get_atom("CA");

                if (NULL == atom_u)  continue;

                if (atom_t->dist_sq(*atom_u) <= threshold_sq)
                    add_edge(i, j, atom_t->dist(*atom_u));
                break;

            case Chain::C_BETA:
                if ("GLY" == ch->residues[j].res_name)
                    atom_u = ch->residues[j].get_atom("CA");
                else
                    atom_u = ch->residues[j].get_atom("CB");

                if (NULL == atom_u)  continue;

                if (atom_t->dist_sq(*atom_u) <= threshold_sq)
                    add_edge(i, j, atom_t->dist(*atom_u));
                break;
		    
            case Chain::ALL_ATOMS:
                min_sq = 4 * threshold_sq;

                for (unsigned k=0; k<ch->residues[i].atoms.size(); k++)  {
                    for (unsigned l=0; l<ch->residues[j].atoms.size(); l++)  {
                        float d = ch->residues[i].atoms[k].dist_sq(ch->residues[j].atoms[l]);

                        if (d < min_sq)  min_sq = d;
                    }
                }
                if (min_sq <= threshold_sq)
                    add_edge(i, j, sqrt(min_sq));
                break;

            case Chain::ALL_VDW_RADIUS:
                min_dist = 2 * threshold;

                for (unsigned k=0; k<ch->residues[i].atoms.size(); k++)  {
                    for (unsigned l=0; l<ch->residues[j].atoms.size(); l++)  {
                        float d = ch->residues[i].atoms[k].dist(ch->residues[j].atoms[l]) -
                            ch->residues[i].atoms[k].radius() - ch->residues[j].atoms[l].radius();

                        if (d < min_dist)  min_dist = d;
                    }
                }

                if (min_dist <= threshold)
                    add_edge(i, j, min_dist);
            }

        }
    }
}


ProteinContactGraph ProteinContactGraph::filter_sphere(const string &seqres, float radius)  {
    Residue *residue_t = NULL;
    residue_t = ch->get_residue(seqres);

    if (NULL == residue_t)  {
        cerr << "ERROR: Residue " << seqres << " not found in chain " << ch->chain_id << ".\n" << endl;
        exit(1);
    }

    ProteinContactGraph temp;

    // Find all residues in the sphere.
    float radius_sq = radius * radius;
    Atom *center = residue_t->get_atom("CA");

    if (NULL == center)  {
        cerr << "ERROR: Residue " + seqres + " does not have a CA atom.\n" << endl;
        exit(1);
    } 

    queue<unsigned> neighbors;
    set<unsigned> seen;
    vector<unsigned> vertex_map(ch->residues.size(), UINT_MAX);
    unsigned current(0);

    // All entries in the adjacency lists are indexes in the ch->residues
    // vector, so we need to find the index of the center residue.   
    for (unsigned i=0; i<ch->residues.size(); i++)  {
        if (seqres == ch->residues[i].residue_seq())  {
            neighbors.push(i);
            break;
        }
    }

    while (!neighbors.empty())  {
        unsigned i = neighbors.front();
            
        if (seen.find(i) == seen.end())  {
            seen.insert(i);

            if (NULL != ch->residues[i].get_atom("CA") && center->dist_sq(*ch->residues[i].get_atom("CA")) <= radius_sq)  {
                vertex_map[i] = current;

                temp.nodes.push_back(ch->residues[i].residue_id()); 
                temp.edges.resize(current + 1);
                temp.distances.resize(current + 1);

                for (unsigned j=0; j<edges[i].size(); j++)  {
                    unsigned k = edges[i][j];

                    if (vertex_map[k] < UINT_MAX)
                        temp.add_edge(current, vertex_map[k], distances[i][j]);

                    if (seen.find(k) == seen.end())
                        neighbors.push(k);
                }
                current++;
            }
        }

        neighbors.pop();
    }

    return temp;
}


vector<int> ProteinContactGraph::describe_shell(const string& seqres, float r1, float r2)  {
    vector<int> temp(32,0);

    Residue *residue_t = NULL;
    residue_t = ch->get_residue(seqres);

    if (NULL == residue_t)  {
        cerr << "ERROR: Residue " << seqres << " not found in chain " << ch->chain_id << ".\n" << endl;
        exit(1);
    }

    // Find all residues in the sphere.
    float r1_sq = r1 * r1;
    float r2_sq = r2 * r2;
    Atom *center = residue_t->get_atom("CA");

    if (NULL == center)  {
        cerr << "ERROR: Residue " + seqres + " does not have a CA atom.\n" << endl;
        exit(1);
    } 

    for (unsigned i=0; i < ch->residues.size(); i++)  {
        if (NULL != ch->residues[i].get_atom("CA"))  { 
            float dist = center->dist_sq(*ch->residues[i].get_atom("CA"));

            if (dist > r1_sq && dist <= r2_sq)  {
                // The standard 20 amino acids ordered in the descending order of their 
                // frequency in PDB Select 25 and (approximetely) SwissProt. 

                if ("LEU" == ch->residues[i].res_name)  {
                    temp[9]++;
                    temp[20]++;  // Aliphatic: A, V, L, I
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                } else if ("ALA" == ch->residues[i].res_name) {
                    temp[0]++;
                    temp[20]++;  // Aliphatic: A, V, L, I
                    temp[31]++;  // Small: A, G, C, S
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                } else if ("GLY" == ch->residues[i].res_name) {
                    temp[5]++;
                    temp[31]++;  // Small: A, G, C, S
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                } else if ("VAL" == ch->residues[i].res_name) {
                    temp[17]++;
                    temp[20]++;  // Aliphatic: A, V, L, I
                    temp[30]++;  // Hydrophilic: R, N, D, E, K, S, T, V 
                } else if ("GLU" == ch->residues[i].res_name) {
                    temp[3]++;
                    temp[24]++;  // Acidic: D, E
                    temp[26]++;  // Charged: D, E, R, K, H
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[30]++;  // Hydrophilic: R, N, D, E, K, S, T, V 
                } else if ("LYS" == ch->residues[i].res_name) {
                    temp[8]++;
                    temp[25]++;  // Basic: K, R, H
                    temp[26]++;  // Charged: D, E, R, K, H
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[30]++;  // Hydrophilic: R, N, D, E, K, S, T, V 
                } else if ("SER" == ch->residues[i].res_name) {
                    temp[15]++;
                    temp[21]++;  // Hydroxyl-containing: S, T, Y
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[31]++;  // Small: A, G, C, S
                    temp[30]++;  // Hydrophilic: R, N, D, E, K, S, T, V 
                } else if ("ASP" == ch->residues[i].res_name) {
                    temp[2]++;
                    temp[24]++;  // Acidic: D, E
                    temp[26]++;  // Charged: D, E, R, K, H
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[30]++;  // Hydrophilic: R, N, D, E, K, S, T, V 
                } else if ("THR" == ch->residues[i].res_name) {
                    temp[16]++;
                    temp[21]++;  // Hydroxyl-containing: S, T, Y
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[30]++;  // Hydrophilic: R, N, D, E, K, S, T, V 
                } else if ("ILE" == ch->residues[i].res_name) {
                    temp[7]++;
                    temp[20]++;  // Aliphatic: A, V, L, I
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                } else if ("ARG" == ch->residues[i].res_name) {
                    temp[14]++;
                    temp[25]++;  // Basic: K, R, H
                    temp[26]++;  // Charged: D, E, R, K, H
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[30]++;  // Hydrophilic: R, N, D, E, K, S, T, V 
                } else if ("ASN" == ch->residues[i].res_name) {
                    temp[11]++;
                    temp[22]++;  // Amide-containing: N, Q
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[30]++;  // Hydrophilic: R, N, D, E, K, S, T, V 
                } else if ("PRO" == ch->residues[i].res_name) {
                    temp[12]++;
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                } else if ("PHE" == ch->residues[i].res_name) {
                    temp[4]++;
                    temp[27]++;  // Aromatic: F, Y, W
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                } else if ("GLN" == ch->residues[i].res_name) {
                    temp[13]++;
                    temp[22]++;  // Amide-containing: N, Q
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                } else if ("TYR" == ch->residues[i].res_name) {
                    temp[19]++;
                    temp[21]++;  // Hydroxyl-containing: S, T, Y
                    temp[27]++;  // Aromatic: F, Y, W
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                } else if ("HIS" == ch->residues[i].res_name) {
                    temp[6]++;
                    temp[25]++;  // Basic: K, R, H
                    temp[26]++;  // Charged: D, E, R, K, H
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                } else if ("MET" == ch->residues[i].res_name) {
                    temp[10]++;
                    temp[23]++;  // Sulfur-containing: C, M
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                } else if ("CYS" == ch->residues[i].res_name)  {
                    temp[1]++;
                    temp[23]++;  // Sulfur-containing: C, M
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[31]++;  // Small: A, G, C, S
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                } else if ("TRP" == ch->residues[i].res_name)  { 
                    temp[18]++;
                    temp[27]++;  // Aromatic: F, Y, W
                    temp[28]++;  // Polar: R, N, D, C, E, Q, H, K, S, T, W, Y
                    temp[29]++;  // Hydrophobic: A, C, G, I, L, M, F, P, W, Y
                }
            }
        }
    }

    return temp;
}


void ProteinContactGraph::write_dot_file(const char *file)  {
    ofstream out(file, ios::out);

    out << "graph PCG {\n";
    out << "    node [shape=circle,style=filled,color=lightgray]; {node [label=\"" << nodes[0] << "\"] n0; }\n";
    out << "\n";

    for (unsigned i=0; i < edges.size(); i++)
        for (unsigned j=0; j < edges[i].size(); j++)
            if (i < edges[i][j])
                out << "n" << i << " -- n" << edges[i][j] << ";\n";

    out << "}";

    out.close();
}


void ProteinContactGraph::write_dot_file(const char *file, const vector<unsigned> &dist, unsigned bound)  {
    ofstream out(file, ios::out);

    out << "graph PCG {\n";
    out << "    graph [splines=true overlap=false]\n";
    out << "    node [shape=doublecircle]; {node [label=\"" << nodes[0] << "\"] n0; }\n";
    out << "    node [shape=circle];\n";
    out << endl;

    for (unsigned i=0; i<nodes.size(); i++)
        if (dist[i]<bound)
            out << "    n" << i << " [label=\"" << nodes[i] << "\"];\n";
    out << endl;

    for (unsigned i=0; i<bound; i++)  {
        out << "{rank=same;"; 
        for (unsigned j=0; j<dist.size(); j++)
            if (i == dist[j])
                out << " n" << j; 
        out << ";}\n";  
    }
    out << endl;

    for (unsigned i=0; i<edges.size(); i++)
        for (unsigned j=0; j<edges[i].size(); j++)
            if (i < edges[i][j] && dist[i] < bound && dist[edges[i][j]] < bound)  {
                out << "    n" << i << " -- n" << edges[i][j] << ";\n";
            }

    out << "}";

    out.close();
}


void ProteinContactGraph::write_edge_list(const char *file)  {
    ofstream out(file, ios::out);

    for (unsigned i=0; i < edges.size(); i++)  {
        for (unsigned j=0; j < edges[i].size(); j++)  {
            if (i < edges[i][j])  {
                out << nodes[i] << "\t" << nodes[edges[i][j]] << "\t" << distances[i][j] << "\t";

                if (ch->has_residues() && edges[i][j]==i-1)
                    out << "backbone" << endl;
                else
                    out << endl;
            }
        }
    }

    out.close();
}    


void ProteinContactGraph::write_compact(const char *file)  {
    ofstream out(file, ios::out);

    for (unsigned i=0; i<nodes.size(); i++)
        out << AminoAcid::three_to_one(nodes[i].substr(0,3)); 
    out << endl;

    for (unsigned i=0; i<edges.size(); i++)  {
        out << i;

        for (unsigned j=0; j<edges[i].size(); j++)  {
            out << "\t" << edges[i][j];
        }
        out << endl;
    }

    out.close();
}


void ProteinContactGraph::add_edge(unsigned i, unsigned j, float d)  {
    edges[i].push_back(j);
    distances[i].push_back(d);

    edges[j].push_back(i);
    distances[j].push_back(d);
}
