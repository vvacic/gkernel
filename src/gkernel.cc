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

#include "gkernel.h"
#include "string.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>


bool GraphletKernel::VERBOSE = false;
bool GraphletKernel::DEBUG = false;
bool GraphletKernel::NORMALIZE = false;
AminoAcid::ReductionScheme GraphletKernel::REDUCTION = AminoAcid::NO_REDUCTION;


void GraphletKernel::read_graphs(string path, const vector<string> &sites)  {
    if (GraphletKernel::VERBOSE)  cerr << "Reading:";

    for (unsigned i=0; i<sites.size(); i++)  {
        if (GraphletKernel::VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;

        graphs.push_back(ProteinContactGraph::read_compact( (path + "/" + sites[i] + ".graph").c_str() ));

        switch (GraphletKernel::REDUCTION)  {
            case AminoAcid::UNLABELED:
            case AminoAcid::BLOSUM_2:
            case AminoAcid::BLOSUM_3:
            case AminoAcid::BLOSUM_4:
            case AminoAcid::BLOSUM_5:
            case AminoAcid::BLOSUM_6:
            case AminoAcid::BLOSUM_8:
            case AminoAcid::BLOSUM_10:
            case AminoAcid::BLOSUM_15:
                graphs.back().nodes = AminoAcid::reduce_alphabet(graphs.back().nodes, GraphletKernel::REDUCTION);
                break;   
            case AminoAcid::NO_REDUCTION:
            default:
                ;
        }
    }

    if (GraphletKernel::VERBOSE)  cerr << endl;
}


void GraphletKernel::compute_matrix()  {
    if (GraphletKernel::VERBOSE)  cerr << "Computing kernel matrix:";

    kernel.resize(graphs.size());
    for (unsigned i=0; i < graphs.size(); i++)
        kernel[i].resize(i+1);

    for (unsigned i=0; i < graphs.size(); i++)  {
        if (GraphletKernel::VERBOSE && (i+1) % 10 == 0)  cerr << " " << i+1;

        hashes.push_back(count(graphs[i]));
        kernel[i][i] = distance_hash_join(hashes[i], hashes[i]);
        for (unsigned j=0; j<i; j++)  {
            kernel[i][j] = distance_hash_join(hashes[i], hashes[j]);
        }
    }

    if (GraphletKernel::VERBOSE)  cerr << endl;
}


void GraphletKernel::write_binary_kernel(const char *file)  {
    if (GraphletKernel::VERBOSE)  cerr << "Writing kernel matrix in binary format:";

    ofstream out(file, ios::out | ios::binary);
    
    unsigned g_size = graphs.size();

    out.write((char*) &g_size, sizeof(unsigned)); 

    for (unsigned i=0; i < g_size; i++)  {
        for (unsigned j=0; j <= i; j++)  {
            out.write((char*) &kernel[i][j], sizeof(float));
        }
    }

    out.close();

    if (GraphletKernel::VERBOSE)  cerr << endl;
}


void GraphletKernel::write_sparse_matlab(const char *file)  {
    if (GraphletKernel::VERBOSE)  cerr << "Writing kernel matrix in Matlab format:";

    ofstream out(file, ios::out);

    for (unsigned i=0; i<graphs.size(); i++)  {
        if (GraphletKernel::VERBOSE && (i+1) % 10 == 0)  cerr << " " << i+1;

        map<Key,float> g_hash = count(graphs[i]);

        for (map<Key,float>::iterator it = g_hash.begin(); it != g_hash.end(); it++)
            out << i << "\t" << it->first << "\t" << it->second << "\t" << print_key(it->first) << endl;
    }
    out.close();

    if (GraphletKernel::VERBOSE)  cerr << endl;
}


void GraphletKernel::write_sparse_svml(const char *file)  {
    if (GraphletKernel::VERBOSE)  cerr << "Writing kernel matrix in SVM^Light format:";

    ofstream out(file, ios::out);

    for (unsigned i=0; i<graphs.size(); i++)  {
        if (GraphletKernel::VERBOSE && (i+1) % 10 == 0)  cerr << " " << i+1;

        out << labels[i];

        map<Key,float> g_hash = count(graphs[i]);
        for (map<Key,float>::iterator it = g_hash.begin(); it != g_hash.end(); it++)
            out << " " << it->first << ":" << it->second;

        out << " #" << i << endl;
    }
    out.close();

    if (GraphletKernel::VERBOSE)  cerr << endl;
}


void GraphletKernel::write_text_matrix(const char *file)  {
    if (GraphletKernel::VERBOSE)  cerr << "Writing kernel matrix in plain text format:";

    ofstream out(file, ios::out);

    for (unsigned i=0; i < graphs.size(); i++)  {
        for (unsigned j=0; j <= i; j++)   { 
            out  << kernel[i][j]  <<  "\t"  ;
        }
        out  <<  endl;
    }

    out.close();

    if (GraphletKernel::VERBOSE)  cerr << endl;
}


void GraphletKernel::write_labels(const char *file)  {
    unsigned num_pos(0), num_neg(0);

    ofstream out(file, ios::out);

    out << labels[0];
    labels[0]==1 ? num_pos++ : num_neg++;

    for (unsigned i=1; i<labels.size(); i++)  {
        out << "\n" << labels[i];
        labels[i]==1 ? num_pos++ : num_neg++;
    }
    out.close();

    if (GraphletKernel::VERBOSE)  {
        cerr << "Positive examples: " << num_pos << endl;
        cerr << "Negative examples: " << num_neg << endl;
    }
}


string GraphletKernel::print_key(Key k)  {
    unsigned orbit = k & 15;  k = k >> 4;

    #if PIVOT
        string temp(4,ZERO_CHAR);
        temp[3] = (k & 31) + ZERO_CHAR;  k = k >> 5;
    #else
        string temp(3,ZERO_CHAR);
    #endif

    temp[2] = (k & 31) + ZERO_CHAR;  k = k >> 5;
    temp[1] = (k & 31) + ZERO_CHAR;  k = k >> 5;
    temp[0] = (k & 31) + ZERO_CHAR;

    ostringstream s;
    s << temp << "," << orbit;
    return s.str();
}


#if !PIVOT
    Key GraphletKernel::make_key(char a, char b, char c, unsigned orbit, char pivot)  {
        Key k(0);

        switch (orbit)  {
            case 1:
                // A
                k = ((unsigned) (a-ZERO_CHAR)) << 10;
                break;

            case 3: 
            case 4:
                // A-A
                if (a<b)  {
                    k = ((unsigned) (a-ZERO_CHAR)) << 5;
                    k = (k+b-ZERO_CHAR) << 5;
                } else  {
                    k = ((unsigned) (b-ZERO_CHAR)) << 5;
                    k = (k+a-ZERO_CHAR) << 5;
                }
                break;

            case 2:
                // A-B
                k = ((unsigned) (a-ZERO_CHAR)) << 5;
                k = (k+b-ZERO_CHAR) << 5;
                break;

            case 5:
            case 6:
            case 11:
                // A-B-C
                k = ((unsigned) (a-ZERO_CHAR)) << 5;
                k = (k+b-ZERO_CHAR) << 5;
                k = k+c-ZERO_CHAR;
                break;

            case 12:
            case 13:
            case 14:
                // A-A-B
                if (a<b)  {
                    k = ((unsigned) (a-ZERO_CHAR)) << 5;
                    k = (k+b-ZERO_CHAR) << 5;
                } else  {
                    k = ((unsigned) (b-ZERO_CHAR)) << 5;
                    k = (k+a-ZERO_CHAR) << 5;
                }
                k = k+c-ZERO_CHAR;
                break;

            case 7:
            case 9:
            case 10:
                // A-B-B
                k = ((unsigned) (a-ZERO_CHAR)) << 5;
                if (b<c)  {
                    k = (k+b-ZERO_CHAR) << 5;
                    k = k+c-ZERO_CHAR;
                } else  {
                    k = (k+c-ZERO_CHAR) << 5;
                    k = k+b-ZERO_CHAR;
                }
                break;

            case 8:
            case 15:
                // A-A-A
                if (a<b)  {
                    if (b<c)  {
                        // a,b,c
                        k = ((unsigned) (a-ZERO_CHAR)) << 5;
                        k = (k+b-ZERO_CHAR) << 5;
                        k = k+c-ZERO_CHAR;
                    } else if (a<c)  { //JLM: small bug > instead of <
                        //  a,c,b
                        k = ((unsigned) (a-ZERO_CHAR)) << 5;
                        k = (k+c-ZERO_CHAR) << 5;
                        k = k+b-ZERO_CHAR;
                    } else  {
                        //  c,a,b
                        k = ((unsigned) (c-ZERO_CHAR)) << 5;
                        k = (k+a-ZERO_CHAR) << 5;
                        k = k+b-ZERO_CHAR;
                    }
                } else  {
                    if (a<c)  {
                        // b,a,c
                        k = ((unsigned) (b-ZERO_CHAR)) << 5;
                        k = (k+a-ZERO_CHAR) << 5;
                        k = k+c-ZERO_CHAR;
                    } else if (b<c)  {
                        // b,c,a
                        k = ((unsigned) (b-ZERO_CHAR)) << 5;
                        k = (k+c-ZERO_CHAR) << 5;
                        k = k+a-ZERO_CHAR;
                    } else  {
                        // c,b,a
                        k = ((unsigned) (c-ZERO_CHAR)) << 5;
                        k = (k+b-ZERO_CHAR) << 5;
                        k = k+a-ZERO_CHAR;
                    }
                }
                break;
        }
        k = (k << 4) + orbit;

        if (GraphletKernel::DEBUG)  cout << endl << print_key(k);
    
        return k;
    }
#else
    //JLM: Modified all cases such that they include the pivot as part of their label. 
    Key GraphletKernel::make_key(char a, char b, char c, unsigned orbit, char pivot)  {
        Key k(0);

        switch (orbit)  {
            case 0:
                // Pivot
                k = ((unsigned) (pivot-ZERO_CHAR)) << 15;
                break;

            case 1:
                // A
                k = ((unsigned) (pivot-ZERO_CHAR)) << 5;
                k = (k+a-ZERO_CHAR) << 10;
                break;

            case 3:
            case 4:
                // A-A
                k = ((unsigned) (pivot-ZERO_CHAR)) << 5;
                if (a<b)  {
                    k = (k+a-ZERO_CHAR) << 5;
                    k = (k+b-ZERO_CHAR) << 5;
                } else  {
                    k = (k+b-ZERO_CHAR) << 5;
                    k = (k+a-ZERO_CHAR) << 5;   
                }
                break;

            case 2:
                // A-B
                k = ((unsigned) (pivot-ZERO_CHAR)) << 5;
                k = (k+a-ZERO_CHAR) << 5;
                k = (k+b-ZERO_CHAR) << 5;
                break;

            case 5:
            case 6:
            case 11:
                // A-B-C
                k = ((unsigned) (pivot-ZERO_CHAR)) << 5;
                k = (k+a-ZERO_CHAR) << 5;
                k = (k+b-ZERO_CHAR) << 5;
                k = k+c-ZERO_CHAR;
                break;

            case 12:
            case 13:
            case 14:
                // A-A-B
                k = ((unsigned) (pivot-ZERO_CHAR)) << 5;
                if (a<b)  {
                    k = (k+a-ZERO_CHAR) << 5;
                    k = (k+b-ZERO_CHAR) << 5;
                } else  {
                    k = (k+b-ZERO_CHAR) << 5;
                    k = (k+a-ZERO_CHAR) << 5;   
                }
                k = k+c-ZERO_CHAR;
                break;

            case 7:
            case 9:
            case 10:
                // A-B-B
                k = ((unsigned) (pivot-ZERO_CHAR)) << 5;
                k = (k+a-ZERO_CHAR) << 5;
                if (b<c)  {
                    k = (k+b-ZERO_CHAR) << 5;
                    k = k+c-ZERO_CHAR;
                } else  {
                    k = (k+c-ZERO_CHAR) << 5;
                    k = k+b-ZERO_CHAR;
                }
                break;

            case 8:
            case 15:  
                // A-A-A
                k = ((unsigned) (pivot-ZERO_CHAR)) << 5;
                if (a<b)  {
                    if (b<c)  {
                        // a,b,c
                        k = (k+a-ZERO_CHAR) << 5;
                        k = (k+b-ZERO_CHAR) << 5;
                        k = k+c-ZERO_CHAR;
                    } else if (a<c)  { //JLM: small bug > instead of <
                        // a,c,b
                        k = (k+a-ZERO_CHAR) << 5;
                        k = (k+c-ZERO_CHAR) << 5;
                        k = k+b-ZERO_CHAR;
                    } else  {
                        //  c,a,b
                        k = (k+c-ZERO_CHAR) << 5;
                        k = (k+a-ZERO_CHAR) << 5;
                        k = k+b-ZERO_CHAR;
                    }
                } else  {
                    if (a<c)  {
                        // b,a,c
                        k = (k+b-ZERO_CHAR) << 5;
                        k = (k+a-ZERO_CHAR) << 5;
                        k = k+c-ZERO_CHAR;
                    } else if (b<c)  {
                        // b,c,a
                        k = (k+b-ZERO_CHAR) << 5;
                        k = (k+c-ZERO_CHAR) << 5;
                        k = k+a-ZERO_CHAR;
                    } else  {
                        // c,b,a
                        k = (k+c-ZERO_CHAR) << 5;
                        k = (k+b-ZERO_CHAR) << 5;
                        k = k+a-ZERO_CHAR;
                    }  
                }
                break;
        }
        k = (k << 4) + orbit;

        if (GraphletKernel::DEBUG)
            cout << endl << GraphletKernel::print_key(k);

        return k;
    }
#endif


void GraphletKernel::increment_hash(map<Key,float> &hash, const Key &k)  {
    map<Key,float>::iterator it;

    if ((it = hash.find(k)) == hash.end())
        hash[k] = 1;
    else
        hash[k] = it->second + 1;
}


map<Key,float> GraphletKernel::count(ProteinContactGraph &g)  {
    vector<unsigned> dist = g.breadth_first_sort();
    map<Key,float> hash;

    unsigned i, j, k;
    string a, b, c, pivot; //JLM: Added pivot variable for 1-graphlets (i.e. case 0).

    //JLM: Added code for 1-graphlets (i.e. orbit or case 0) key generation and hash table insetion.
    pivot = g.nodes[0];

    #if PIVOT
        // 1-graphlets, case 0
        increment_hash(hash, make_key(ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, 0, pivot[0]));
    #endif

    for (unsigned i_=0; i_ < g.edges[0].size(); i_++)  {
        i = g.edges[0][i_];
        a = g.nodes[i];

        // 2-graphlets, case 01
        increment_hash(hash, make_key(a[0], ZERO_CHAR, ZERO_CHAR, 1, pivot[0]));

        for (unsigned j_=0; j_<i_; j_++)  {
            // 3-graphlets, case 011
            j = g.edges[0][j_];
            b = g.nodes[j];

            bool found(false);
            unsigned t(0);
            while (!found && t < g.edges[i].size())  {
                found = (j == g.edges[i][t++]);
            }

            if (found)  {
                increment_hash(hash, make_key(a[0], b[0], ZERO_CHAR, 4, pivot[0]));
            } else  {
                increment_hash(hash, make_key(a[0], b[0], ZERO_CHAR, 3, pivot[0]));
            }

            // 4-graphlets, case 0111
            for (unsigned k_=0; k_ < j_; k_++)  {
                k = g.edges[0][k_];
                c = g.nodes[k];

                bool found_ij(false), found_ik(false), found_jk(false);
                unsigned t(0);
                while(!found_ij && t < g.edges[i].size())
                    found_ij = (j == g.edges[i][t++]);

                t = 0;
                while (!found_ik && t < g.edges[i].size())
                    found_ik = (k == g.edges[i][t++]);

                t = 0;
                while (!found_jk && t < g.edges[j].size())
                    found_jk = (k == g.edges[j][t++]);
 
                if (found_ij)  {
                    if (found_ik)  {
                        if (found_jk)  {
                            increment_hash(hash, make_key(a[0], b[0], c[0], 15, pivot[0]));
                        } else  {
                            increment_hash(hash, make_key(b[0], c[0], a[0], 14, pivot[0]));
                        }
                    }
                    else  {
                        if (found_jk)  {
                            increment_hash(hash, make_key(a[0], c[0], b[0], 14, pivot[0]));
                        } else  {
                            increment_hash(hash, make_key(c[0], a[0], b[0], 10, pivot[0]));
                        }
                    }
                }
                else  {
                    if (found_ik)  {
                        if (found_jk)  {
                            increment_hash(hash, make_key(a[0], b[0], c[0], 14, pivot[0]));
                        } else  {
                            increment_hash(hash, make_key(b[0], a[0], c[0], 10, pivot[0]));
                        }
                    }
                    else  {
                        if (found_jk)  {
                            increment_hash(hash, make_key(a[0], b[0], c[0], 10, pivot[0]));
                        } else  {
                            increment_hash(hash, make_key(a[0], b[0], c[0], 8, pivot[0]));
                        }
                    }
                }
            }
        }

        // 4-graphlets, case 0112
        for (unsigned j_=0; j_ < g.edges[0].size(); j_++)  {
            if (i_ == j_)  continue;

            j = g.edges[0][j_];
            b = g.nodes[j];

            for (unsigned k_=0; k_ < g.edges[i].size(); k_++)  {
                k = g.edges[i][k_];
                if (2 != dist[k])  continue;
                c = g.nodes[k];

                bool found_ij(false), found_jk(false);
                unsigned t(0);
                while(!found_ij && t < g.edges[i].size())
                    found_ij = (j == g.edges[i][t++]);

                t = 0;
                while (!found_jk && t < g.edges[j].size())
                    found_jk = (k == g.edges[j][t++]);
 
                if (found_ij)  {
                    if (found_jk)  {
                        //JLM: Added this extra check to remove double counting for orbit 13.  
                        if(i_ < j_)  {
                            increment_hash(hash, make_key(a[0], b[0], c[0], 13, pivot[0]));
                        }
                        //JLM: Doubling counting quick fix. 
                        else  {
                            continue;
                        }
                    } else  {
                        increment_hash(hash, make_key(b[0], a[0], c[0], 11, pivot[0]));
                    }
                } else  { 
                    if (found_jk)  {
                        //JLM: Added this extra check to remove double counting for orbit 12.
                        if(i_ < j_)  {
                            increment_hash(hash, make_key(a[0], b[0], c[0], 12, pivot[0]));
                        }
                        //JLM: Doubling counting quick fix.
                        else  {
                            continue;
                        }
                    } else  {
                        increment_hash(hash, make_key(b[0], a[0], c[0], 6, pivot[0]));
                    }
                }
            }
        }
 
        for (unsigned j_=0; j_ < g.edges[i].size(); j_++)  {
            // 3-graphlets, case 012
            j = g.edges[i][j_];
            if (2 != dist[j])  continue;
            b = g.nodes[j];

            increment_hash(hash, make_key(a[0], b[0], ZERO_CHAR, 2, pivot[0]));

            // 4-graphlets, case 0122
            for (unsigned k_=0; k_<j_; k_++)  {
                k = g.edges[i][k_];
                if (2 != dist[k])  continue; //JLM: There was a typo here, j instead of k.
                c = g.nodes[k];

                bool found(false);
                unsigned t(0);
                while (!found && t < g.edges[j].size())
                    found = (k == g.edges[j][t++]);

                if (found)  {
                    increment_hash(hash, make_key(a[0], b[0], c[0], 9, pivot[0]));
                } else  {
                    increment_hash(hash, make_key(a[0], b[0], c[0], 7, pivot[0]));
                }
            }

            // 4-graphlets, case 0123 
            for (unsigned k_=0; k_ < g.edges[j].size(); k_++)  {
                k = g.edges[j][k_];
                if(dist[k] <= 1 || k == i || k == j)  continue; //JLM: Expanded orbit 5 case to include all "local" graphlets whose path satisfies the 0123 property. 
                //if (3 != dist[k])  continue; //JLM: Previous line.
                c = g.nodes[k];

                //JLM: Verify that local path doesn't classify under previous orbits
                bool not_found(true);
            	unsigned t(0);
            	while (not_found && t < g.edges[i].size())  {
                    if(k == g.edges[i][t++])
                        not_found = false;
                    }
            	if (not_found)  {
                    increment_hash(hash, make_key(a[0], b[0], c[0], 5, pivot[0]));
            	}
            }
       }
    }

    if (GraphletKernel::NORMALIZE)
        normalize_spectral(hash);

    return hash;
}


void GraphletKernel::normalize_spectral(map<Key,float> &hash)  {
    float norm = sqrt(distance_hash_join(hash, hash));

    for (map<Key,float>::iterator it = hash.begin(); it != hash.end(); it++)
        it->second /= norm;
}


float GraphletKernel::distance_hash_join(map<Key,float> g_hash, map<Key,float> h_hash)  {
    float sum(0);

    for (map<Key,float>::iterator git = g_hash.begin(); git != g_hash.end(); git++)  {
        map<Key,float>::iterator hit = h_hash.find(git->first);

        if (hit != h_hash.end())
            sum += git->second * hit->second; 
    }
    return sum;
}
