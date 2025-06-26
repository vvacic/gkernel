/* Makes a graphlet kernel matrix in SVM^Light or Matlab format
 *
 * Vladimir Vacic
 * Algorithms and Computational Biology Lab
 * Department of Computer Science and Engineering
 * University of California, Riverside
 * Riverside, CA 92521, USA
 *
 * May-2-2008
 */

#include "gkernel.h"
#include "string.h"
#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;


void print_help()  {
    cout << "Usage: make_kernel -p FILE -n FILE -a PATH -[k|m|s] OUTPUT [...]\n";
    cout << "Options:\n\n";

    cout << "  -h         Displays this message.\n\n";

    cout << "  -p FILE    List of positives.\n";
    cout << "  -n FILE    List of negatives.\n";
    cout << "  -a PATH    Path to subgraph files.\n\n";

    cout << "  -r REDUCT  Alphabet reduction scheme. Can be NO_REDUCTION,\n";
    cout << "             UNLABELED, BLOSUM_2, BLOSUM_3, BLOSUM_4, BLOSUM_5,\n";
    cout << "             BLOSUM_6, BLOSUM_8, BLOSUM_10 or BLOSUM_15.\n";
    cout << "             Defaults to NO_REDUCTION.\n\n";

    cout << "  -N         Normalize the kernel matrix.\n";
    cout << "             Defaults to false.\n\n";

    cout << "  -k KERNEL  Output file for the binary kernel matrix.\n";
    cout << "   or\n";
    cout << "  -m MATLAB  Output file for the sparse attribute matrix (Matlab).\n";
    cout << "   or\n";
    cout << "  -s SVML    Output file for the sparse attribute matrix (SVM^Light).\n";
    cout << "   or\n";
    cout << "  -t TEXT    Output file for the plain text kernel matrix.\n";
    cout << "             Defaults to binary kernel matrix.\n\n";

    cout << "  -l LABELS  Output file for the example labels.\n\n";

    cout << "  -V         Verbose (prints progress messages).\n"; 
    cout << "             Defaults to false.\n\n";

    cout << "  -D         Debug (prints debug messages).\n"; 
    cout << "             Defaults to false.\n\n";
}

int main(int argc, char* argv[])  {
    typedef enum outformat  {
        KERNEL,
        MATLAB,
        SVML,
        TEXT
     } OutputFormat;

    string pos_file;
    string neg_file;
    string path;

    OutputFormat format(KERNEL);
    string output_file;
    string labels_file;

    bool VERBOSE(false);
    bool DEBUG(false);
    bool NORMALIZE(false);
    AminoAcid::ReductionScheme REDUCTION(AminoAcid::NO_REDUCTION);
    
    // Parse command line arguments.
    for (int i=1; i<argc && (argv[i])[0] == '-'; i++)  {
        switch ((argv[i])[1])  {
            case 'h': print_help(); exit(0);
            case 'p': i++; pos_file=argv[i]; break;
            case 'n': i++; neg_file=argv[i]; break;
            case 'a': i++; path=argv[i]; break;
            case 'r': 
                i++;
                if (0 == strcmp(argv[i], "NO_REDUCTION"))
                    REDUCTION = AminoAcid::NO_REDUCTION;
                else if (0 == strcmp(argv[i], "UNLABELED"))
                    REDUCTION = AminoAcid::UNLABELED;
                else if (0 == strcmp(argv[i], "BLOSUM_2"))
                    REDUCTION = AminoAcid::BLOSUM_2;
                else if (0 == strcmp(argv[i], "BLOSUM_3"))
                    REDUCTION = AminoAcid::BLOSUM_3;
                else if (0 == strcmp(argv[i], "BLOSUM_4"))
                    REDUCTION = AminoAcid::BLOSUM_4;
                else if (0 == strcmp(argv[i], "BLOSUM_5"))
                    REDUCTION = AminoAcid::BLOSUM_5;
                else if (0 == strcmp(argv[i], "BLOSUM_6"))
                    REDUCTION = AminoAcid::BLOSUM_6;
                else if (0 == strcmp(argv[i], "BLOSUM_8"))
                    REDUCTION = AminoAcid::BLOSUM_8;
                else if (0 == strcmp(argv[i], "BLOSUM_10"))
                    REDUCTION = AminoAcid::BLOSUM_10;
                else if (0 == strcmp(argv[i], "BLOSUM_15"))
                    REDUCTION = AminoAcid::BLOSUM_15;
                else  {
                    cerr << "ERROR: Unknown alphabet reduction " << argv[i] << "." << endl;
                    print_help();
                    exit(1);
                }
                break;
            case 'k': i++; format=KERNEL; output_file=argv[i]; break;
            case 'm': i++; format=MATLAB; output_file=argv[i]; break;
            case 's': i++; format=SVML; output_file=argv[i]; break;
            case 't': i++; format=TEXT; output_file=argv[i]; break;
            case 'l': i++; labels_file=argv[i]; break;
            case 'V': VERBOSE=true; break;
            case 'D': DEBUG=true; break;
            case 'N': NORMALIZE=true; break;
            default:
                cerr << "ERROR: Unknown option " << argv[i] << endl;
                print_help();
                exit(1);
        }
    }

    if (0 == output_file.size())  {
        cerr << "ERROR: Output file name not specified." << endl;
        print_help();
        exit(1);
    }

    GraphletKernel gk;

    string line;
    vector<string> sites;
    vector<int> labels;

    // Read lists of graphs.
    ifstream p(pos_file.c_str(), ios::in);
    if (p.fail()) {
        cerr << "ERROR: Positive file " << pos_file << " cannot be opened." << endl;
        exit(1);
    }
    else  {
        while(getline(p, line))  {
            vector<string> tokens = split(line, "\t");
            sites.push_back(strip(tokens[0]) + "_" + strip(tokens[1]) + "_" + strip(tokens[2]));
            labels.push_back(1);
        }
    }
    p.close();

    ifstream n(neg_file.c_str(), ios::in);
    if (n.fail())  {
        cerr << "ERROR: Negative file " << neg_file << " cannot be opened." << endl;
        exit(1);
    }
    else  {
        while(getline(n, line))  {
            vector<string> tokens = split(line, "\t");
            sites.push_back(strip(tokens[0]) + "_" + strip(tokens[1]) + "_" + strip(tokens[2]));
            labels.push_back(-1);
        }
    }
    n.close();

    if (sites.size() < 2)  {
        cerr << "ERROR: Too few graphs." << endl << endl;
        print_help();
        exit(1);
    }

    if (VERBOSE)
        GraphletKernel::set_verbose();

    if (DEBUG)
        GraphletKernel::set_debug();

    if (NORMALIZE)
        GraphletKernel::set_normalize();

    if (REDUCTION != AminoAcid::NO_REDUCTION)
        GraphletKernel::set_reduction(REDUCTION);

    gk.read_graphs(path, sites);
    gk.set_labels(labels);
    gk.compute_matrix();

    switch (format)  {
        case KERNEL:
            gk.write_binary_kernel(output_file.c_str());
            break;
        case MATLAB:
            gk.write_sparse_matlab(output_file.c_str());
            break;
        case SVML:
            gk.write_sparse_svml(output_file.c_str());
            break;
        case TEXT:
            gk.write_text_matrix(output_file.c_str());
    }

    if (labels_file.size() > 0)
        gk.write_labels(labels_file.c_str());

    exit(0);
}

