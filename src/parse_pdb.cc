/* Parses a PDB file and outputs a ProteinContactGraph  
 *
 * Vladimir Vacic
 * Computer Science and Engineering
 * University of California, Riverside
 *
 * Apr-2-2008
 */

#include "protgraph.h"
#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;


void print_help()  {
    cout << "Usage: parse_pdb [options] PDB_FILE\n";
    cout << "Options:\n\n";
    cout << "  -h           Displays this message.\n\n";

    cout << "  -c CHAIN     PDB chain.\n";
    cout << "               Defaults to the first chain found in the PDB file.\n";
    cout << "  -m METHOD    Method used for determining interacting residues.\n";
    cout << "               Can be C_ALPHA, C_BETA, ALL_ATOMS, ALL_VDW_RADIUS.\n";
    cout << "               Defaults to C_ALPHA.\n";
    cout << "  -d DISTANCE  Threshold distance in Angstroms.\n";
    cout << "               Defaults to 6A.\n\n";

    cout << "  -r RESIDUE   Only the subgraph centered at this residue.\n";
    cout << "               Residues should be specified by their residue sequence\n";
    cout << "               number and an insertion code, if one exists (resSeq and\n"; 
    cout << "               iCode, columns 23-26 and 27 in the ATOM/HETATM lines).\n";
    cout << "  -s RADIUS    Includes only the residues with C alphas within the \n";
    cout << "               sphere of the given RADIUS, centered at the C alpha\n";
    cout << "               of the -r residue.\n\n";

    cout << "  -E or -C     Output file format, edge list or compact.\n"; 
    cout << "               Defaults to edge list.\n";
    cout << "  -o OUTFILE   Output file for the results.\n";
}

int main(int argc, char* argv[])  {
    // Defaults.
    char chain(0);
    Chain::ConnectionMethod method(Chain::C_ALPHA);
    double distance(6.0);

    string residue;
    double radius(0);

    bool compact_output(false);
    string out_file;

    // Parse command line arguments.
    int i(1);
    for (; i < argc && (argv[i])[0] == '-'; i++)  {
        switch ((argv[i])[1])  {
        case 'h':   print_help();  exit(0);
        case 'c':   i++;  chain = argv[i][0];  break;
        case 'm':
            i++;
            if (0 == strcmp(argv[i], "C_ALPHA"))
                method = Chain::C_ALPHA;
            else if (0 == strcmp(argv[i], "C_BETA"))
                method = Chain::C_BETA;
            else if (0 == strcmp(argv[i], "ALL_ATOMS"))
                method = Chain::ALL_ATOMS;
            else if (0 == strcmp(argv[i], "ALL_VDW_RADIUS"))
                method = Chain::ALL_VDW_RADIUS;
            else  {
                cerr << "ERROR: Unknown method " << argv[i] << endl;
                print_help();
                exit(1);
            }
            break;
        case 'd':
            i++;  distance = to_f(argv[i]);
            if (distance <= 0)  {
                cerr << "ERROR: Distance must be a positive real, but is " << argv[i] << endl;
                print_help();
                exit(1);
            }
            break;
        case 'r':  i++;  residue = argv[i];  break;
        case 's':
            i++; radius = to_f(argv[i]);
            if (radius <= 0)  {
                cerr << "ERROR: Sphere radius must be a positive real, but is " << argv[i] << endl;
                print_help();
                exit(1);
            }
            break;
        case 'E':  compact_output=false;  break;
        case 'C':  compact_output=true;  break;
        case 'o':  i++;  out_file = argv[i];  break;
        default: 
            cerr << "ERROR: Unknown option " << argv[i] << endl;
            print_help();
            exit(1);
        }
    }

    if (out_file.size() == 0)  {
        cerr << "ERROR: Output file was not specified." << endl;
        print_help();
        exit(1);
    }

    if ("" != residue && 0 == radius)  {
        cerr << "ERROR: Switch -r requires switch -s." << endl;
        print_help();
        exit(1);
    }

    if (i >= argc)  {
        cerr << "ERROR: PDB file was not specified." << endl;
        print_help();
        exit(1);
    }

    PDB pdb;
    pdb.parse(argv[i]);

    Chain *ch;
    
    // Build protein structure graph and pass it to on
    if (0 == chain)  {
        chain = pdb.chains[0].chain_id;
        ch    = &pdb.chains[0]; 
    }
    else  {
        ch = pdb.get_chain(chain);
        if (NULL == ch)  {
            cerr << "ERROR: Chain " << chain << " could not be found in " << argv[i] << "." << endl;
            exit(1);
        }
    }

    ProteinContactGraph g(ch, method, distance); 

    // Filter subgraph if necessary 
    if ("" != residue)  {
        if (NULL == ch->get_residue(residue))  { 
            cerr << "ERROR: Residue " + residue + " not found in chain " << ch->chain_id << "." << endl;
            exit(1); 
        }

        if (radius > 0)  {
            g = g.filter_sphere(residue, radius);
        }
    }

    if (compact_output)
        g.write_compact(out_file.c_str());
    else
        g.write_edge_list(out_file.c_str());	

    exit(0);
}
