/* A simplified PDB parser based on the BioRuby 1.2.1 (http://www.bioruby.org)
 * PDB classes. Parses only the following types of lines:
 *
 *   SEQRES
 *   ATOM
 *   HETATM
 *
 * All other lines are ignored. For structures obtained via NMR, by default 
 * only the first MODEL is loaded and all others are discarded. Solvent 
 * molecules (HOH) are alse excluded. 
 *
 * For format description see http://www.wwpdb.org/docs.html
 *
 * Vladimir Vacic
 * Algorithms and Computational Biology Lab
 * Department of Computer Science and Engineering
 * University of California, Riverside
 * Riverside, CA 92521, USA
 *
 * Apr-2-2008
 */

#include "aminoacid.h"
#include "pdb.h"
#include <fstream>
#include <cstdio>


/*
 * Point - utility class to store a point in 3D space 
 */

Point Point::geometric_center(const vector<Point> &coords)  {
    Point t(0, 0, 0);

    for (unsigned i=0; i<coords.size(); i++)  {
        t.x += coords[i].x;
        t.y += coords[i].y;
        t.z += coords[i].z;
    }

    t.x /= coords.size();
    t.y /= coords.size();
    t.z /= coords.size();

    return t;
}

double Point::dihedral_angle(const Point &coord1, const Point &coord2, const Point &coord3, const Point &coord4)  {
/*
   (a1,b1,c1,d) = calculatePlane(coord1,coord2,coord3)
   (a2,b2,c2)   = calculatePlane(coord2,coord3,coord4)
      
    double torsion = acos((a1*a2 + b1*b2 + c1*c2)/(Math.sqrt(a1**2 + b1**2 + c1**2) * Math.sqrt(a2**2 + b2**2 + c2**2)))
      
    if ((a1*coord4.x + b1*coord4.y + c1*coord4.z + d) < 0)
        -torsion
    else
        return torsion;
*/
    return 0;
}
      
vector<double> Point::calculate_plane(const Point &coord1, const Point &coord2, const Point &coord3)  {
    double a, b, c, d;
    vector<double> temp;

    a = coord1.y * (coord2.z - coord3.z) +
        coord2.y * (coord3.z - coord1.z) + 
        coord3.y * (coord1.z - coord2.z);
    b = coord1.z * (coord2.x - coord3.x) +
        coord2.z * (coord3.x - coord1.x) + 
        coord3.z * (coord1.x - coord2.x);
    c = coord1.x * (coord2.y - coord3.y) +
        coord2.x * (coord3.y - coord1.y) + 
        coord3.x * (coord1.y - coord2.y);
    d = -1 * ((coord1.x * (coord2.y * coord3.z - coord3.y * coord2.z)) +
             (coord2.x * (coord3.y * coord1.z - coord1.y * coord3.z)) +
             (coord3.x * (coord1.y * coord2.z - coord2.y * coord1.z)));

    temp.push_back(a);
    temp.push_back(b);
    temp.push_back(c);
    temp.push_back(d);

    return temp;
}


/** 
 *  Atom - wrapper class for atoms and hetatoms 
 */

void Atom::parse(const string &line)  {
    serial   = to_i(line.substr(6,5));
    name     = strip(line.substr(12,4));
    alt_loc  = line.substr(16,1);
    res_name = strip(line.substr(17,3));
    chain_id = line.at(21);
    res_seq  = to_i(line.substr(22,4));
    i_code   = (' '==line.at(26)) ? "" : line.substr(26,1);

    center.x = to_f(line.substr(30,8));
    center.y = to_f(line.substr(38,8));
    center.z = to_f(line.substr(46,8));

    occupancy = to_f(line.substr(54,6));
    b_factor  = to_f(line.substr(60,6));

    seg_id  = strip(line.substr(72,4));  // 72 - 75 rstrip
    element = strip(line.substr(76,2));  // 76 - 77 lstrip
    charge  = strip(line.substr(78,2));
}

string Atom::unparse()  {
    string atomname;

    switch (name.size())  {
    case 0:
        atomname = "    ";
        break;
    case 1:
        atomname = " " + name + "  ";
        break;
    case 2:
        atomname = " " + name + " ";
        break;
    case 3:
        atomname = name;
        break;
    default:
        atomname = name.substr(0,4);
    }


    // There are inconsistencies with the atom names which the following code 
    // is trying to rectify.  

/*
          case atomname.length
          when 2
            if /\A[0-9]/ =~ atomname then
              return sprintf('%-4s', atomname)
            elsif /[0-9]\z/ =~ atomname then
              return sprintf(' %-3s', atomname)
            end
          when 3
            if /\A[0-9]/ =~ atomname then
              return sprintf('%-4s', atomname)
            end
          end
*/

    if (0 == atomname.size())  {
/*
          # ambiguous case for two- or three-letter name
          elem = self.element.to_s.strip
          if elem.size > 0 and i = atomname.index(elem) then
            if i == 0 and elem.size == 1 then
              return sprintf(' %-3s', atomname)
            else
              return sprintf('%-4s', atomname)
            end
          end

          if self.kind_of?(HETATM) then
            if /\A(B[^AEHIKR]|C[^ADEFLMORSU]|F[^EMR]|H[^EFGOS]|I[^NR]|K[^R]|N[^ABDEIOP]|O[^S]|P[^ABDMORTU]|S[^BCEGIMNR]|V|W|Y[^B])/ =~
                atomname then
              return sprintf(' %-3s', atomname)
            else
              return sprintf('%-4s', atomname)
            end
          else # ATOM
            if /\A[CHONSP]/ =~ atomname then
              return sprintf(' %-3s', atomname)
            else
              return sprintf('%-4s', atomname)
            end
          end
*/
    }

    char buffer[100];

    snprintf(buffer, 
        sizeof(buffer),
        "%-6s%5d %-4s%-1s%3s %-1c%4d%-1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%-2s",
        type.c_str(),
        serial, 
        atomname.c_str(),
        alt_loc.c_str(),
        res_name.c_str(),
        chain_id,
        res_seq,
        i_code.c_str(),
        center.x, center.y, center.z,
        occupancy,
        b_factor,
        seg_id.c_str(),
        element.c_str(),
        charge.c_str()
    );
    
    return string(buffer);
}

float Atom::radius() const  {
    // NOTE: Check if 1,2,3 are H, not C !

    switch (name.at(0))  {
    case '1':
    case '2':
    case '3':
    case 'H':
        return 1.20;
    case 'C':
        return 1.70;
    case 'N':
        return 1.55;
    case 'O':
        return 1.52;
    case 'S':
    case 'P':
        return 1.80;
    case 'F':
        return 1.47;
    default:
        return 1.70;
    }
}

float Atom::mass() const  {
    switch (name.at(0))  {
        case 'H':
            return 1;
        case 'C':
            return 12;
        case 'N':
            return 14;
        case 'O':
            return 16;
        case 'S':
            return 32;
        case 'P':
            return 31;
        default:
            return 0;
    }
}

Point Atom::center_mass(const vector<Atom> &atoms)  { 
    Point t(0, 0, 0);
    float total(0);

    for (unsigned i=0; i<atoms.size(); i++)  {
        float m = atoms[i].mass();
        total += m;

        t.x += atoms[i].center.x * m;
        t.y += atoms[i].center.y * m;
        t.z += atoms[i].center.z * m;
    }

    t.x /= total;
    t.y /= total;
    t.z /= total;

    return t;
}


/*
 * Residue
 */

Atom* Residue::get_atom(const string &name)  {
    for(unsigned i=0; i<atoms.size(); i++) 
        if (name == atoms[i].name)
            return &atoms[i];

    cerr << "ERROR: Atom " << name << " not found in " << residue_seq() << "." << endl;
    return NULL;
}


/*
 * Chain
 */

Residue* Chain::get_residue(const string &residue_seq)  {

// Not elegant, but it works, unlike the alternative...

//    map<string,Residue*>::iterator it = residue_map.find(residue_seq);
//
//    if (it == residue_map.end())
//        return NULL;
//    else
//        return it->second;

    for (unsigned i=0; i<residues.size(); i++)
        if (residue_seq == residues[i].residue_seq())
            return &residues[i];

    return NULL;
}

Residue* Chain::get_heterogen(const string &residue_seq)  {
    map<string,Residue*>::iterator it = heterogen_map.find(residue_seq);

    if (it == heterogen_map.end())
        return NULL;
    else
        return it->second;
}

string Chain::seq_atom()  {
    string t(residues.size(), ' ');

    for (unsigned i=0; i<residues.size(); i++)  {
        t[i] = AminoAcid::three_to_one(residues[i].res_name);
    }
    return t;
}

vector<vector<float> > Chain::contact_map(const Chain &partner, ConnectionMethod method)  {
    vector<vector<float> > temp(residues.size());

    for (unsigned i=0; i<residues.size(); i++)  {
        temp[i].resize(partner.residues.size());

        for (unsigned j=0; j<partner.residues.size(); j++)
            temp[i][j] = -1;
    }
    return temp;
}

vector<bool> Chain::interacting(const Chain &partner, ConnectionMethod method, float threshold)  {
    vector<vector<float> > map = contact_map(partner, method);
    vector<bool> temp;

    for (unsigned i=0; i < residues.size(); i++)  { 
        bool interacts(false);

        unsigned j(0);
        while (j < partner.residues.size() && !interacts)  {
            if (map[i][j] <= threshold)
                interacts = true;
 
            j++;
        }

        temp.push_back(interacts);
    }
    return temp;
}

vector<vector<string> > Chain::interaction_map(ConnectionMethod, float)  {
    vector<vector<string> > temp;

    return temp;
}


/*
 * PDB
 */

void PDB::parse(const char *pdb_file)  {
    ifstream pdb(pdb_file, ios::in);                                       

    if (!pdb)  {
        cerr << "ERROR: PDB file " << pdb_file << " could not be opened." << endl;
        exit(1);
    }

    Atom atom_t;
    string residue_seq;

    Residue *residue_t = NULL;
    Residue *ligand_t  = NULL;
    Chain *chain_t = NULL;

    string line;
    bool model_ended(false);

    while(getline(pdb, line) && !model_ended)  {
        string key = strip(line.substr(0,6)); 

        if ("ATOM" == key || "HETATM" == key)  {
            atom_t.parse(line);
            atom_t.type = key;

            if (chain_t == NULL || chain_t->chain_id != atom_t.chain_id)  {
                chain_t = get_chain(atom_t.chain_id);

                if (chain_t == NULL)  {
                    chains.push_back( Chain(atom_t.chain_id) );
                    chain_t = &chains.back();
                } 

                residue_t = NULL;
                ligand_t  = NULL;
	    }
 
            residue_seq = atom_t.residue_seq();
        }

        // Residue or modified residue
        if ("ATOM" == key || ("HETATM" == key && 'X' != AminoAcid::three_to_one(atom_t.res_name)))  {
            if (residue_t == NULL || (residue_t->residue_seq() != residue_seq))  {
                residue_t = chain_t->get_residue(residue_seq);

                if (residue_t == NULL)  {
                    chain_t->residues.push_back( Residue(atom_t.res_name, atom_t.res_seq, atom_t.i_code) );
                    chain_t->residue_map[ residue_seq ] = &chain_t->residues.back(); 

                    residue_t = chain_t->get_residue(residue_seq);
                }
            }

	    residue_t->atoms.push_back( atom_t );
        }
        // True heterogen, but not solvent 
        else if ("HETATM"== key && "HOH" != atom_t.res_name && 'X' == AminoAcid::three_to_one(atom_t.res_name))  {
            if (ligand_t == NULL || (ligand_t->residue_seq() != residue_seq))  {
                ligand_t = chain_t->get_heterogen(residue_seq);

                if (ligand_t == NULL)  {
                    chain_t->heterogens.push_back( Residue(atom_t.res_name, atom_t.res_seq, atom_t.i_code) );
                    chain_t->heterogen_map[ residue_seq ] = &chain_t->heterogens.back(); 
			
                    ligand_t = chain_t->get_heterogen(residue_seq);
                }
            }

            ligand_t->atoms.push_back( atom_t );
        }
        else if ("TER" == key)  {
            residue_t = NULL;
            ligand_t = NULL;
            chain_t = NULL;
        }
        else if ("HEADER" == key)  {
            classification = strip(line.substr(10,40));
            deposit_date = strip(line.substr(50,9));
            id = strip(line.substr(62,4));
        }
        else if ("SEQRES" == key)  {
            char chain_id = line.at(11);

            Chain *ch = get_chain(chain_id);

            if (NULL == ch)  {
                chains.push_back(Chain(chain_id));
                ch = &chains.back(); 
            } 

            for (int i=19; i<70; i+=4)  {
                // No more residues
                if ("   " == line.substr(i,3))
                    break;

                ch->sequence.push_back(AminoAcid::three_to_one(line.substr(i,3)));
            }
        }
        else if ("ENDMDL" == key)  {
            model_ended = true;
        }
    }

    pdb.close();
}

Chain* PDB::get_chain(char ch)  {
    for (unsigned i=0; i<chains.size(); i++)
        if (chains[i].chain_id==ch)
            return &chains[i];

    return NULL;
}
