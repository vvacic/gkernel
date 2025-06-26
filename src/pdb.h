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

#ifndef __PDB_H__
#define __PDB_H__

#include "string.h"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
using namespace std;


/*
 * Point - utility class to store a point in 3D space 
 */
class Point  {
public:
    Point()  {}
    Point(double i, double j, double k) : x(i), y(j), z(k)  {}
    ~Point()  {}

    /** Euclidean distance between two points. */
    inline double dist_sq(const Point &t)  {
        return (x-t.x)*(x-t.x) + (y-t.y)*(y-t.y) + (z-t.z)*(z-t.z);
    }

    inline double dist(const Point &t)  {
        return sqrt(dist_sq(t));
    }

    inline Point operator-(const Point &t)  {
        return Point(x-t.x, y-t.y, z-t.z);
    }

    /** Dot product. */
    inline double dot(const Point &t)  {
        return x*t.x + y*t.y + z*t.z; 
    }

    /** Cross product. */
    inline Point cross(const Point &t)  {
        return Point(y*t.z - t.y*z, z*t.x-t.z*x, x*t.y - t.x*y);
    }

    /** Coordinates of the geometric center (average co-ord). */
    static Point geometric_center(const vector<Point>&);

    /** Angle between the lnes defined by two pairs of points. */
    static double dihedral_angle(const Point&, const Point&, const Point&, const Point&);

    /** Plane defined by three points. */
    static vector<double> calculate_plane(const Point&, const Point&, const Point&);

    inline string to_s()  {
        ostringstream s;
        s << x << "\t" << y << "\t" << z;
        return s.str();
    }

    double x;
    double y;
    double z;
};

/** 
 *  Atom - wrapper class for atoms and hetatoms 
 */
class Atom  {
public:
    Atom()  {}
    Atom(string t) : type(t)  {}
    ~Atom()  {}

    /** Parse a PDB ATOM line and create an Atom object. */
    void parse(const string&);

    /** Create a PDB ATOM/HETATM line. */
    string unparse();
    
    /** Creates residue id from an ATOM (or HETATM) object. */
    inline string residue_id()  {
        return res_name + to_s(res_seq) + i_code;
    }

    /** Creates residue sequence number from an ATOM (or HETATM) object. */
    inline string residue_seq()  {
        return to_s(res_seq) + i_code;
    }

    /** Van Der Waals radii of atoms from IMB Jena Image Library - 
     *  Structural Biology Glossary
     *  http://www.fli-leibniz.de/ImgLibDoc/glossary/IMAGE_VDWR.html */
    float radius() const;
    float mass() const;

    inline Point get_coords()  { return center; }

    inline float dist(const Atom &a)  {
        return center.dist(a.center);
    }

    inline float dist_sq(const Atom &a)  {
        return center.dist_sq(a.center);
    }

    /** Coordinates of the center of mass. */
    static Point center_mass(const vector<Atom>&); 



    string type; //ATOM or HETATM
    int    serial;
    string name;
    string alt_loc;
    string res_name;
    char   chain_id;
    int    res_seq;
    string i_code;
    Point  center;
    float occupancy;
    float b_factor;
    string seg_id;
    string element;
    string charge;
};



/** 
 *  Wrapper class for amino acids, nucleotides and ligands. 
 */
class Residue  {
public:
    Residue()  {}
    Residue(string n, int s, string i) : res_name(n), res_seq(s), i_code(i)  {}
    ~Residue()  {}

    /** Residue id. */
    inline string residue_id()  {
        return res_name + to_s(res_seq) + i_code;
    }

    /** Residue sequence number. */
    inline string residue_seq()  {
        return to_s(res_seq) + i_code;
    }

    /** Get atom by atom name. Returns NULL is such atom not present. */
    Atom* get_atom(const string&);


    // Data.

    string res_name;
    int    res_seq;
    string i_code;
    char   chain;

    vector<Atom> atoms;
};



/**
 *  Wrapper class for PDB chains.
 */
class Chain  {
public:
    typedef enum ConnectionMethod  {
        C_ALPHA,
        C_BETA,
        ALL_ATOMS,
        ALL_VDW_RADIUS
    } ConnectionMethod;


    Chain()  {}
    Chain(char ch) : chain_id(ch)  {}
    ~Chain()  {}


    inline bool has_residues()  {
        return residues.size() > 0;
    }
    inline bool has_heterogens()  {
        return heterogens.size() > 0;
    }

    Residue* get_residue(const string&);
    Residue* get_heterogen(const string&);

    /** Amino acid sequence from the SEQRES lines. */
    inline string seq_res()  {
        return sequence;
    }

    /** Amino acid sequence from ATOM lines. */
    string seq_atom();


    /** Contact map for the two chains; matrix entries are smallest distances 
     *  between residues under the connection method.  
     *  
     *  NOTE: Remains to be implemented. */
    vector<vector<float> > contact_map(const Chain&, ConnectionMethod);

    /** All residues from the chain which interact with the other chain,
     *  based on the closeness criterium.
     *  
     *  NOTE: Remains to be implemented. */
    vector<bool> interacting(const Chain&, ConnectionMethod, float);

    /** For each residue in this chain returns an array of chain or heterogen 
     *  identifiers with which the residue interacts.
     *  
     *  NOTE: Remains to be implemented. */
    vector<vector<string> > interaction_map(ConnectionMethod, float);


    char   chain_id;
    string sequence;  // From the SEQRES lines. 

    vector<Residue>       residues;
    map<string, Residue*> residue_map;
    vector<Residue>       heterogens;
    map<string, Residue*> heterogen_map;    
};


/*
 * PDB - wrapper class for the PDB object
 */
class PDB  {
public:
    PDB()  {}
    ~PDB()  {}

    /** Parses a PDB flat file and creates a PDB object. */ 
    void parse(const char*);

    Chain* get_chain(char);
    void print_chain(char);

    string id;
    vector<Chain> chains;

    string classification;
    string deposit_date;
};

#endif

