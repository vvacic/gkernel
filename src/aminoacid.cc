/* AminoAcid class. Contains utility functions for validating,
 * converting and reducing amino acid symbols.
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
#include <cctype>
#include <cstdlib>
#include <iostream>


char AminoAcid::THREE_LETTER_CODE[][4] = {
    "ALA", "XXX", "CYS", "ASP", "GLU", "PHE", "GLY",
    "HIS", "ILE", "XXX", "LYS", "LEU", "MET", "ASN", 
    "XXX", "PRO", "GLN", "ARG", "SER", "THR", "XXX",
    "VAL", "TRP", "XXX", "TYR", "XXX"
};

char AminoAcid::BLOSUM_REDUCTION[][27] = {
   //ABCDEFGHIJKLMNOPQRSTUVWXYZ
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ", //0
    "XXXXXXXXXXXXXXXXXXXXXXXXXX", //1
    "LXLEELLELXELLEXLEELLXLLXLX", //2
    "LXLEEFLELXELLEXLEELLXLFXFX", //3
    "AXLEEFAELXELLEXAEEAAXLFXFX", //4
    "AXLEEFAKLXKLLEXAEKAAXLFXFX", //5
    "AXLEEFAKLXKLLEXPEKAAXLFXFX", //6
    "AXLEEFAHLXKLLEXPEKSSXLFXFX", //8
    "AXCEEFGHLXKLLEXPEKSSXLFXFX", //10
    "AXCDEFGHLXKLLNXPQKSTXLWXFX"  //15
};

void AminoAcid::validate_or_fail(char aa)  {
    if (!isalpha(aa))  {
        cerr << "ERROR: Unknown amino acid " << aa << endl;
        exit(1);
    }

    char up = toupper(aa);

    if ('B'==up || 'J'==up || 'O'==up || 'U'==up || 'X'==up || 'Z'==up)  {
        cerr << "ERROR: Unknown amino acid " << aa << endl;
        exit(1);
    }
}

char AminoAcid::three_to_one(const string &t)  {
    // The standard 20 amino acids ordered in the descending order of their 
    // frequency in PDB Select 25 and (approximetely) SwissProt. 

    if ("LEU"==t) return 'L';	
    if ("ALA"==t) return 'A';
    if ("GLY"==t) return 'G';
    if ("VAL"==t) return 'V';
    if ("GLU"==t) return 'E';
    if ("LYS"==t) return 'K';
    if ("SER"==t) return 'S';
    if ("ASP"==t) return 'D';
    if ("THR"==t) return 'T';
    if ("ILE"==t) return 'I';
    if ("ARG"==t) return 'R';
    if ("ASN"==t) return 'N';
    if ("PRO"==t) return 'P';
    if ("PHE"==t) return 'F';
    if ("GLN"==t) return 'Q';
    if ("TYR"==t) return 'Y';
    if ("HIS"==t) return 'H';
    if ("MET"==t) return 'M';
    if ("CYS"==t) return 'C';
    if ("TRP"==t) return 'W';

    // Non-standard amino acids.

    if ("CSE"==t) return 'U';	// SELENOCYSTEINE
    if ("SEC"==t) return 'U';	// SELENOCYSTEINE
    if ("PYL"==t) return 'O';   // PYRROLYSINE	

    // Ambiguous amino acids.

    if ("ASX"==t) return 'B';	// ASPARAGINE OR ASPARTATE
    if ("GLX"==t) return 'Z';	// GLUTAMINE OR GLUTAMATE
    if ("XLE"==t) return 'J';	// LEUCINE OR ISOLEUCINE
    if ("XAA"==t) return 'X';	// UNSPECIFIED OR UNKNOWN
    if ("UNK"==t) return 'X';	// UNSPECIFIED OR UNKNOWN

    // Modified amino acids.

    if ("SEP"==t) return 'S';	// PHOSPHONOSERINE
    if ("TPO"==t) return 'T';	// PHOSPHOTHREONINE
    if ("PTR"==t) return 'Y';	// O-PHOSPHOTYROSINE
    if ("NEP"==t) return 'H';	// N1-PHOSPHONOHISTIDINE	
    if ("PHD"==t) return 'D';	// ASPARTYL PHOSPHATE	
    if ("143"==t) return 'C';	// S-2,3-DIHYDRO-5-GLYCIN-2-YL-ISOXAZOL-3-YL-CYSTEINE
    if ("CAS"==t) return 'C';	// S-(DIMETHYLARSENIC)CYSTEINE	
    if ("CGA"==t) return 'E';	// CARBOXYMETHYLATED GLUTAMIC ACID	
    if ("CME"==t) return 'C';	// S,S-(2-HYDROXYETHYL)THIOCYSTEINE	
    if ("CSD"==t) return 'A';	// 3-SULFINOALANINE
    if ("CSO"==t) return 'C';	// S-HYDROXYCYSTEINE	
    if ("CSS"==t) return 'C';	// S-MERCAPTOCYSTEINE	
    if ("CYD"==t) return 'C';	// 2-AMINO-6-(CYSTEIN-S-YL)-5-OXO-HEXANOIC ACID	
    if ("CYJ"==t) return 'K';	// (Z)-N~6~-[(4R,5S)-5-(2-CARBOXYETHYL)-4-(CARBOXYMETHYL)PIPERIDIN-3-YLIDENE]-L-LYSINE
    if ("FHL"==t) return 'K';	// (E)-N~6~-[3-CARBOXY-1-(HYDROXYMETHYL)PROPYLIDENE]-L-LYSINE
    if ("HIP"==t) return 'H';	// ND1-PHOSPHONOHISTIDINE	
    if ("KCX"==t) return 'K';	// LYSINE NZ-CARBOXYLIC ACID	
    if ("LET"==t) return 'K';	// (Z)-N^6-{3-CARBOXY-1-[(4-CARBOXY-2-OXOBUTOXY)METHYL]PROPYLIDENE}-L-LYSINE	
    if ("LLP"==t) return 'K';	// 2-LYSINE(3-HYDROXY-2-METHYL-5-PHOSPHONOOXYMETHYL-PYRIDIN-4-YLMETHANE)
    if ("LSO"==t) return 'K';	// (Z)-N~6~-(3-CARBOXY-1-{[(4-CARBOXY-2-OXOBUTYL)SULFONYL]METHYL}PROPYLIDENE)-L-LYSINE
    if ("MIS"==t) return 'S';	// MONOISOPROPYLPHOSPHORYLSERINE	
    if ("NC2"==t) return 'S';	// NITROCEFIN ACYL-SERINE	
    if ("OBS"==t) return 'K';	// (Z)-N^6-[(4S,5R)-5-(2-CARBOXYETHYL)-4-(CARBOXYMETHYL)-1-HYDROXYDIHYDRO-2H-THIOPYRANIUM-3(4H)-YLIDENE]-L-LYSINE
    if ("ONL"==t) return 'L';	// 5-OXO-L-NORLEUCINE	
    if ("S1H"==t) return 'S';	// 1-HEXADECANOSULFONYL-O-L-SERINE	
    if ("SBG"==t) return 'S';	// O-[(S)-HYDROXY(METHYL)PHOSPHORYL]-L-SERINE	
    if ("SCY"==t) return 'C';	// S-ACETYL-CYSTEINE	
    if ("SEB"==t) return 'S';	// O-BENZYLSULFONYL-SERINE	
    if ("SGB"==t) return 'S';	// O-[(S)-METHYL(1-METHYLETHOXY)PHOSPHORYL]-L-SERINE	
    if ("SGX"==t) return 'S';	// O-[(S)-AMINO(METHOXY)PHOSPHORYL]-L-SERINE	
    if ("SUN"==t) return 'S';	// O-[(R)-(DIMETHYLAMINO)(ETHOXY)PHOSPHORYL]-L-SERINE
    if ("SVW"==t) return 'S';	// O-[(R)-AMINO(HYDROXY)PHOSPHORYL]-L-SERINE	
    if ("SVX"==t) return 'S';	// O-[(R)-ETHOXY(METHYL)PHOSPHORYL]-L-SERINE	
    if ("SVY"==t) return 'S';	// O-[BIS(1-METHYLETHOXY)PHOSPHORYL]-L-SERINE	
    if ("SVZ"==t) return 'S';	// O-[(S)-HYDROXY(ISOPROPYLAMINO)PHOSPHORYL]-D-SERINE	
    if ("SXE"==t) return 'S';	// O-{(S)-ETHOXY[(1-METHYLETHYL)AMINO]PHOSPHORYL}-L-SERINE	
    if ("TYX"==t) return 'C';	// S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE	

    // Everything else: metals, true heterogens, ets. 

    return 'X';
}

string AminoAcid::one_to_three(char aa)  {
    if (isalpha(aa))
        return string(THREE_LETTER_CODE[toupper(aa) - 'A']);
    
    return string("XXX");
} 

vector<string> AminoAcid::reduce_alphabet(const vector<string> &seq, ReductionScheme scheme)  {
    vector<string> temp; 

    switch (scheme)  {
        case NO_REDUCTION: break;
        case UNLABELED:
        case BLOSUM_2:
        case BLOSUM_3:
        case BLOSUM_4:
        case BLOSUM_5:
        case BLOSUM_6:
        case BLOSUM_8:
        case BLOSUM_10:
        case BLOSUM_15:
            for (unsigned i=0; i<seq.size(); i++)  {
                if (seq[i].size() > 1)  {
                    cerr << "ERROR: Reduction " << scheme << " attepted on not single letter amino acid " << seq[i] << "." << endl;
                    exit(1);
                }
                temp.push_back(string(1, BLOSUM_REDUCTION[scheme][seq[i][0] - 'A']));
            }
            break;
        default:
            cerr << "ERROR: Reduction scheme " << scheme << " not recognized." << endl;
            exit(1);
    }
    return temp;
}