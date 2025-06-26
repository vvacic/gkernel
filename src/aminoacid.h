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

#ifndef __AMINO_H__
#define __AMINO_H__

#include <string>
#include <vector>
using namespace std;


class AminoAcid  {
public:
    typedef enum ReductionScheme  {
        NO_REDUCTION,
        UNLABELED,
        BLOSUM_2,
        BLOSUM_3,
        BLOSUM_4,
        BLOSUM_5,
        BLOSUM_6,
        BLOSUM_8,
        BLOSUM_10,
        BLOSUM_15,
        REDUCTION_SIZE
     } ReductionScheme;

    /** Exits program with an error message if the parameter is not a 
     *  valid amino acid character. */
    static void validate_or_fail(char);

    /** One to three letter amino acide code conversion. */
    static string one_to_three(char);
    
    /** Three to one letter amino acide code conversion. */
    static char three_to_one(const string&);

    /** Simplify amino acid alphabet based on reduction scheme that follows
     *  hierarchical clustering based on BLOSUM50 distances, as described in:
     *  
     *  Beckstette M, Homann R, Giegrich R, and Kurtz S. Fast index based 
     *  algorithms and software for matching position specific scoring matrices.
     */
    static vector<string> reduce_alphabet(const vector<string>&, ReductionScheme);

    static char BLOSUM_REDUCTION[][27];
    static char THREE_LETTER_CODE[][4];
};

#endif
