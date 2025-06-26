/* Auxiliary string functions
 *
 * Vladimir Vacic
 * Algorithms and Computational Biology Lab
 * Department of Computer Science and Engineering
 * University of California, Riverside
 * Riverside, CA 92521, USA
 *
 * Apr-2-2008
 */

#ifndef __STRING_H__
#define __STRING_H__

#include <cctype>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
using namespace std;


/** String tokenizer. Splits the argument string using the delimiter
 *  parameter. */
vector<string> split(const string&, const string&);

/** Removes white space characters from begining and end of the string. */
string strip(const string &s);

inline int to_i(const string &s)  {
    return atoi(s.c_str());
}

inline double to_f(const string &s)  {
    return atof(s.c_str());
}

inline string to_s(int t) {
    ostringstream s;
    s << t;
    return s.str();
}

#endif