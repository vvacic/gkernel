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

#include "string.h"


vector<string> split(const string &input, const string &delimiter)  {
    vector<string> tokens;

    if (!input.empty())  {
        if (delimiter.empty()) {
            for (char c : input)  {
                tokens.emplace_back(1, c);  // single character string
            }
        }
        else {
            string::size_type start(0), end(0);

            while ((end = input.find(delimiter, start)) != std::string::npos) {
                tokens.emplace_back(input.substr(start, end - start));
                start = end + delimiter.length();
            }
            tokens.emplace_back(input.substr(start));
        }
    }

    return tokens;
}

string strip(const string &input)  {
    bool found_space(false);
    unsigned start(0);

    while (start < input.size() && !found_space)  {
        if (isspace(input.at(start))) 
            start++;
        else
            found_space = true;
    } 
    
    found_space = false;
    unsigned end(input.size()-1);

    while (end > start && !found_space)  {
        if (isspace(input.at(end)))
            end--;
        else
            found_space = true;
    }

    return input.substr(start, end - start + 1);
}
