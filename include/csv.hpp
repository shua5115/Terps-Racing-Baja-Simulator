#pragma once

#include <vector>
#include <array>
#include <string>
#include <istream>
#include <map>

std::vector<std::vector<std::string>> read_csv(std::istream &stream);

// Reads a racecapture log for a set of named columns.
// Only logs rows when all columns have data
template<int N>
std::vector<std::array<double, N>> read_rc_log(std::istream &stream, const std::array<const char *, N> &colnames) {
    std::vector<std::array<double, N>> data;
    std::string line;

    std::map<size_t, size_t> indices; // maps index of each column from the source file to index in the data array
    std::getline(stream, line); // get the first line of the file (headers)
    // headers are in form: "Name"|"units"|min|max|hz,
    // we only care about the "Name" part
    std::stringstream ss(line); // this copies line, so we can overwrite it after this
    for(size_t head_index = 0; std::getline(ss, line, ','); head_index++) {
        auto i = line.find('"', 1); // find second quote
        if (i == std::string::npos) continue; // if second quote doesn't exist, skip
        line.resize(i); // truncate to: "Name
        line.erase(line.begin()); // remove first quote (shifting only length(name)-1 chars rather than the whole string)
        // line should hold: Name
        for (size_t j = 0; j < N; j++) {
            const auto testname = colnames.at(j);
            if (line == testname) {
                indices.at(head_index) = j;
                break;
            }
        }
    }
    
    // Read all remaining lines of the file
    while (std::getline(stream, line)) {
        ss = std::stringstream(line); // new stringstream for this line
        std::array<double, N> row;
        row.fill(NAN);
        for(size_t head_index = 0; std::getline(ss, line, ','); head_index++) {
            // line now contains the value as a string, parseable as a float
            if (line.size() > 0 && indices.contains(head_index)) {
                double val = atof(line.c_str()); // Note: will return 0 on failure
                row[indices[head_index]] = val;
            }
        }
        // check if any values are NaN
        if (std::find(row.begin(), row.end(), NAN) == row.end()) {
            data.push_back(row);
        }
    }
    
    return data;
}