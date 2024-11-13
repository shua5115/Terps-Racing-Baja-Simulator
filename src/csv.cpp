#include "csv.hpp"

// trim from start (in place)
static void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

static void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}

std::vector<std::vector<std::string>> read_csv(std::istream &stream, std::map<std::string, size_t> *header) {
    std::vector<std::vector<std::string>> data;
    std::string line;
    while(std::getline(stream, line)) {
        auto &row = data.emplace_back();
        for(auto it = line.begin(); it != line.end(); it++) {
            it = std::find_if(it, line.end(), [](unsigned char c){ return !isspace(c); });
            if (it == line.end()) break;
            if (*it == '"') {
                it++;
                auto elem_end = std::find(it, line.end(), '"');
                row.emplace_back(it, elem_end);
                it = std::find(elem_end, line.end(), ',');
            } else {
                auto elem_end = std::find(it, line.end(), ',');
                row.emplace_back(it, elem_end);
                it = elem_end;
            }
            if (it == line.end()) break;
        }
    }
    if (header != nullptr && data.size() > 0) {
        header->clear();
        const auto &row = data.at(0);
        for(size_t i = 0; i < row.size(); i++) {
            header->insert({row.at(i), i});
        }
        data.erase(data.begin());
    }
    return data;
}