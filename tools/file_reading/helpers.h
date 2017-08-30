#ifndef __lfdem_file_reading_helpers__
#define __lfdem_file_reading_helpers__

#include <string>
#include <sstream>
#include <set>
#include <vector>


namespace lfdem_helper{

inline std::pair<std::string, std::string> split_first(const std::string& str, const std::string& sep) {
  auto found = str.find(sep);
  if (found != std::string::npos) {
    return std::make_pair(str.substr(0, found), str.substr(found+sep.size(), std::string::npos));
  } else {
    return std::make_pair(str.substr(0, found), "");
  }
}

inline std::vector<std::string> split(const std::string& str, const std::string& sep) {
  std::vector<std::string> substr;
  substr.push_back(str);
  std::string last_substr(str);
  if (sep.empty()) {
    return substr;
  }
  do {
    auto split = split_first(last_substr, sep);
    substr[substr.size()-1] = split.first;
    substr.push_back(split.second);
    last_substr = substr[substr.size()-1];
  } while(!last_substr.empty());
  substr.pop_back();
  return substr;
}

inline std::string join(const std::vector<std::string>& str, const std::string& sep) {
  std::ostringstream joined;
  bool first = true;
  for(const auto &s: str) {
    if(!first) {
      joined << sep;
    } else {
      first = false;
    }
    joined << s;
  }
  return joined.str();
}

inline std::vector<int> slice_to_cols(const std::pair<int, int> &slice)
{
  std::vector<int> cols;
  for (auto i=slice.first; i<slice.second; i++) {
    cols.push_back(i);
  }
  return cols;
}

inline std::vector<int> slices_to_cols(const std::vector<std::pair<int, int>> &slices)
{
  std::set<int> cols;
  for (auto s: slices) {
    auto new_cols = slice_to_cols(s);
    cols.insert(new_cols.begin(), new_cols.end());
  }
  return std::vector<int>(cols.begin(), cols.end());
}

inline void to_double(const std::string &s, double &d) {
  try
  {
    d = std::stod(s);
  }
  catch (std::invalid_argument e)
  {
    throw std::runtime_error("Invalid double "+s);
  }
  catch (std::out_of_range e)
  {
    auto ld = std::stold(s);
    if (ld < (long double)std::numeric_limits<double>::min()) {
      d = 0;
    }
    if (ld > (long double)std::numeric_limits<double>::max()) {
      throw std::runtime_error("double overflow "+s);;
    }
  }
}

}

#endif  // __lfdem_file_reading_helpers__
