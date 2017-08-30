#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <stdexcept>
#include <limits>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <set>
#include <assert.h>
#include "helpers.h"
#include "lfdem_files.h"


/**** Public lf_file methods ***************/

inline lf_file::lf_file(std::string fname):
_col_nb(0),
post_header_pos(0)
{
  open_file(fname);
  parse_header();
}

inline lf_file::lf_file(std::string fname,
                        std::vector<int> columns_to_read):
cols_to_read(columns_to_read),
_col_nb(0),
post_header_pos(0)
{
  std::sort(cols_to_read.begin(), cols_to_read.end());
  open_file(fname);
  parse_header();
  filter_coldef_by_number();
}

inline lf_file::lf_file(std::string fname,
                        std::vector<std::string> fields_to_read):
_col_nb(0),
post_header_pos(0)
{
  open_file(fname);
  parse_header();
  check_fields(fields_to_read);
  auto slices = get_field_slices(fields_to_read);
  cols_to_read = lfdem_helper::slices_to_cols(slices);
  filter_coldef_by_field(fields_to_read);
}

inline lf_file::~lf_file(){
  if (file_stream.is_open()) {
    file_stream.close();
  }
}


inline std::map < std::string, std::string > lf_file::get_meta_data() const {
  return meta_data;
}

inline std::map < std::string, std::pair<int, int> > lf_file::get_column_def() const {
  return columns;
}

inline int lf_file::col_nb() const
{
  if (cols_to_read.empty())
    return _col_nb;
  else
    return cols_to_read.size();
}

inline void lf_file::rewind()
{
  if (post_header_pos>0) {
    file_stream.clear();
    file_stream.seekg(post_header_pos);
  }
}

/**** Private lf_file methods ***************/
inline std::pair<int, int> lf_file::str2cols(const std::string &str) {
  std::string sep = "-";
  auto substr = lfdem_helper::split(str, sep);
  if (substr.size() == 1) {
    return std::make_pair(std::stoi(substr[0])-1, std::stoi(substr[0]));
  } else {
    return std::make_pair(std::stoi(substr[0])-1, std::stoi(substr[1]));
  }
}

inline void lf_file::parse_metadata_field(std::vector<std::string> &split_str)
{
  auto meta_name = split_str[1];
  split_str.erase(split_str.begin(), split_str.begin()+2);
  meta_data[meta_name] = lfdem_helper::join(split_str, " ");
}

inline void lf_file::filter_coldef_by_number() {
  // std::map < std::string, std::pair<int, int>  columns
  decltype(columns) filtered_cols;

  for (auto &field: columns) {
    auto field_name = field.first;
    auto slice = field.second;

    for (unsigned i=0; i<cols_to_read.size(); i++) {
      auto cd = cols_to_read[i];
      if (slice.first <= cd && cd < slice.second) { // we keep part of this field
        if (slice.second - slice.first > 1) { // spans multiple cols
          field_name += "-"+std::to_string(cd-slice.first);
        }
        filtered_cols[field_name] = std::make_pair(i, i+1);
      }
    }
  }
  columns = filtered_cols;
}


inline void lf_file::check_fields(std::vector<std::string> &fields_to_read) {
  // check validity
  for (auto &f: fields_to_read) {
    if (!columns.count(f))
      throw std::runtime_error(" No \""+f+"\" field found.");
  }

  // remove multiple entries
  std::set<std::string> unique_fields (fields_to_read.begin(), fields_to_read.end());
  fields_to_read = std::vector<std::string>(unique_fields.begin(), unique_fields.end());

  // sort by column order
  std::sort(fields_to_read.begin(), fields_to_read.end(), [this] (std::string f1, std::string f2) -> bool {return this->columns[f1].first < this->columns[f2].first;});
}


inline std::vector<std::pair<int, int>> lf_file::get_field_slices(const std::vector<std::string> &fields_to_read) {
  std::set<std::pair<int, int>> slices;
  for (auto &field: columns) {
    auto field_name = field.first;
    if (std::find(fields_to_read.begin(), fields_to_read.end(), field_name) != fields_to_read.end()) {
      auto slice = field.second;
      slices.insert(slice);
    }
  }
  return std::vector<std::pair<int, int>>(slices.begin(), slices.end());
}


inline void lf_file::filter_coldef_by_field(const std::vector<std::string> &fields_to_read) {
  decltype(columns) filtered_cols;
  int col = 0;
  for (auto &field_name: fields_to_read) {
    auto slice_width = columns[field_name].second - columns[field_name].first;
    filtered_cols[field_name] = std::make_pair(col, col+slice_width);
    col += slice_width;
  }
  columns = filtered_cols;
}


inline void lf_file::parse_coldef_field(std::vector<std::string> &split_str)
{
  auto col_range_str = split_str[0].substr(1, std::string::npos);
  if (col_range_str.back() != ':') {
    throw std::runtime_error(" Ill-formed header for column definitions.");
  }
  col_range_str = col_range_str.substr(0, col_range_str.size()-1);
  auto col_range = str2cols(col_range_str);
  if (col_range.second > _col_nb) {
    _col_nb = col_range.second;
  }
  split_str.erase(split_str.begin());
  auto col_name = lfdem_helper::join(split_str, " ");
  columns[col_name] = col_range;
}


inline void lf_file::parse_header()
{
  file_stream.seekg(file_stream.beg);
  for (std::string line; std::getline(file_stream, line); ) {
    if (line.empty() || file_stream.eof()) {
      break;
    }
    std::string sep = " ";
    auto split_str = lfdem_helper::split(line, sep);
    if (split_str.size() > 1) {
      if (split_str[0] == "#") {
        parse_metadata_field(split_str);
      } else {
        parse_coldef_field(split_str);
      }
    }
  }

  if (!cols_to_read.empty()) {
    auto max_to_read = std::max_element(cols_to_read.begin(), cols_to_read.end());
    if (_col_nb <  *max_to_read) {
      throw std::runtime_error(" Asking to read col "+std::to_string(*max_to_read)+\
                               " but file contains only "+std::to_string(_col_nb)+" columns.");
    }
    if (*std::min_element(cols_to_read.begin(), cols_to_read.end()) < 1) {
      throw std::runtime_error(" Asking to read col < 0");
    }
  }

  post_header_pos = file_stream.tellg();
}


inline void lf_file::open_file(std::string fname) {
  file_stream.open(fname.c_str(), std::ifstream::in);
  if (!file_stream.is_open()) {
    throw std::runtime_error("Could not open file.\n ");
  }
}



/********** Public lf_data_file methods *************/
inline lf_data_file::lf_data_file(std::string fname,
                                  std::vector<int> columns_to_read):
lf_file(fname, columns_to_read),
is_read(false) {}

inline lf_data_file::lf_data_file(std::string fname,
                                  std::vector<std::string> fields_to_read):
lf_file(fname, fields_to_read),
is_read(false) {}

inline lf_data_file::lf_data_file(std::string fname):
lf_file(fname),
is_read(false) {}

inline std::vector<std::vector<double>> lf_data_file::get_data() {
  read_data();
  return data;
}

/********** Private lf_data_file methods *************/
inline void lf_data_file::read_data() {
  if (cols_to_read.empty()) {
    read_data_full();
  } else {
    read_data_cols();
  }
}


inline void lf_data_file::read_data_full() {
  assert(cols_to_read.empty());
  std::vector<double> record (col_nb());
  while (true) {
    for (auto &r: record) {
      file_stream >> r;
    }
    if (!file_stream.eof()) {
      data.push_back(record);
    } else {
      is_read = true;
      return;
    }
  };
}


inline void lf_data_file::read_data_cols() {
  assert(!cols_to_read.empty());
  std::vector<double> record (col_nb());
  std::istream_iterator<double> iit (file_stream);
  std::vector<int> col_diff (col_nb());
  std::adjacent_difference(cols_to_read.begin(), cols_to_read.end(), col_diff.begin());
  while (true) {
    for (unsigned i=0; i<col_diff.size(); i++){
      std::advance(iit, col_diff[i]);
      record[i] = *iit;
    }
    std::advance(iit, _col_nb - cols_to_read.back());
    if (!file_stream.eof()) {
        data.push_back(record);
    } else {
        is_read = true;
        return;
    }
  };
}


/**** Public lf_snapshot_file methods ***************/
inline lf_snapshot_file::lf_snapshot_file(std::string fname): lf_file(fname)
{}

inline struct Frame lf_snapshot_file::next_frame() {
  std::streampos pos = file_stream.tellg();
  struct Frame frame;
  if (read_frame_meta(frame)) {
    read_frame_data(frame);
    if (frame_locations.empty() || pos > frame_locations[frame_locations.size()-1]) {
      frame_locations.push_back(pos);
    }
  }
  return frame;
}

inline struct Frame lf_snapshot_file::get_frame(std::size_t frame_nb) {
  if (frame_locations.empty() || frame_nb > frame_locations.size()-1) {
    struct Frame frame;
    do {
      frame = next_frame();
    }  while (frame_nb > frame_locations.size()-1 && !frame.meta_data.empty());
    return frame;
  } else {
    file_stream.seekg(frame_locations[frame_nb]);
    return next_frame();
  }
}


/********** Private lf_snapshot_file methods *************/
inline bool lf_snapshot_file::read_frame_meta(struct Frame &frame) {
  frame.meta_data.clear();
  std::string line;
  while(getline(file_stream, line)) {
    if (!line.empty()) {
      break;
    }
  };
  if (line.empty()) {
    return false;
  }
  auto frame_meta = lfdem_helper::split(line, " ");
  if (frame_meta[0] != "#") {
    throw std::runtime_error(" Ill-formed frame header.");
  }
  lfdem_helper::to_double(frame_meta[1], frame.meta_data["curvilinear strain"]);
  lfdem_helper::to_double(frame_meta[2], frame.meta_data["shear disp x"]);
  lfdem_helper::to_double(frame_meta[3], frame.meta_data["rate"]);
  lfdem_helper::to_double(frame_meta[4], frame.meta_data["target stress"]);
  lfdem_helper::to_double(frame_meta[5], frame.meta_data["time"]);
  return true;
}

inline void lf_snapshot_file::read_frame_data(struct Frame &frame) {
  frame.data.clear();
  std::string line;
  std::vector<double> record (col_nb());
  std::streampos pos;

  // first line
  while(getline(file_stream, line)) {
    if (!line.empty()) {
      break;
    }
  };
  if (line.empty()) { // eof
    return;
  }
  while (!file_stream.eof()) {
    if (line[0]=='#') {
      break;
    }
    std::istringstream ss (line);
    for (auto &r: record) {
      ss >> r;
      if (ss.fail()) {
        if (r == std::numeric_limits<double>::max()) {
          throw std::runtime_error(" double overflow in line\n"+line);
        } else if (r == std::numeric_limits<double>::min()) {
          r = 0;
        } else {
          throw std::runtime_error(" unknown parsing error in line:\n"+line);
        }
      }
    }
    frame.data.push_back(record);
    pos = file_stream.tellg();
    getline(file_stream, line);
  }
  file_stream.seekg(pos);
}
