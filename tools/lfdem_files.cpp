#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <stdexcept>
#include <limits>


class lf_file {
public:
  lf_file(std::string fname);
  ~lf_file();

  std::map < std::string, std::string > get_meta_data() const;
  std::map < std::string, std::pair<int, int> > get_column_def() const;
  int col_nb();
  void rewind();

protected:
  std::ifstream f;

private:
  std::map < std::string, std::string > meta_data;
  std::map < std::string, std::pair<int, int> > columns;
  int _col_nb;
  std::streampos post_header_pos;

  std::pair<int, int> str2cols(const std::string &str);
  void parse_header();
};

class lf_data_file: public lf_file {
public:
  lf_data_file(std::string fname);
  std::vector<std::vector<double>> get_data() const;

private:
  std::vector<std::vector<double>> data;
  void read_data();
};

struct Frame {
  std::map<std::string, double>  meta_data;
  std::vector<std::vector<double>> data;
};

class lf_snapshot_file: public lf_file {
public:
  lf_snapshot_file(std::string fname);

  struct Frame get_frame(std::size_t frame_nb);
  struct Frame next_frame();

private:
  std::vector <std::streampos> frame_locations;
  bool read_frame_meta(struct Frame &frame);
  void read_frame_data(struct Frame &frame);
};


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




/**** Public lf_file methods ***************/

inline lf_file::lf_file(std::string fname):
post_header_pos(0)
{
  f.open(fname.c_str(), std::ifstream::in);
  if (!f.is_open()) {
    throw std::runtime_error("Could not open file.\n ");
  }
  parse_header();
}

inline lf_file::~lf_file(){
  if (f.is_open()) {
    f.close();
  }
}

inline std::map < std::string, std::string > lf_file::get_meta_data() const {
  return meta_data;
}

inline std::map < std::string, std::pair<int, int> > lf_file::get_column_def() const {
  return columns;
}

inline int lf_file::col_nb()
{
 return _col_nb;
}

inline void lf_file::rewind()
{
  if (post_header_pos>0) {
    f.clear();
    f.seekg(post_header_pos);
  }
}
/**** Private lf_file methods ***************/
inline std::pair<int, int> lf_file::str2cols(const std::string &str) {
  std::string sep = "-";
  auto substr = split(str, sep);
  if (substr.size() == 1) {
    return std::make_pair(std::stoi(substr[0]), std::stoi(substr[0]));
  } else {
    return std::make_pair(std::stoi(substr[0]), std::stoi(substr[1]));
  }
}

inline void lf_file::parse_header()
{
  f.seekg(f.beg);
  for (std::string line; std::getline(f, line); ) {
    if (line.empty() || f.eof()) {
      break;
    }
    std::string sep = " ";
    auto fields = split(line, sep);
    if (fields.size() > 1) {
      if (fields[0] == "#") {
        auto meta_name = fields[1];
        fields.erase(fields.begin(), fields.begin()+2);
        meta_data[meta_name] = join(fields, " ");
      } else {
        auto col_range_str = fields[0].substr(1, std::string::npos);
        if (col_range_str.back() != ':') {
          throw std::runtime_error(" Ill-formed header for column definitions.");
        }
        col_range_str = col_range_str.substr(0, col_range_str.size()-1);
        auto col_range = str2cols(col_range_str);

        fields.erase(fields.begin());
        auto col_name = join(fields, " ");
        columns[col_name] = col_range;
      }
    }
  }
  _col_nb = 0;
  for (const auto &elm: columns) {
    const auto &col_slice = elm.second;
    if (col_slice.second > _col_nb) {
      _col_nb = col_slice.second;
    }
  }
  post_header_pos = f.tellg();
}




/********** Public lf_data_file methods *************/
inline lf_data_file::lf_data_file(std::string fname): lf_file(fname) {
  read_data();
}

inline std::vector<std::vector<double>> lf_data_file::get_data() const {
  return data;
}

/********** Private lf_data_file methods *************/
inline void lf_data_file::read_data() {
  std::vector<double> record (col_nb());
  while (true) {
    for (auto &r: record) {
      f >> r;
    }
    if (!f.eof()) {
        data.push_back(record);
    } else {
        return;
    }
  };
}


/**** Public lf_snapshot_file methods ***************/
inline lf_snapshot_file::lf_snapshot_file(std::string fname): lf_file(fname)
{}

inline struct Frame lf_snapshot_file::next_frame() {
  std::streampos pos = f.tellg();
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
    f.seekg(frame_locations[frame_nb]);
    return next_frame();
  }
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

/********** Private lf_snapshot_file methods *************/
inline bool lf_snapshot_file::read_frame_meta(struct Frame &frame) {
  frame.meta_data.clear();
  std::string line;
  while(getline(f, line)) {
    if (!line.empty()) {
      break;
    }
  };
  if (line.empty()) {
    return false;
  }
  auto frame_meta = split(line, " ");
  if (frame_meta[0] != "#") {
    throw std::runtime_error(" Ill-formed frame header.");
  }
  to_double(frame_meta[1], frame.meta_data["curvilinear strain"]);
  to_double(frame_meta[2], frame.meta_data["shear disp x"]);
  to_double(frame_meta[3], frame.meta_data["rate"]);
  to_double(frame_meta[4], frame.meta_data["target stress"]);
  to_double(frame_meta[5], frame.meta_data["time"]);
  return true;
}

inline void lf_snapshot_file::read_frame_data(struct Frame &frame) {
  frame.data.clear();
  std::string line;
  std::vector<double> record (col_nb());
  std::streampos pos;

  // first line
  while(getline(f, line)) {
    if (!line.empty()) {
      break;
    }
  };
  if (line.empty()) { // eof
    return;
  }
  while (!f.eof()) {
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
    pos = f.tellg();
    getline(f, line);
  }
  f.seekg(pos);
}
