#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <stdexcept>


class lf_file {
public:
  lf_file(std::string fname);
  ~lf_file();

  std::map < std::string, std::string > get_meta_data() const;
  std::map < std::string, std::pair<int, int> > get_column_def() const;

  int col_nb();


protected:
  std::ifstream f;

private:
  std::map < std::string, std::string > meta_data;
  std::map < std::string, std::pair<int, int> > columns;
  int _col_nb;

  std::pair<int, int> str2cols(const std::string &str);
  void parse_header();
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

struct Frame {
  std::map < std::string, double >  meta_data;
  std::vector < std::vector < double > > data;
};




/**** Public lf_file methods ***************/

inline lf_file::lf_file(std::string fname) {
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
  f.seekg(0);
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
}



class lf_data_file: public lf_file {
public:
  lf_data_file(std::string fname): lf_file(fname) {
    read_data();
  }
  std::vector < std::vector < double > > get_data() const {
    return data;
  }

private:
  std::vector < std::vector < double > > data;
  void read_data() {
    std::string first_line;
    do {
      getline(f, first_line);
    } while (first_line.empty());
    auto fields = split(first_line, " ");
    std::vector < double > record;
    for (const auto &field: fields) {
      record.push_back(stod(field));
    }
    do {
      data.push_back(record);
      for (auto &r: record) {
        f >> r;
      }
    } while (!f.eof());
  }
};


class lf_snapshot_file: public lf_file {
public:
  lf_snapshot_file(std::string fname): lf_file(fname) {}

  struct Frame get_frame(std::size_t frame_nb) {
    if(read_frame(frame_nb)) {
      return frame;
    } else {
      struct Frame empty_frame;
      return empty_frame;
    }
  }

  struct Frame next_frame() {
    if(read_next_frame()) {
      return frame;
    } else {
      struct Frame empty_frame;
      return empty_frame;
    }
  }

private:
  std::vector < std::streampos > frame_locations;
  struct Frame frame;

  bool read_frame_meta() {
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
    frame.meta_data["curvilinear strain"] = stod(frame_meta[1]);
    frame.meta_data["shear disp x"] = stod(frame_meta[2]);
    frame.meta_data["rate"] = stod(frame_meta[3]);
    frame.meta_data["target stress"] = stod(frame_meta[4]);
    frame.meta_data["time"] = stod(frame_meta[5]);
    return true;
  };

  void read_frame_data() {
    frame.data.clear();
    std::string line;
    std::vector < double > record;
    std::streampos pos;

    // first line
    while(getline(f, line)) {
      if (!line.empty()) {
        break;
      }
    };
    if (line.empty()) {
      return;
    }
    std::istringstream ss (line);
    while(!ss.eof()) {
      double field;
      ss >> field;
      record.push_back(field);
    }
    frame.data.push_back(record);

    pos = f.tellg();
    getline(f, line);
    while (line[0]!='#' && !f.eof()) {
      std::istringstream ss2 (line);
      ss2.str(line);
      for (auto &r: record) {
        ss2 >> r;
      }
      frame.data.push_back(record);
      pos = f.tellg();
      getline(f, line);
    }
    f.seekg(pos);
  };

  bool read_next_frame() {
    std::streampos pos = f.tellg();
    if (read_frame_meta()) {
      read_frame_data();
      if (frame_locations.empty() || pos > frame_locations[frame_locations.size()-1]) {
        frame_locations.push_back(pos);
      }
      return true;
    } else {
      return false;
    }
  };

  bool read_frame(std::size_t frame_nb) {
    if (frame_locations.empty() || frame_nb > frame_locations.size()-1) {
      bool read = false;
      do {
        read = read_next_frame();
      }  while (frame_nb > frame_locations.size()-1 && read);
      return read;
    } else {
      f.seekg(frame_locations[frame_nb]);
      return read_next_frame();
    }
  };
};
