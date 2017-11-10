#ifndef __lfdem_file_reading_files__
#define __lfdem_file_reading_files__
#include <map>
#include <vector>
#include <string>
#include <fstream>

class lf_file {
public:
  lf_file(std::string fname);
  lf_file(std::string fname,
          std::vector<int> columns_to_read);
  lf_file(std::string fname,
          std::vector<std::string> fields_to_read);
  ~lf_file();

  std::map < std::string, std::string > get_meta_data() const;
  std::map < std::string, std::pair<int, int> > get_column_def() const;
  int col_nb() const;
  void rewind();

protected:
  std::ifstream file_stream;
  std::vector<int> cols_to_read;
  int _col_nb;

private:
  std::map < std::string, std::string > meta_data;
  std::map < std::string, std::pair<int, int> > columns;
  std::streampos post_header_pos;

  std::pair<int, int> str2cols(const std::string &str);
  void parse_file_header();
  void open_file(std::string fname);
  void parse_coldef_field(std::vector<std::string> &split_str);
  void parse_metadata_field(std::vector<std::string> &split_str);
  void filter_coldef_by_field(const std::vector<std::string> &fields_to_read);
  void filter_coldef_by_number();
  void check_fields(std::vector<std::string> &fields_to_read);
  std::vector<std::pair<int, int>> get_field_slices(const std::vector<std::string> &fields_to_read);
};

class lf_data_file: public lf_file {
public:
  lf_data_file(std::string fname);
  lf_data_file(std::string fname,
               std::vector<int> columns_to_read);
  lf_data_file(std::string fname,
               std::vector<std::string> fields_to_read);
  std::vector<std::vector<double>> get_data();

private:
  std::vector<std::vector<double>> data;
  void read_data();
  void read_data_full();
  void read_data_cols();
  bool is_read;
};

struct Frame {
    std::map<std::string, double>  meta_data;
    std::vector<std::vector<double>> data;
};

enum class SnapshotFileVersion {
    oneline_hdr,
    keyword_hdr
};

class lf_snapshot_file: public lf_file {
public:
    lf_snapshot_file(std::string fname);
    lf_snapshot_file(std::string fname,
                    std::vector<int> columns_to_read);
    lf_snapshot_file(std::string fname,
                     std::vector<std::string> fields_to_read);
    struct Frame get_frame(std::size_t frame_nb);
    struct Frame next_frame();

private:
    SnapshotFileVersion version;
    std::vector <std::streampos> frame_locations;
    bool read_frame_meta(struct Frame &frame);
    bool parse_frame_header(struct Frame &frame);
    void read_frame_data_full(struct Frame &frame);
    void getFileVersion();
};

#endif //__lfdem_file_reading_files__
