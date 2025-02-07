#include <fstream>
#include <string>

class Logger {
public:                            // constructors
  Logger(const Logger &) = delete; // forbid copying
  Logger(Logger &&);               // move constructor
  Logger(const std::string &);     // constructor from
                                   // file name
  ~Logger();

public: // public interface methods
  void log_warning(const std::string &) const;
  void log_message(const std::string &) const;
  void log_error(const std::string &) const;

private: // implementation details
  void log_with_type(const std::string &) const;

private: // fields
  std::ofstream log_stream;
  static bool is_initialized;
};
