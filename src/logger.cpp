#include "logger.hpp"
#include <chrono>
#include <ctime>
#include <string>

Logger::Logger(Logger &&other) { log_stream = std::move(other.log_stream); }

Logger::Logger(const std::string &log_filename) {
  log_stream.open(log_filename, std::ios::app);
}

Logger::~Logger() { log_stream.close(); }

void Logger::log_with_type(const std::string &msg) const {

  std::time_t time =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  log_stream <<
}
