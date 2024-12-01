#ifndef UTILS_H
#define UTILS_H
#include <iostream>

void display_progress(int progress, const std::string& leading_str)
{
  int progress_bar_width = 30;
  int cursor_position = progress_bar_width * progress / 100;
  std::cout << leading_str;
  for (int i = 0; i < progress_bar_width; ++i) {
    if (i <= cursor_position) {
      std::cout << "█";
    } else {
      std::cout << "▒";
    }
  }

  std::cout << " " << int(progress) << " %\r";
  std::cout.flush();
}

#endif // UTILS_H
