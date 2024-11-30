#include <iostream>

void display_progress(int progress)
{
  int progress_bar_width = 70;
  std::cout << "[";
  int cursor_position = progress_bar_width * progress;

  for (int i = 0; i < progress_bar_width; ++i) {
    if (i < cursor_position) {
      std::cout << "=";
    } else if (i == cursor_position) {
      std::cout << ">";
    } else {
      std::cout << " ";
    }
  }

  std::cout << "] " << int(progress) << " %\r";
  std::cout.flush();
}
