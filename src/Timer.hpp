//    This file is part of ASMC, developed by Pier Francesco Palamara.
//
//    ASMC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ASMC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ASMC.  If not, see <https://www.gnu.org/licenses/>.

#include <chrono>

#ifndef TIMER_HPP
#define TIMER_HPP

class Timer
{
private:
  using timer_t = std::chrono::system_clock;
  using sys_time = std::chrono::time_point<timer_t>;

  sys_time prevtime, curtime;

public:
  /// constructs a timer, recording the initial time
  Timer();

  /// updates the current time and returns the time since the last update in seconds
  double update_time();
};

#endif
