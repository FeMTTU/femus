// Timer.h - Utility class for timing different parts of a program

#ifndef Timer_H
#define Timer_H 1

#include <sys/time.h>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>

// The Timer class: all functions are inline
class Timer {
public:
  // Constructor can specify a number of unnamed activities
  Timer(unsigned int n_activities = 0) 
    : current_activity_(-1), timer_on_(false) {
    last_time_.tv_sec = 0;
    last_time_.tv_usec = 0;
    timings_.reserve(100);
    names_.reserve(100);
    for (unsigned int i = 0; i < n_activities; i++) {
      std::stringstream s;
      s << "Activity " << i;
      timings_.push_back(0.0);
      names_.push_back(s.str());
    }
  }

  // When the timer is destructed (typically at program exit), print
  // out the times spent in each activity
  ~Timer() {
    print();
  }

  // Print out the times spent in each activity
  void print() {
    double sum = 0.0;
    std::cerr << timings_.size() << " activities:\n";
    for (unsigned int i = 0; i < timings_.size(); i++) {
      std::cerr.width(10);
      std::cerr << std::right << timings_[i] << " s: " << names_[i] << "\n";
      sum += timings_[i];
    }
    std::cerr.width(10);
    std::cerr << std::right << sum << " s: Total\n";
  }

  // Register a new activity with the specified name, returning the
  // tag to be used to specify it in future, as an unsigned int
  unsigned int new_activity(const std::string& name) {
    unsigned int tag = timings_.size();
    names_.push_back(name);
    timings_.push_back(0.0);
    return tag;
  }

  // Stop timing current activity
  void stop() {
    if (timer_on_) {
      timings_[current_activity_] += split_();
    }
    timer_on_ = false;
  };

  // Start timing specified activity
  void start(unsigned int activity) {
    if (timer_on_) {
      timings_[current_activity_] += split_();
    }
    else {
      split_();
    }

    if (activity < timings_.size()) {
      current_activity_ = activity;
      timer_on_ = true;
    }
    else {
      // Activity out of range - to keep this inline function fast we
      // don't throw an exception but just don't record the time for
      // this event
      timer_on_ = false;
    }
  };

  // Return the list of timings in seconds as a constant reference to
  // a vector of doubles
  const std::vector<double>& timings() { return timings_; }

private:
  // Use Unix system call to get the time accurately
  double split_() {
    struct timeval time;
    gettimeofday(&time, NULL);
    double dsec = time.tv_sec - last_time_.tv_sec
      + 0.000001 * (time.tv_usec - last_time_.tv_usec);
    last_time_ = time;
    return dsec;
  }
  // Data
  std::vector<double> timings_;
  std::vector<std::string> names_;
  unsigned int current_activity_;
  timeval last_time_;
  bool timer_on_;
};

#endif
