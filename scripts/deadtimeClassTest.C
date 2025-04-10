#include<iostream>
#include<vector>
#include<algorithm>

class TimeRangeManager {

  public:
    struct TimeRange {
      int startTime;
      int endTime;

      bool contains(int time) const {
        return time >= startTime && time <= endTime;
      }

      bool operator<(const TimeRange& other) const {
        return startTime < other.startTime;
      }
    };

    // Add a new time range
    void addRange(int startTime, int endTime) {
      if (startTime > endTime) std::swap(startTime, endTime);  // ensure correct order
      ranges.push_back({startTime, endTime});
      isMerged = true;
    }

    // Check if a time is within any merged range
    bool contains(int time) {
      if (isMerged) mergeRanges();
        
      for (const auto& range : mergedRanges) {
        if (range.contains(time)) return true;
      }
      return false;
    }

    // Get merged ranges (e.g., for debugging or output)
    const std::vector<TimeRange>& getMergedRanges() {
      if (isMerged) mergeRanges();
      return mergedRanges;
    }

  private:
    std::vector<TimeRange> ranges;
    std::vector<TimeRange> mergedRanges;
    bool isMerged = false;

    // Merge overlapping and adjacent ranges
    void mergeRanges() {
      if (ranges.empty()) {
        mergedRanges.clear();
        isMerged = false;
        return;
      }

      std::sort(ranges.begin(), ranges.end());
      mergedRanges.clear();
      mergedRanges.push_back(ranges[0]);

      for (size_t i = 1; i < ranges.size(); ++i) {
        TimeRange& last = mergedRanges.back();
        const TimeRange& current = ranges[i];

        if (current.startTime <= last.endTime) {
          last.endTime = std::max(last.endTime, current.endTime);
        } else {
          mergedRanges.push_back(current);
        }
      }

      isMerged = false;
    }
};


int deadtimeClassTest() {
  TimeRangeManager manager;

  manager.addRange(100, 200);
  manager.addRange(200, 300);
  manager.addRange(290, 400);
  manager.addRange(500, 600);
  manager.addRange(601, 620);

  int timeToCheck = 50;
  if (manager.contains(timeToCheck)) {
    std::cout << "Time " << timeToCheck << " is within a range.\n";
  } else {
    std::cout << "Time " << timeToCheck << " is NOT within any range.\n";
  }

  std::cout << "Merged Ranges:\n";
  for (const auto& r : manager.getMergedRanges()) {
    std::cout << "[" << r.startTime << ", " << r.endTime << "]\n";
  }

  return 0;
}
