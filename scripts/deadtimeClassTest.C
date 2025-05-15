#include<iostream>
#include<vector>
#include<algorithm>
#include<cstdint>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>  // for uint64_t
#include <TMath.h>

class TimeRangeManager {
  public:
    struct TimeRange {
      uint64_t start;
      uint64_t end;

      bool contains(uint64_t time) const {
        return time >= start && time <= end;
      }

      bool operator<(const TimeRange& other) const {
        return start < other.start;
      }
    };

    // Add a new time range
    void addRange(uint64_t start, uint64_t end) {
      if (start > end) std::swap(start, end);  // ensure proper order
      ranges.push_back({start, end});
      isMerged = false;
    }

    // Check if a time is within any merged range
    bool contains(uint64_t time) {
      if (!isMerged) mergeRanges();
      for (const auto& range : mergedRanges) {
        if (range.contains(time)) return true;
      }
      return false;
    }

    // Access the merged ranges
    const std::vector<TimeRange>& getMergedRanges() {
      if (!isMerged) mergeRanges();
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
        isMerged = true;
        return;
      }

      std::sort(ranges.begin(), ranges.end());
      mergedRanges.clear();
      mergedRanges.push_back(ranges[0]);

      for (size_t i = 1; i < ranges.size(); ++i) {
        TimeRange& last = mergedRanges.back();
        const TimeRange& current = ranges[i];

        if (current.start <= last.end) {
          last.end = std::max(last.end, current.end);
        } else {
          mergedRanges.push_back(current);
        }
      }

      isMerged = true;
    } 
};



class TimeRangeManagerLocal {
  public:

    struct TimeRange {
      uint64_t  start;
      uint64_t  end;
      double    posX;
      double    posY;

      bool contains(uint64_t time, double x, double y, double posRange) const {
        return time >= start && time <= end && TMath::Abs(x-posX) <= posRange && TMath::Abs(y-posY) <= posRange;
      }

      bool operator<(const TimeRange& other) const {
        return start < other.start;
      }
    };

    // Constructor
    TimeRangeManagerLocal(double positionThreshold=5){
      posRange = positionThreshold;
    }

    // Add a new time range
    void addRange(uint64_t start, uint64_t end, double posX, double posY) {
      if (start > end) std::swap(start, end);  // ensure proper order
      ranges.push_back({start, end, posX, posY});
    }

    // Check if a time is within any merged range
    bool contains(uint64_t time, double x, double y) {
      for (const auto& range : ranges) {
        if (range.contains(time, x, y, posRange)) return true;
      }
      return false;
    }

    // Access the ranges
    const std::vector<TimeRange>& getRanges() {
      return ranges;
    }

  private:
    double                  posRange;
    std::vector<TimeRange>  ranges;

};


void deadtimeClassTest(const int64_t timeToCheck, const double posXToCheck, const double posYToCheck, const double posRangeToCheck) {
  TimeRangeManagerLocal manager(posRangeToCheck);

  manager.addRange(1740183203410465830, 1740183203410465830+50e3, 245, 46);
  manager.addRange(1740183203516077830, 1740183203516077830+50e3, 123, 23);
  manager.addRange(1740183203610865830, 1740183203610865830+50e3, 165, 11);

  if (manager.contains(timeToCheck, posXToCheck, posYToCheck)) {
    std::cout << "Time " << timeToCheck << " is inside manager.\n";
  } else {
    std::cout << "Time " << timeToCheck << " is NOT inside manager.\n";
  }

  std::cout << "Merged Ranges:\n";
  for (const auto& r : manager.getRanges()) {
    std::cout << "[" << r.start << ", " << r.end  << ", " << r.posX << ", " << r.posY << "]\n";
  }

  std::exit(0);
}
