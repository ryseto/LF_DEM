#ifndef __LF_DEM__StressComponent__
#define __LF_DEM__StressComponent__
#include <vector>
#include <string>
#include "StressTensor.h"

struct StressComponent
{
  unsigned int type;
  unsigned int rate_dependence;
  std::string group;
  std::vector<StressTensor> particle_stress;

  StressComponent(){};
  StressComponent(unsigned int _type,
                  std::size_t size,
                  unsigned int _rate_dependence,
                  const std::string &_group):
                  type(_type),
                  rate_dependence(_rate_dependence),
                  group(_group) {
    particle_stress.resize(size);
    reset();
  }

  StressTensor getTotalStress() const
  {
    StressTensor total_stress;
    for(const auto &s: particle_stress) {
      total_stress += s;
    }
    return total_stress;
  }

  void reset()
  {
    for(auto &s: particle_stress) {
      s.reset();
    }
  }
};

#endif/* defined(__LF_DEM__StressComponent__) */
