//
//  DimensionalQty.h
//  LF_DEM
//
//  Copyright (c) 2017 Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__Dimensional__
#define __LF_DEM__Dimensional__
#define _USE_MATH_DEFINES
#include <cmath> // for M_PI
#include <stdexcept>
#include <map>
#include <assert.h>
#include <iostream>

namespace Dimensional {

enum class Dimension {
  Force,
  Time,
  Viscosity,
  Stress,
  Rate,
  Velocity,
  TimeOrStrain,
  Strain,
  none
};

enum class Unit {
  hydro,
  repulsion,
  brownian,
  cohesion,
  critical_load,
  ft_max,
  kn,
  kt,
  kr,
  stress,
  sigma_zz,
  none
};

inline Unit suffix2unit(std::string s) {
  if (s=="h") {
    return Unit::hydro;
  }
  if (s=="r" || s=="repulsion") {
    return Unit::repulsion;
  }
  if (s=="b" || s=="brownian") {
    return Unit::brownian;
  }
  if (s=="c" || s=="cohesion") {
    return Unit::cohesion;
  }
  if (s=="cl" || s=="critical_load") {
    return Unit::critical_load;
  }
  if (s=="ft" || s=="ft_max") {
    return Unit::ft_max;
  }
  if (s=="kn") {
    return Unit::kn;
  }
  if (s=="kt") {
    return Unit::kt;
  }
  if (s=="kr") {
    return Unit::kr;
  }
  if (s=="s") {
    return Unit::stress;
  }
  if (s=="sz" || s=="sigma_zz") {
    return Unit::sigma_zz;
  }
  return Unit::none;
}

inline std::string unit2suffix(Unit unit) {
  if (unit==Unit::hydro) {
    return "h";
  }
  if (unit==Unit::repulsion) {
    return "r";
  }
  if (unit==Unit::brownian) {
    return "b";
  }
  if (unit==Unit::cohesion) {
    return "c";
  }
  if (unit==Unit::critical_load) {
    return "cl";
  }
  if (unit==Unit::ft_max) {
    return "ft";
  }
  if (unit==Unit::kn) {
    return "kn";
  }
  if (unit==Unit::kt) {
    return "kt";
  }
  if (unit==Unit::kr) {
    return "kr";
  }
  if (unit==Unit::stress) {
    return "s";
  }
  if (unit==Unit::sigma_zz) {
    return "sz";
  }
  if (unit==Unit::none) {
    return "";
  }
  return "";
}

template<typename T>
struct DimensionalQty {
  Dimension dimension;
  T value;
  Unit unit;
};


inline bool getSuffix(const std::string& str,
                      std::string& value,
                      std::string& suffix)
{
	size_t suffix_pos = str.find_first_of("abcdfghijklmnopqrstuvwxyz"); // omission of "e" is intended, to allow for scientific notation like "1e5h"
	value = str.substr(0, suffix_pos);
	if (suffix_pos != str.npos) {
		suffix = str.substr(suffix_pos, str.length());
		return true;
	} else {
		return false;
	}
}

inline void errorNoSuffix(std::string quantity)
{
	std::cerr << "Error : no unit scale (suffix) provided for " << quantity << std::endl; exit(1);
}


inline DimensionalQty<double> str2DimensionalQty(Dimension dimension,
                                                 std::string value_str,
                                                 std::string name)
{
	DimensionalQty<double> inv;
	inv.dimension = dimension;

	std::string numeral, suffix;
	bool caught_suffix = true;
	caught_suffix = getSuffix(value_str, numeral, suffix);
	if (!caught_suffix) {
		errorNoSuffix(name);
	}
	inv.value = stod(numeral);
	inv.unit = suffix2unit(suffix);

  if (inv.dimension == Dimension::TimeOrStrain) {
    if (inv.unit == Unit::hydro) {
      inv.dimension = Dimension::Strain;
      inv.unit = Unit::none;
    } else {
      inv.dimension = Dimension::Time;
    }
  }
	return inv;
}

class UnitSystem {
public:
  // void add(Param::Parameter param, DimensionalQty value);
  void add(Unit unit, DimensionalQty<double> quantity);
  void setInternalUnit(Unit unit);
  template<typename T> void convertToInternalUnit(DimensionalQty<T> &quantity) const;
  template<typename T> void convertFromInternalUnit(DimensionalQty<T> &quantity, Unit unit) const;
  const std::map<Unit, DimensionalQty<double>> getForceTree() const {return unit_nodes;};
  Unit getLargestUnit() const;
private:
  std::map<Unit, DimensionalQty<double>> unit_nodes;
  void convertToParentUnit(DimensionalQty<double> &node);
  void flipDependency(Unit node_name);
  void convertNodeUnit(DimensionalQty<double> &node, Unit unit);
  template<typename T> void convertUnit(DimensionalQty<T> &quantity,
                                        Unit new_unit) const; // for arbitrary Dimension
};

template<typename T>
void UnitSystem::convertUnit(DimensionalQty<T> &quantity, Unit new_unit) const
{
  /**
    \brief Convert quantity to new_unit.
    */

  // special cases
  if (quantity.dimension == Dimension::none || quantity.dimension == Dimension::Strain) {
    return;
  }
  if (quantity.dimension == Dimension::TimeOrStrain) {
    throw std::runtime_error("UnitSystem::convertUnits : cannot convert units of a TimeOrStrain quantity.");
  }

  if (!unit_nodes.count(quantity.unit)) {
    throw std::runtime_error("UnitSystem::convertUnits : cannot convert from "+unit2suffix(quantity.unit)+" unit.");
  }
  if (!unit_nodes.count(new_unit)) {
    throw std::runtime_error("UnitSystem::convertUnits : cannot convert to "+unit2suffix(new_unit)+" unit.");
  }

  // value of the old force unit in new_unit
  double force_unit_ratio = unit_nodes.at(quantity.unit).value/unit_nodes.at(new_unit).value;
  switch (quantity.dimension) {
    case Dimension::Force:
      quantity.value *= force_unit_ratio;
      break;
    case Dimension::Time:
      quantity.value /= force_unit_ratio;
      break;
    case Dimension::Rate:
      quantity.value *= force_unit_ratio;
      break;
    case Dimension::Viscosity:
      quantity.value *= 6*M_PI;  // viscosity in solvent viscosity units irrespective the unit system chosen
      break;
    case Dimension::Stress:
      quantity.value *= 6*M_PI*force_unit_ratio;
      break;
    case Dimension::Velocity:
      quantity.value *= force_unit_ratio;
      break;
    case Dimension::none: case Dimension::Strain: case Dimension::TimeOrStrain:
      // Keep these cases here even if unreachable code (see above if statements)
      // in order to avoid gcc -Wswitch warnings.
      // While default: case would achieve a similar results,
      // it would also switch off ALL -Wswitch warnings, which is dangerous.
      break;
  }
  quantity.unit = new_unit;
}

template<typename T>
void UnitSystem::convertToInternalUnit(DimensionalQty<T> &quantity) const
{
  auto internal_unit = unit_nodes.cbegin()->second.unit;
  convertUnit(quantity, internal_unit);
}

template<typename T>
void UnitSystem::convertFromInternalUnit(DimensionalQty<T> &quantity, Unit unit) const
{
  auto internal_unit = unit_nodes.cbegin()->second.unit;
  assert(quantity.unit == internal_unit);

  convertUnit(quantity, unit);
}

} // namespace Dimensional
#endif // #ifndef __LF_DEM__Dimensional__
