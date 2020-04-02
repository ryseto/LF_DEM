//
//  DimensionalQty.cpp
//  LF_DEM
//
//  Copyright (c) 2017 Romain Mari. All rights reserved.
//

#include "DimensionalQty.h"
#include <set>

namespace Dimensional {

	bool getSuffix(const std::string& str,
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
	
	void errorNoSuffix(std::string quantity)
	{
		std::cerr << "Error : no unit scale (suffix) provided for " << quantity << std::endl; exit(1);
	}
	
	template<>
	DimensionalQty<double> & DimensionalQty<double>::operator=(std::string value_str)
	{
		std::string numeral, suffix;
		bool caught_suffix = true;
		caught_suffix = getSuffix(value_str, numeral, suffix);
		assert(caught_suffix);
		value = std::stod(numeral);
		unit = suffix2unit(suffix);
		return *this;
	}
	
	DimensionalQty<double> str2DimensionalQty(Dimension dimension,
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

	void UnitSystem::add(Unit unit, DimensionalQty<double> quantity)
	{
		assert(quantity.dimension == Dimension::Force || quantity.dimension == Dimension::Stress);

		if (quantity.value == 0) {
			return;
		}
		if (quantity.dimension == Dimension::Stress) {
			quantity.value /= 6*M_PI; // at some point we have to get rid of this weird unit choice
		}
		std::cout << " UnitSystem : adding force scale " << Dimensional::unit2suffix(unit) << std::endl;
		unit_nodes[unit] = quantity;
		auto parent_node_name = quantity.unit;
		// no orphans!
		if (unit_nodes.count(parent_node_name) == 0) {
			unit_nodes[parent_node_name] = {Dimension::Force, 1, parent_node_name};
		}
	}

	void UnitSystem::convertToParentUnit(DimensionalQty<double> &node)
	{
		auto &parent_node = unit_nodes[node.unit];
		node.value *= parent_node.value;
		node.unit = parent_node.unit;
	}

	void UnitSystem::convertNodeUnit(DimensionalQty<double> &node, Unit unit)
	{
		if (node.unit != unit) {
			auto &parent_node = unit_nodes[node.unit];
			if (parent_node.unit != node.unit) {
				convertNodeUnit(parent_node, unit);
				convertToParentUnit(node);
			} else { // parent_node is a root, but is not unit: the unit system is not closed
				throw std::runtime_error(" UnitSystem:: cannot express "
										 +unit2suffix(node.unit)+" in "
										 +unit2suffix(unit)+" units ");
			}
		}
	}

	void UnitSystem::flipDependency(Unit node_name)
	{
		auto &node = unit_nodes[node_name];
		const auto &parent_node_name = node.unit;
		auto &parent_node = unit_nodes[parent_node_name];
		if (node_name == parent_node_name) {
			return;
		}
		flipDependency(parent_node_name);
		if (node.value == 0) {
			throw std::runtime_error(" UnitSystem:: could not flip the unit system: "
									 +unit2suffix(node_name)+" has value 0, but value of "
									 +unit2suffix(parent_node_name)+" depends on it.");
		}
		parent_node.value = 1/node.value;
		parent_node.unit = node_name;
	}

	void UnitSystem::setInternalUnit(Unit unit)
	{
		/**
		 \brief Check force units consistency, expresses all input forces in the unit "unit".
		 */
		checkCycles();
		// the unit has a value of 1*unit (says captain obvious)
		if (unit_nodes.at(unit).unit != unit) { // if the unit force is expressed in other units than itself
			flipDependency(unit);
		}
		unit_nodes[unit] = {Dimension::Force, 1, unit};
		for (auto &node: unit_nodes) {
			convertNodeUnit(node.second, unit);
		}
		has_internal_unit = true;
	}

	void UnitSystem::checkCycles()
	{
		/**
		 \brief Checks if there is no circular dependency in the input forces.
		 */
		// We are looking for loops in a unidirectional graph defined by unit_nodes
		std::set<Unit> visited;
		for (auto node: unit_nodes) {
			auto force = node.first;
			//std::cerr << "force = " << node.first << std::endl;
			auto value = node.second;
			do {
				visited.insert(force);
				force = value.unit;
				value = unit_nodes[force];
			} while (visited.count(force) == 0);
			// now we are either at a leaf or we came back at some already visited node
			if (unit_nodes[force].unit != force) {
				throw std::runtime_error("UnitSystem::checkCycles:: Circular dependency in the forces.");
			}
			visited.clear();
		}
	}
} // namespace Dimensional
