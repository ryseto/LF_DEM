#include <functional>
#include <vector>
#include <memory>
#include "DimensionalQty.h"
#include "ParameterSet.h"

#ifndef __LF_DEM__ParameterSetFactory__
#define __LF_DEM__ParameterSetFactory__

namespace Parameters {
template<typename T>
struct InputParameter
{
	std::string name_str;
	std::function<void(ParameterSet &, InputParameter<T>)> exportToParameterSet;
	T value;
};

class ParameterSetFactory {
public:
	ParameterSetFactory(Dimensional::Unit guarranted_unit);
	void setFromFile(const std::string& filename_parameters);
	void setFromStringStream(std::stringstream& ss_initial_setup);
	void setFromLine(std::string& line);
	void setParameterFromKeyValue(const std::string &keyword,
								  const std::string &value);
	ParameterSet getParameterSet() const;
	std::vector<Dimensional::ForceScale> getForceScales() const;
	void setSystemOfUnits(const Dimensional::UnitSystem &unit_system);
private:
	std::vector<InputParameter<bool>> BoolParams;
	std::vector<InputParameter<double>> DoubleParams;
	std::vector<InputParameter<int>> IntParams;
	std::vector<InputParameter<unsigned>> UIntParams;
	std::vector<InputParameter<std::string>> StrParams;
	std::vector<InputParameter<Dimensional::ForceScale>> ForceScaleParams;                 // things that define a force scale in the sense of Dimension::Unit (i.e. are associated to a suffixa and can be used as unit)
	std::vector<InputParameter<Dimensional::DimensionalQty<double>>> DimValDblParams;      // things that are defined as double in the code but input with suffix
	std::vector<InputParameter<Dimensional::DimensionalQty<double>>> TrueDimValDblParams;  // things that are defined as DimensionalValue<double> in the code (and as such have input with suffix)
	void setDefaultValues(Dimensional::Unit guarranted_unit);
	void convertParameterUnit(const Dimensional::UnitSystem &unit_system, 
							  InputParameter<Dimensional::DimensionalQty<double>> &param);
};

} // namespace Parameters
#endif // #define __LF_DEM__ParameterSetFactory__
