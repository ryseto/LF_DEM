#include <functional>
#include <vector>
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
	ParameterSetFactory();
	void setFromFile(const std::string& filename_parameters);
	void setParameterFromKeyValue(const std::string &keyword, 
								  const std::string &value);
	ParameterSet getParameterSet() const;
	std::vector<Dimensional::ForceScale> getForceScales() const;
	void setSystemOfUnits(const Dimensional::UnitSystem &unit_system);
private:
	std::vector<InputParameter<bool>> BoolParams;
	std::vector<InputParameter<double>> DoubleParams;
	std::vector<InputParameter<int>> IntParams;
	std::vector<InputParameter<std::string>> StrParams;
	std::vector<InputParameter<Dimensional::ForceScale>> ForceScaleParams;
	std::vector<InputParameter<Dimensional::DimensionalQty<double>>> DimValDblParams;
	std::vector<InputParameter<Dimensional::DimensionalQty<double>>> TrueDimValDblParams;
	void setDefaultValues();
};



} // namespace Parameters
#endif // #define __LF_DEM__ParameterSetFactory__