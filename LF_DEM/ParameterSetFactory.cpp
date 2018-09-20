#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "global.h"
#include "ParameterSetFactory.h"
#define PARAM_INIT(name, default_value) {#name,  [](ParameterSet &p, InputParameter<decltype(ParameterSet::name)> in) {p.name = in.value;}, default_value}
#define PARAM_INIT_DIMQTY(name, default_value) {#name,  [](ParameterSet &p, InputParameter<Dimensional::DimensionalQty<decltype(ParameterSet::name)>> in) {p.name = in.value.value;}, default_value}
#define PARAM_INIT_FORCESCALE(name, default_value) {#name,  [](ParameterSet &p, InputParameter<Dimensional::ForceScale> in) {p.name = in.value.dim_qty.value;}, default_value}

namespace Parameters {

void Str2KeyValue(const std::string& str_parameter,
				  std::string& keyword,
				  std::string& value)
{
	std::string::size_type pos_equal = str_parameter.find("=");
	keyword = str_parameter.substr(0, pos_equal);
	value = str_parameter.substr(pos_equal+1);
	return;
}

ParameterSetFactory::ParameterSetFactory() 
{
	setDefaultValues();
}

void ParameterSetFactory::setDefaultValues() { 

/*================================================
=            DEFAULT PARAMETER VALUES            =
=================================================*/

	/*================================
	=            BOOLEANS            =
	================================*/
	BoolParams = \
	{
		PARAM_INIT(fixed_dt, false),
		PARAM_INIT(keep_input_strain, false),
		PARAM_INIT(monolayer, false),
		PARAM_INIT(auto_determine_knkt, false),
		PARAM_INIT(output.out_bond_order_parameter6, false),
        PARAM_INIT(output.new_material_functions, false),
		PARAM_INIT(output.out_data_vel_components, false),
		PARAM_INIT(output.out_binary_conf, false),
		PARAM_INIT(output.out_data_interaction, true),
		PARAM_INIT(output.out_data_particle, true),
		PARAM_INIT(output.origin_zero_flow, true),
		PARAM_INIT(output.log_time_interval, false),
		PARAM_INIT(output.out_na_vel, false),
		PARAM_INIT(output.out_na_disp, false),
		PARAM_INIT(output.recording_interaction_history, false)
	};

	/*===========================================
	=            DIMENSIONLESS REALS            =
	===========================================*/
	DoubleParams = \
	{
		PARAM_INIT(brownian_relaxation_time, 1),
		PARAM_INIT(strain_reversal, 1),
		PARAM_INIT(theta_shear, 0),
		PARAM_INIT(rest_threshold, 1e-4),
		PARAM_INIT(Pe_switch, 5),
		PARAM_INIT(start_adjust, 0.2),
		PARAM_INIT(max_dt_auto_det, 1e-3),
		PARAM_INIT(min_dt_auto_det, 1e-7),
		PARAM_INIT(memory_strain_avg, 0.01),
		PARAM_INIT(memory_strain_k, 0.02),
		PARAM_INIT(disp_tan_target, 0.05),
		PARAM_INIT(overlap_target, 0.05),
		PARAM_INIT(disp_max, 2e-3),
		PARAM_INIT(lub_max_gap, 0.5),
		PARAM_INIT(lub_reduce_parameter, 1e-3),
		PARAM_INIT(sd_coeff, 1),
		PARAM_INIT(interaction_range, -1),
		PARAM_INIT(repulsive_length, 0.05),
		PARAM_INIT(repulsive_max_length, -1),
		PARAM_INIT(vdW_coeffient, -1),
		PARAM_INIT(vdW_singularity_cutoff, 0.1),
		PARAM_INIT(magic_angle, 0),
		PARAM_INIT(mu_static, 1),
		PARAM_INIT(mu_dynamic, -1),
		PARAM_INIT(mu_rolling, 0),
		PARAM_INIT(TA_adhesion.adhesion_range, 1e-2),
		PARAM_INIT(output.recording_start, 1),
		PARAM_INIT(shear_jamming_rate, 1e-6)
	};

	/*================================
	=            INTEGERS            =
	================================*/
	IntParams = \
	{
		PARAM_INIT(output.nb_output_data_log_time, 100),
		PARAM_INIT(output.nb_output_config_log_time, 100),
		PARAM_INIT(integration_method, 1),
		PARAM_INIT(friction_model, 1),
		PARAM_INIT(np_fixed, 0),
		PARAM_INIT(simulation_mode, 0),
		PARAM_INIT(shear_jamming_max_count, 30)
	};

	/*===============================
	=            STRINGS            =
	===============================*/
	StrParams = \
	{
		PARAM_INIT(flow_type, ""),
		PARAM_INIT(event_handler, ""),
		PARAM_INIT(output.out_particle_stress, ""),
		PARAM_INIT(lubrication_model, "tangential")
	};

	/*====================================
	=            FORCE SCALES            =
	====================================*/
	Dimensional::ForceScale default_val;

	default_val = {Dimensional::Unit::kn, {Dimensional::Dimension::Force, 0, Dimensional::Unit::hydro}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(kn, default_val));

	default_val = {Dimensional::Unit::kt, {Dimensional::Dimension::Force, 0, Dimensional::Unit::kn}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(kt, default_val));

	default_val = {Dimensional::Unit::kr, {Dimensional::Dimension::Force, 0, Dimensional::Unit::kn}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(kr, default_val));

	default_val = {Dimensional::Unit::delayed_adhesion, {Dimensional::Dimension::Force, 0, Dimensional::Unit::hydro}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(TA_adhesion.adhesion_max_force, default_val));

	default_val = {Dimensional::Unit::repulsion, {Dimensional::Dimension::Force, 0, Dimensional::Unit::hydro}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(repulsion, default_val));

	default_val = {Dimensional::Unit::critical_load, {Dimensional::Dimension::Force, 0, Dimensional::Unit::hydro}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(critical_load, default_val));

	default_val = {Dimensional::Unit::cohesion, {Dimensional::Dimension::Force, 0, Dimensional::Unit::hydro}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(cohesion, default_val));

	default_val = {Dimensional::Unit::brownian, {Dimensional::Dimension::Force, 0, Dimensional::Unit::hydro}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(brownian, default_val));

	/*==============================================
	=            Dimensional Quantities            =
	==============================================*/

	/*----------  Double dim vals  ----------*/

	Dimensional::DimensionalQty<double> default_qty;

	default_qty = {Dimensional::Dimension::Time, 1e-4, Dimensional::Unit::hydro};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(dt, default_qty));
	
	default_qty = {Dimensional::Dimension::Time, 1e-3, Dimensional::Unit::hydro};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(contact_relaxation_time, default_qty));

	default_qty = {Dimensional::Dimension::Time, 0, Dimensional::Unit::hydro};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(contact_relaxation_time_tan, default_qty));
	
	default_qty = {Dimensional::Dimension::Force, 0.1, Dimensional::Unit::kn};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(min_kn_auto_det, default_qty));

	default_qty = {Dimensional::Dimension::Force, 1e3, Dimensional::Unit::kn};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(max_kn_auto_det, default_qty));

	default_qty ={Dimensional::Dimension::Force, 0.1, Dimensional::Unit::kt};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(min_kt_auto_det, default_qty));

	default_qty = {Dimensional::Dimension::Force, 1e3, Dimensional::Unit::kt};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(max_kt_auto_det, default_qty));

	default_qty = {Dimensional::Dimension::Time, 0, Dimensional::Unit::hydro};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(TA_adhesion.activation_time, default_qty));

	/*----------  True dim vals  ----------*/
	
	default_qty = {Dimensional::Dimension::TimeOrStrain, 10, Dimensional::Unit::hydro};
	TrueDimValDblParams.push_back(PARAM_INIT(time_end, default_qty));

	default_qty = {Dimensional::Dimension::TimeOrStrain, 1e-2, Dimensional::Unit::hydro};
	TrueDimValDblParams.push_back(PARAM_INIT(output.time_interval_output_data, default_qty));

	default_qty = {Dimensional::Dimension::TimeOrStrain, 1e-1, Dimensional::Unit::hydro};
	TrueDimValDblParams.push_back(PARAM_INIT(output.time_interval_output_config, default_qty));

	default_qty = {Dimensional::Dimension::TimeOrStrain, 1e-3, Dimensional::Unit::hydro};
	TrueDimValDblParams.push_back(PARAM_INIT(output.initial_log_time, default_qty));
}

void ParameterSetFactory::setFromFile(const std::string& filename_parameters)
{
	/**
	 \brief Read and parse the parameter file
	 */
	std::ifstream fin;
	fin.open(filename_parameters.c_str());
	if (!fin) {
		std::ostringstream error_str;
		error_str << " Parameter file '" << filename_parameters << "' not found." << std::endl;
		throw std::runtime_error(error_str.str());
	}
	std::string keyword, value;
	while (!fin.eof()) {
		std::string line;
		if (!getline(fin, line, ';')) {
			break;
		}
		if (fin.eof()) {
			break;
		}
		std::string str_parameter;
		removeBlank(line);
		str_parameter = line;
		std::string::size_type begin_comment;
		std::string::size_type end_comment;
		do {
			begin_comment = str_parameter.find("/*");
			end_comment = str_parameter.find("*/");
			if (begin_comment > 10000) {
				break;
			}
			str_parameter = str_parameter.substr(end_comment+2);
		} while (true);
		if (begin_comment > end_comment) {
			std::cerr << str_parameter.find("/*") << std::endl;
			std::cerr << str_parameter.find("*/") << std::endl;
			throw std::runtime_error("syntax error in the parameter file.");
		}
		std::string::size_type pos_slashslash = str_parameter.find("//");
		if (pos_slashslash != std::string::npos) {
			throw std::runtime_error(" // is not the syntax to comment out. Use /* comment */");
		}
		Str2KeyValue(str_parameter, keyword, value);
		setParameterFromKeyValue(keyword, value);
	}
	fin.close();
}

void ParameterSetFactory::setParameterFromKeyValue(const std::string &keyword, 
												   const std::string &value)
{
	for (auto &inp: BoolParams) {
		if (inp.name_str == keyword) {
			inp.value = str2bool(value);
			return;
		}
	}
	for (auto &inp: DoubleParams) {
		if (inp.name_str == keyword) {
			inp.value = stod(value);
			return;
		}
	}
	for (auto &inp: IntParams) {
		if (inp.name_str == keyword) {
			inp.value = stoi(value);
			return;
		}
	}
	for (auto &inp: StrParams) {
		if (inp.name_str == keyword) {
			inp.value = value;
			inp.value.erase(remove(inp.value.begin(), inp.value.end(), '\"' ), inp.value.end());
			return;
		}
	}
	for (auto &inp: DimValDblParams) {
		if (inp.name_str == keyword) {
			inp.value = value;
			return;
		}
	}
	for (auto &inp: TrueDimValDblParams) {
		if (inp.name_str == keyword) {
			inp.value = value;
			return;
		}
	}
	for (auto &inp: ForceScaleParams) {
		if (inp.name_str == keyword) {
			inp.value.dim_qty = value;
			return;
		}
	}
	std::ostringstream error_str;
	error_str << "keyword " << keyword << " is not associated with an parameter" << std::endl;
	throw std::runtime_error(error_str.str());
}

std::vector<Dimensional::ForceScale> ParameterSetFactory::getForceScales() const
{
	std::vector<Dimensional::ForceScale> fs;
	for (auto in_fs: ForceScaleParams) {
		fs.push_back(in_fs.value);
	}
	return fs;
}

void ParameterSetFactory::convertParameterUnit(const Dimensional::UnitSystem &unit_system, 
											   InputParameter<Dimensional::DimensionalQty<double>> &param)
{
	if (param.value.dimension == Dimensional::Dimension::TimeOrStrain) {
		if (param.value.unit == Dimensional::Unit::hydro) {
			param.value.dimension = Dimensional::Dimension::Strain;
			param.value.unit = Dimensional::Unit::none;
		} else {
			param.value.dimension = Dimensional::Dimension::Time;
		}
	}
	if (param.value.value == 0) {
		param.value.unit = unit_system.getInternalUnit();
	} else {
		unit_system.convertToInternalUnit(param.value);
	}
}

void ParameterSetFactory::setSystemOfUnits(const Dimensional::UnitSystem &unit_system)
{
	auto force_scales = unit_system.getForceScales();
	for (auto &force: ForceScaleParams) {
		if (force_scales.count(force.value.type)) {
			force.value.dim_qty = force_scales.at(force.value.type);
		}
	}
	for (auto &param: DimValDblParams) {
		convertParameterUnit(unit_system, param);
	}
	for (auto &param: TrueDimValDblParams) {
		convertParameterUnit(unit_system, param);
	}
}

ParameterSet ParameterSetFactory::getParameterSet() const
{
	ParameterSet p;
	for (auto &inp: BoolParams) {
		inp.exportToParameterSet(p, inp);
	}
	for (auto &inp: DoubleParams) {
		inp.exportToParameterSet(p, inp);
	}
	for (auto &inp: IntParams) {
		inp.exportToParameterSet(p, inp);
	}
	for (auto &inp: StrParams) {
		inp.exportToParameterSet(p, inp);
	}
	for (auto &inp: DimValDblParams) {
		inp.exportToParameterSet(p, inp);
	}
	for (auto &inp: TrueDimValDblParams) {
		inp.exportToParameterSet(p, inp);
	}
	for (auto &inp: ForceScaleParams) {
		inp.exportToParameterSet(p, inp);
	}
	return p;
}

} // namespace Parameters
