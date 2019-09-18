#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "ParameterSetFactory.h"
#define PARAM_INIT(name, default_value) {#name,  [](ParameterSet &p, InputParameter<decltype(ParameterSet::name)> in) {p.name = in.value;}, default_value}
#define PARAM_INIT_ENUMCLASS(name, default_value) {#name,\
												   [](ParameterSet &p, InputParameter<std::underlying_type<decltype(ParameterSet::name)>::type> in) {p.name = static_cast<decltype(ParameterSet::name)>(in.value);},\
													static_cast<std::underlying_type<decltype(ParameterSet::name)>::type>(default_value)}

#define PARAM_INIT_DIMQTY(name, default_value) {#name,  [](ParameterSet &p, InputParameter<Dimensional::DimensionalQty<decltype(ParameterSet::name)>> in) {p.name = in.value.value;}, default_value}
#define PARAM_INIT_FORCESCALE(name, default_value) {#name,  [](ParameterSet &p, InputParameter<Dimensional::ForceScale> in) {p.name = in.value.dim_qty.value;}, default_value}

inline void removeBlank(std::string& str)
{
	str.erase(std::remove_if(str.begin(), str.end(), (int(*)(int))isspace), str.end());
}

bool str2bool(const std::string& value)
{
	if (value == "true") {
		return true;
	} else if (value == "false") {
		return false;
	} else {
		std::cerr << "The value should be true or false" << std::endl;
		exit(1);
	}
}

void Str2KeyValue(const std::string& str_parameter,
				  std::string& keyword,
				  std::string& value)
{
	std::string::size_type pos_equal = str_parameter.find("=");
	keyword = str_parameter.substr(0, pos_equal);
	value = str_parameter.substr(pos_equal+1);
}

namespace Parameters {

ParameterSetFactory::ParameterSetFactory(Dimensional::Unit guarranted_unit) 
{
	setDefaultValues(guarranted_unit);
}

void ParameterSetFactory::setDefaultValues(Dimensional::Unit guarranted_unit)
{

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
		PARAM_INIT(output.new_material_functions, false),
		PARAM_INIT(output.out_data_vel_components, false),
		PARAM_INIT(output.out_binary_conf, false),
		PARAM_INIT(output.out_data_interaction, true),
		PARAM_INIT(output.out_data_particle, true),
		PARAM_INIT(output.out_gsd, true),
		PARAM_INIT(output.out_gsd_na_velocity, true),
		PARAM_INIT(output.origin_zero_flow, true),
		PARAM_INIT(output.log_time_interval, false),
		PARAM_INIT(output.out_na_vel, false),
		PARAM_INIT(output.out_na_disp, false),
		PARAM_INIT(output.effective_coordination_number, false),
		PARAM_INIT(check_static_force_balance, false),
		PARAM_INIT(solvent_flow, false),
		PARAM_INIT(confinement.on, false),
		PARAM_INIT(lub.smooth, false)
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
		PARAM_INIT(lub.max_gap, 0.5),
		PARAM_INIT(lub.regularization_length, 1e-3),
		PARAM_INIT(sd_coeff, 1),
		PARAM_INIT(repulsion.screening_length, 0.05),
		PARAM_INIT(repulsion.max_length, -1),
		PARAM_INIT(repulsion.smoothing, 5e-3),
		PARAM_INIT(vdw.coefficient, -1),
		PARAM_INIT(vdw.singularity_cutoff, 0.1),
		PARAM_INIT(vdw.amplitude, 0),
		PARAM_INIT(magic_angle, 0),
		PARAM_INIT(contact.mu_static, 1),
		PARAM_INIT(contact.mu_dynamic, -1),
		PARAM_INIT(contact.mu_rolling, 0),
		PARAM_INIT(TA_adhesion.adhesion_range, 1e-2),
		PARAM_INIT(output.recording_start, 1),
		PARAM_INIT(sj_disp_max_shrink_factor, 1.1),
		PARAM_INIT(sj_disp_max_goal, 1e-6),
		PARAM_INIT(sj_shear_rate, 0),
		PARAM_INIT(sj_velocity, 1e-3),
		PARAM_INIT(body_force_angle, 0),
		PARAM_INIT(sflow_dx, 5),
		PARAM_INIT(sflow_smooth_length, 3),
		PARAM_INIT(sflow_ReNum_p, 0.001),
		//PARAM_INIT(sflow_ReNum, 0.1),
		PARAM_INIT(sflow_Darcy_coeff, 1),
		PARAM_INIT(sflow_pcontrol_increment, 1e-4),
		PARAM_INIT(sflow_pcontrol_rtime, 0.1),
		PARAM_INIT(sflow_pcontrol_damper, 100),
		PARAM_INIT(sflow_target_flux, 0),
		PARAM_INIT(confinement.y_min, 0),
		PARAM_INIT(confinement.y_max, 0),
	};

	/*================================
	=            INTEGERS            =
	================================*/
	IntParams = \
	{
		PARAM_INIT(output.nb_output_data_log_time, 100),
		PARAM_INIT(output.nb_output_config_log_time, 100),
		PARAM_INIT(integration_method, 1),
		PARAM_INIT(np_fixed, 0),
		PARAM_INIT(simulation_mode, 0),
		PARAM_INIT(sj_check_count, 500),
		PARAM_INIT(sj_reversal_repetition, 10),
		PARAM_INIT(sflow_boundary_conditions, 0),
		PARAM_INIT(sflow_Darcy_power, 0)
	};

	UIntParams = \
	{
		PARAM_INIT_ENUMCLASS(contact.friction_model, Interactions::FrictionModel::Coulomb)
	};

	/*===============================
	=            STRINGS            =
	===============================*/
	StrParams = \
	{
		PARAM_INIT(flow_type, "shear"),
		PARAM_INIT(event_handler, ""),
		PARAM_INIT(output.out_particle_stress, ""),
		PARAM_INIT(lub.model, "none"),
		PARAM_INIT(sj_program_file, "")
	};

	/*====================================
	=            FORCE SCALES            =
	====================================*/
	Dimensional::ForceScale default_val;

	default_val = {Dimensional::Unit::kn, {Dimensional::Dimension::Force, 0, guarranted_unit}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(contact.kn, default_val));

	default_val = {Dimensional::Unit::kt, {Dimensional::Dimension::Force, 0, Dimensional::Unit::kn}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(contact.kt, default_val));

	default_val = {Dimensional::Unit::kr, {Dimensional::Dimension::Force, 0, Dimensional::Unit::kn}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(contact.kr, default_val));

	default_val = {Dimensional::Unit::kn, {Dimensional::Dimension::Force, 0, guarranted_unit}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(contact.ft_max, default_val));

	default_val = {Dimensional::Unit::critical_load, {Dimensional::Dimension::Force, 0, guarranted_unit}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(contact.critical_load, default_val));

	default_val = {Dimensional::Unit::adhesion, {Dimensional::Dimension::Force, 0, guarranted_unit}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(contact.adhesion, default_val));

	default_val = {Dimensional::Unit::delayed_adhesion, {Dimensional::Dimension::Force, 0, guarranted_unit}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(TA_adhesion.adhesion_max_force, default_val));

	default_val = {Dimensional::Unit::repulsion, {Dimensional::Dimension::Force, 0, guarranted_unit}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(repulsion.repulsion, default_val));

	default_val = {Dimensional::Unit::brownian, {Dimensional::Dimension::Force, 0, guarranted_unit}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(brownian, default_val));

	default_val = {Dimensional::Unit::bodyforce, {Dimensional::Dimension::Force, 0, guarranted_unit}};
	ForceScaleParams.push_back(PARAM_INIT_FORCESCALE(bodyforce, default_val));
	
	
	/*==============================================
	=            Dimensional Quantities            =
	==============================================*/

	/*----------  Double dim vals  ----------*/

	Dimensional::DimensionalQty<double> default_qty;

	default_qty = {Dimensional::Dimension::Time, 1e-4, guarranted_unit};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(dt, default_qty));

	default_qty = {Dimensional::Dimension::Time, -1, guarranted_unit};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(dt_max, default_qty));
	
	default_qty = {Dimensional::Dimension::Time, -1, guarranted_unit};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(dt_jamming, default_qty));
	
	default_qty = {Dimensional::Dimension::Time, 1e-3, guarranted_unit};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(contact.relaxation_time, default_qty));

	default_qty = {Dimensional::Dimension::Time, -1, guarranted_unit};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(contact.relaxation_time_tan, default_qty));
	
	default_qty = {Dimensional::Dimension::Force, 0.1, Dimensional::Unit::kn};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(min_kn_auto_det, default_qty));

	default_qty = {Dimensional::Dimension::Force, 1e3, Dimensional::Unit::kn};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(max_kn_auto_det, default_qty));

	default_qty ={Dimensional::Dimension::Force, 0.1, Dimensional::Unit::kt};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(min_kt_auto_det, default_qty));

	default_qty = {Dimensional::Dimension::Force, 1e3, Dimensional::Unit::kt};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(max_kt_auto_det, default_qty));

	default_qty = {Dimensional::Dimension::Time, 0, guarranted_unit};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(TA_adhesion.activation_time, default_qty));

	default_qty = {Dimensional::Dimension::Force, 0, guarranted_unit};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(dimer.stiffness, default_qty));

	default_qty = {Dimensional::Dimension::Time, 0, guarranted_unit};
	DimValDblParams.push_back(PARAM_INIT_DIMQTY(dimer.relaxation_time, default_qty));


	/*----------  True dim vals  ----------*/
	
	default_qty = {Dimensional::Dimension::TimeOrStrain, 10, guarranted_unit};
	TrueDimValDblParams.push_back(PARAM_INIT(time_end, default_qty));

	default_qty = {Dimensional::Dimension::TimeOrStrain, 1e-2, guarranted_unit};
	TrueDimValDblParams.push_back(PARAM_INIT(output.time_interval_output_data, default_qty));

	default_qty = {Dimensional::Dimension::TimeOrStrain, 1e-1, guarranted_unit};
	TrueDimValDblParams.push_back(PARAM_INIT(output.time_interval_output_config, default_qty));

	default_qty = {Dimensional::Dimension::TimeOrStrain, 1e-3, guarranted_unit};
	TrueDimValDblParams.push_back(PARAM_INIT(output.initial_log_time, default_qty));

	default_qty = {Dimensional::Dimension::Force, 0, guarranted_unit};
	TrueDimValDblParams.push_back(PARAM_INIT(confinement.k, default_qty));

}

void ParameterSetFactory::setFromFile(const std::string& filename_parameters)
{
	/**
	 \brief Read and parse the parameter file
	 */
	std::string indent = "  ParameterSetFactory::\t";
	std::cout << indent << "setFromFile..." << std::endl;
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
	std::cout << indent << "setFromFile...done" << std::endl;
	fin.close();
}

void ParameterSetFactory::setFromStringStream(std::stringstream& ss_initial_setup)
{
	std::string indent = "  ParameterSetFactory::\t";
	std::cout << indent << "setFromStringStream..." << std::endl;
	std::string keyword, value;
	std::string line;
	while (std::getline(ss_initial_setup, line, ';')) {
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
	std::cout << indent << "setFromStringStream...done" << std::endl;
}

void ParameterSetFactory::setFromLine(std::string& line)
{
	std::string indent = "  ParameterSetFactory::\t";
	std::cout << indent << "setFromStringStream..." << std::endl;
	removeBlank(line);
	std::string keyword, value;
	Str2KeyValue(line, keyword, value);
	std::cerr << keyword << ' ' << value << std::endl;
	setParameterFromKeyValue(keyword, value);
	std::cout << indent << "setFromStringStream...done" << std::endl;
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
	for (auto &inp: UIntParams) {
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
	error_str << "keyword " << keyword << " is not associated with any parameter" << std::endl;
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
	for (auto &inp: UIntParams) {
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

	Interactions::Lub::setupLubricationParameters(p.lub);
	Interactions::setupContactParameters(p.contact);
	return p;
}


} // namespace Parameters
