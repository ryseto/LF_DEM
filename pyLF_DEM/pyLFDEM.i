%module pyLFDEM
%include <typemaps.i>
%apply double &INPUT { double &next_output_data };
%apply double &INPUT { double &next_output_config };
%apply int &INPUT { int &binconf_counter };
%{
/* Put headers and other declarations here */
#define SWIG_FILE_WITH_INIT
#include "../LF_DEM/Simulation.h"
#include "../LF_DEM/SystemHelperFunctions.h"
#include "../LF_DEM/Timer.h"
#include "../LF_DEM/global.h"
#include "../LF_DEM/ParameterSet.h"
#include "../LF_DEM/ParticleConfig.h"
%}
%include <std_string.i>
%include <std_set.i>
%include <std_map.i>
%include <std_vector.i>
%include <std_pair.i>
%include <std_shared_ptr.i>
%shared_ptr(ParticleConfig);
%rename(TAAParams) Interactions::TActAdhesion::Parameters;
%rename(ConfinementParams) Confinement::Parameters;

%include "../LF_DEM/vec3d.h"
%include "../LF_DEM/Sym2Tensor.h"
%include "../LF_DEM/DimensionalQty.h"
%include "../LF_DEM/Configuration.h"
%include "../LF_DEM/ParameterSet.h"
%include "../LF_DEM/ImposedDeformation.h"
%include "../LF_DEM/ParticleConfig.h"
%include "../LF_DEM/System.h"
%include "../LF_DEM/ContactParams.h"
%include "../LF_DEM/LubricationParams.h"
%include "../LF_DEM/RepulsiveForceParams.h"
%include "../LF_DEM/ConfinementParams.h"
%include "../LF_DEM/VanDerWaalsParams.h"
%include "../LF_DEM/DimerParams.h"
%include "../LF_DEM/DimerState.h"
%include "../LF_DEM/TimeActivatedAdhesion_Params.h"
%include "../LF_DEM/ParameterSetFactory.h"
%include "../LF_DEM/SystemHelperFunctions.h"
%include "../LF_DEM/ShearType.h"
%include "../LF_DEM/ControlVariable.h"
%include "../LF_DEM/Simulation.h"
%include "../LF_DEM/global.h"
%include "../LF_DEM/Timer.h"

// Instantiate templates used by example
namespace std {
  %template(StringVector) vector<string>;
  %template(Vec3dVector) vector<vec3d>;
  %template(DoubleVector) vector<double>;
  %template(IntVector) vector<int>;
  %template(Sym2Vector) vector<Sym2Tensor>;
  %template(UnloadedDimerStateVector) vector<Interactions::Dimer::UnloadedDimerState>;
  %template(PairUIntUInt) pair<unsigned int, unsigned int>;
  %template(PairIntInt) pair<int, int>;
  %template(PairDoubleString) pair<double, string>;
  %template(SetString) set<string>;
  %template(StringStringMap) map<string,string>;
}
%template(DoubleDimQty) Dimensional::DimensionalQty<double>;
// imperfect error handling. Errors get passed to Python, but sometimes still seg fault
%include "exception.i"
%exception {
    try {
        $action
    } catch(const std::exception& e) {
	  SWIG_exception(SWIG_RuntimeError, e.what());
    }
}
