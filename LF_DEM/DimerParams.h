#ifndef __LF_DEM__DimerParams__
#define __LF_DEM__DimerParams__

namespace Interactions
{
namespace Dimer
{

struct DimerParams {
	double stiffness;  ///< Spring stiffness. [0]
	double relaxation_time; ///< Spring-Dashpot relaxation time. [0]
};

}
}
#endif