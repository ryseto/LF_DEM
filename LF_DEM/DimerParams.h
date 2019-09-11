#ifndef __LF_DEM__DimerParams__
#define __LF_DEM__DimerParams__

namespace Interactions
{
namespace Dimer
{

struct DimerParams {
	double stiffness;  ///< Spring stiffness. [0]
	double resistance; ///< Dashpot damping coefficient. [0]
};

}
}
#endif