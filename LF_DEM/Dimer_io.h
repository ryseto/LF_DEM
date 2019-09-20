#ifndef __LF_DEM__Dimer_io__
#define __LF_DEM__Dimer_io__

#include <string>
#include <iostream>
#include <vector>

namespace Interactions {

namespace Dimer {

class Dimer;
class DimerManager;

namespace io 
{
std::vector <struct DimerState> readStatesBStream(std::istream &input);
std::vector <struct UnloadedDimerState> readTxtDimer(const std::string& filename);
void writeStatesBStream(std::ostream &conf_export, const DimerManager &dm);
}

}
}

#endif