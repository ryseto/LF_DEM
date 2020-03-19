#ifndef __LF_DEM__SystemHF__
#define __LF_DEM__SystemHF__

double evaluateMaxAngVelocity(const System & sys);
double evaluateMinGap(const System & sys);
double evaluateAvgContactGap(const System & sys);
double evaluateMaxContactGap(const System & sys);
double evaluateMaxVelocity(const System &sys);
double evaluateMaxNAVelocityComponent(const System & sys, std::string component);
double evaluateMaxDispTan(const System &sys);
double evaluateMaxDispRolling(const System &sys);
std::pair<unsigned int, unsigned int> countNumberOfContact(const System &sys);
// std::pair<unsigned, double> getTAAdhesionActivityStatistics(const System &sys);
std::pair<unsigned, double> getActAdhesionActivityStatistics(const System &sys);
double getPotentialEnergy(const System &sys);
void isInContact(const System &sys, std::vector<int> &isincontact);

#endif /* defined(__LF_DEM__SystemHF__) */
