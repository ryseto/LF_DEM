struct InputValue{
  std::string type;
  std::string name; // the name of the ParameterSet member
  double *value; // a pointer to the actual ParameterSet member
  std::string unit;
};
