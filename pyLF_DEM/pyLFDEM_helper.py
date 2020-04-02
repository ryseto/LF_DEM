import pyLFDEM as lf
from collections import namedtuple

def getLFControlVar(control_var):
	if control_var == "r" or control_var == "rate":
		return lf.ControlVariable_rate
	if control_var == "s" or control_var == "stress":
		return lf.ControlVariable_stress

	raise RuntimeError("Unknown control variable ")
		# +control var+"\nPossible choices are \"r\"/\"rate\" or \"s\"/\"stress\"")



SimuArgs = namedtuple('SimuArgs', ['ctrl_str', 'conf_file', 'params_file', 'binary_conf', 'identifier'])

def setupBasicSimu(**kwargs):
	"""
		Inputs:
		  	**kwargs :
		  		Accepts following keys: 
		  			- 'ctrl_str': string to specify the control variable and its value, with form "X Yz", 
		  						  where "X"="r" for rate control, "s" for stress control,
		  						  "Y" a float and "z" a unit suffix.
		  						  Example: "r 1.5kn"
		  						  (required)

		  			- 'conf_file': name of file containing the initial configuration, either in text or binary format (required)

		  			- 'binary_conf': boolean, True if the file passed in 'conf_file' is binary (default False)

		  			- 'params_file': name of parameter file (required)

		  			- 'identifier': optional name of the simulation (default "pyLFDEM")	
	"""
	try:
		ctrl_str = kwargs.pop('ctrl_str')
	except KeyError:
		raise KeyError("Please provide 'ctrl_str'")
	try:
		conf_file = kwargs.pop('conf_file')
	except KeyError:
		raise KeyError("Please provide 'conf_file'")
	try:
		params_file = kwargs.pop('params_file')
	except KeyError:
		raise KeyError("Please provide 'params_file'")		
	try:
		identifier = kwargs.pop('identifier')
	except KeyError:
		identifier = "pyLFDEM"
	try:
		binary_conf = bool(kwargs.pop('binary_conf'))
	except KeyError:
		binary_conf = False

	simu = lf.Simulation()
	simu.setupControl(getLFControlVar(ctrl_str.split()[0]), ctrl_str.split()[1])

	input_files = lf.StringStringMap({"config": conf_file, "params": params_file})
	simu.setupSimulation(str(kwargs), input_files, binary_conf, identifier)

	return simu