import pyLFDEM as lf
from collections import namedtuple

def getLFControlVar(control_var):
	if control_var == "r" or control_var == "rate":
		return lf.ControlVariable_rate
	if control_var == "s" or control_var == "stress":
		return lf.ControlVariable_stress

	raise RuntimeError("Unknown control variable ")
		# +control var+"\nPossible choices are \"r\"/\"rate\" or \"s\"/\"stress\"")



def setupBasicSimu(**kwargs):
	"""
		Wrapper around LF_DEM's Simulation::setupSimulation method.

		Inputs:
		  	**kwargs :
		  		Accepts following keys: 
		  			- 'ctrl_str': string to specify the control variable and its value, with form "X Yz", 
		  						  where "X"="r" for rate control, "s" for stress control,
		  						  "Y" a float and "z" a unit suffix.
		  						  Example: "r 1.5kn"
		  						  (required)

		  			- 'conf_file': string, name of file containing the initial configuration, either in text or binary format (required)

		  			- 'binary_conf': boolean, True if the file passed in 'conf_file' is binary (default False)

		  			- 'params_file': string, name of parameter file (required)

		  			- 'identifier': string, optional name of the simulation (default "pyLFDEM")	

		  			- 'call_str': string, call arguments to be printed by LF_DEM in the input_ file of the simulation (default "").
								   setupBasicSimu will append its own call signature to this string

					- 'overwrite': boolean, start simulation even if corresponding output files already present, 
						           and about to be overwritten (default False)

		Output:
			Instance of Simulation class configured according to inputs.

	"""
	in_args = "\nThrough pyLFDEM_helper.setupBasicSimu("+str(kwargs)+")"
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
	try:
		call_str = kwargs.pop('call_str')
	except KeyError:
		call_str = ""
	try:
		overwrite = kwargs.pop('overwrite')
	except KeyError:
		overwrite = False
	
	call_str = call_str + in_args

	simu = lf.Simulation()
	simu.force_to_run = overwrite
	simu.setupControl(getLFControlVar(ctrl_str.split()[0]), ctrl_str.split()[1])

	input_files = lf.StringStringMap({"config": conf_file, "params": params_file})
	simu.setupSimulation(call_str, input_files, binary_conf, identifier)

	return simu