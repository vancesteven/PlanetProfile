% pylib = py.importlib.import_module('seafreeze');
% py.importlib.reload(pylib);

%pylib = py.importlib.import_module('Thermodynamics.FromLiterature.Helgerud_sI');
%py.importlib.reload(pylib);

%pylib = py.importlib.import_module('Thermodynamics.FromLiterature.getK_Andersson2005');
%py.importlib.reload(pylib);

pylib = py.importlib.import_module('Thermodynamics.FromLiterature.conductiveMantleTemperature');
py.importlib.reload(pylib);

%pylib = py.importlib.import_module('Thermodynamics.FromLiterature.ConvectionDeschampsSotin2001');
%py.importlib.reload(pylib);

pylib = py.importlib.import_module('MantlePlot');
py.importlib.reload(pylib);

pylib = py.importlib.import_module('MatToPy');
py.importlib.reload(pylib);