import os
import pathlib
from typing import Any

import kim_edn

import kimvv 

class KIMVVTestDriver:
    @property
    def kimspec(self):
        mypath = pathlib.Path(__file__).parent.resolve()
        myname = self.__class__.__name__
        return kim_edn.load(os.path.join(mypath, myname, "kimspec.edn"))

    def _resolve_dependencies(self, **kwargs):
        '''
        defaults to equilibrium but can be defined within 
        TestDriver for driver specific needs
        '''
        # updates kwargs if needed and returns
        if isinstance(self, kimvv.EquilibriumCrystalStructure):
            return kwargs
        else:
            print ("Resolving dependencies...")
            ecs_test = kimvv.EquilibriumCrystalStructure(self._calc)
            ecs_test(self._get_atoms())
            # is this fragile?
            self._update_nominal_parameter_values(ecs_test._get_atoms())
            return kwargs
          
    @classmethod
    def printdoc(cls):
        print(cls._calculate.__doc__)

# new call decorator
def override_call_method(cls):
    def __call__(self, material: Any = None, **kwargs):
        """
        Taken from kim-tools with added dependency functionality
        Main operation of a Test Driver:

            * Run :func:`~KIMTestDriver._setup` (the base class provides a barebones
              version, derived classes may override)
            * Call :func:`~KIMTestDriver._calculate` (implemented by each individual
              Test Driver)

        Args:
            material:
                Placeholder object for arguments describing the material to run
                the Test Driver on

        Returns:
            The property instances calculated during the current run
        """
        # count how many instances we had before we started
        previous_properties_end = len(self.property_instances)

        os.makedirs("output", exist_ok=True)

        # _setup is likely overridden by an derived class
        self._setup(material, **kwargs)
        # resolve dependencies
        # since input to calculate may depend on output, return kwargs
        kwargs = self._resolve_dependencies(**kwargs)

        # implemented by each individual Test Driver
        self._calculate(**kwargs)

        # The current invocation returns a Python list of dictionaries containing all
        # properties computed during this run
        return self.property_instances[previous_properties_end:]
        
    cls.__call__ = __call__
    return cls
