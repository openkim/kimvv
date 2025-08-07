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

    def _resolve_dependencies(self, material, **kwargs):
        """
        defaults to equilibrium but can be defined within
        TestDriver for driver specific needs
        """
        # updates kwargs if needed and returns
        if isinstance(self, kimvv.EquilibriumCrystalStructure):
            return material, kwargs
        else:
            print("Resolving dependencies...")
            ecs_test = kimvv.EquilibriumCrystalStructure(self._calc)
            ecs_results = ecs_test(material)
            for result in ecs_results:
                if result["property-id"].endswith("crystal-structure-npt"):
                    material_relaxed = result
                    break
            return material_relaxed, kwargs

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

        # resolve dependencies
        # since input to calculate may depend on output, return kwargs
        material_relaxed, kwargs = self._resolve_dependencies(material, **kwargs)

        self._setup(material_relaxed, **kwargs)

        # implemented by each individual Test Driver
        self._calculate(**kwargs)

        # The current invocation returns a Python list of dictionaries containing all
        # properties computed during this run
        return self.property_instances[previous_properties_end:]

    cls.__call__ = __call__
    return cls
