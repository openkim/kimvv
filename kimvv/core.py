import inspect
import os
import pathlib
import shutil
from tempfile import TemporaryDirectory
from typing import Any

import kim_edn

import kimvv


def get_default_args(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }


class KIMVVTestDriver:
    @classmethod
    def get_kimspec(cls):
        mypath = pathlib.Path(__file__).parent.resolve()
        myname = cls.__name__
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
        print("\nDescription of method and non-structure arguments:\n")
        print(cls._calculate.__doc__)
        print("\nDefaults:")
        print(get_default_args(cls._calculate))
        print()


# new call decorator
def override_call_method(cls):
    def __call__(self, material: Any = None, **kwargs):
        """
        Taken from kim-tools with added dependency functionality
        Main operation of a Test Driver:

            * Call :func:`~KIMTestDriver._init_output_dir`
            * Run :func:`~KIMTestDriver._setup` (the base class provides a barebones
              version, derived classes may override)
            * Call :func:`~KIMTestDriver._calculate` (implemented by each individual
              Test Driver)
            * Call :func:`~KIMTestDriver._archive_aux_files`

        Args:
            material:
                Placeholder object for arguments describing the material to run
                the Test Driver on

        Returns:
            The property instances calculated during the current run
        """
        # count how many instances we had before we started
        previous_properties_end = len(self.property_instances)

        # Set up the output directory
        self._init_output_dir()

        with TemporaryDirectory() as d:
            original_output_bak = shutil.move("output", d)
            try:
                # resolve dependencies
                # since input to calculate may depend on output, return kwargs
                material_relaxed, kwargs = self._resolve_dependencies(
                    material, **kwargs
                )
            finally:
                # Dependencies wrote an output directory
                if os.path.exists("output"):
                    i = 0
                    while os.path.exists(f"output.dependency.{i}"):
                        i += 1
                    dep_output_bak = f"output.dependency.{i}"
                    msg = (
                        "Backing up 'output' directory written by dependency to "
                        f"{dep_output_bak}"
                    )
                    print(msg)
                    os.rename("output", dep_output_bak)
                # restore original output
                shutil.move(original_output_bak, ".")

        try:
            # _setup is likely overridden by an derived class
            self._setup(material_relaxed, **kwargs)

            # implemented by each individual Test Driver
            self._calculate(**kwargs)
        finally:
            # Postprocess output directory for this invocation
            self._archive_aux_files()

        # The current invocation returns a Python list of dictionaries containing all
        # properties computed during this run
        return self.property_instances[previous_properties_end:]

    cls.__call__ = __call__
    return cls
