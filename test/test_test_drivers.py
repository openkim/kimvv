import json
import warnings
from typing import Union

import pytest
from ase.build import bulk
from ase.calculators.calculator import Calculator
from ase.calculators.lj import LennardJones
from kim_tools import KIMTestDriver, get_atoms_from_crystal_structure

import kimvv

with open("test_inputs.json") as f:
    DRIVERS = json.load(f)

# Test on FCC Au
MODELS = [
    "LennardJones612_UniversalShifted__MO_959249795837_003",
    "Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu",
    LennardJones(sigma=2.42324, epsilon=2.30580, rc=9.69298),
]

test_tuples = [(driver, model) for driver in DRIVERS for model in MODELS]


@pytest.mark.parametrize("td_name,model", test_tuples)
@pytest.mark.filterwarnings("ignore:(?!WARNING Your )")
def test_test_driver(td_name: str, model: Union[str, Calculator]) -> None:
    """
    Run ``td_name`` with ``model`` on and confirm that it returns at least one
    instance of the properties it claims to return, and no others.

    Args:
        td_name:
            The name of the class of Test Driver to run.
        model:
            The model to use.
    """
    # Start with FCC Au
    atoms_init = bulk("Au")

    # Relax it with EquilibriumCrystalStructure
    crystal_structure_relaxed = kimvv.EquilibriumCrystalStructure(model)(atoms_init)[0]

    TestDriver = getattr(kimvv, td_name)

    td_kwargs = DRIVERS[td_name]

    td = TestDriver(model)

    try:
        # Run the driver with crystal structure dict
        results_from_dict = td(crystal_structure_relaxed, **td_kwargs)

        # Run the atoms with relaxed atoms
        atoms_relaxed = get_atoms_from_crystal_structure(crystal_structure_relaxed)
        results_from_atoms = td(atoms_relaxed, **td_kwargs)

        for results in (results_from_dict, results_from_atoms):
            # Should return at least something
            assert len(results_from_atoms) > 1

            # If we have properties in our kimspec, check that the test driver
            # only reports those
            if "properties" in td.kimspec:
                properties = td.kimspec["properties"]
                for result in results:
                    assert result["property-id"] in properties
            else:
                warnings.warn(
                    "WARNING Your kimspec.edn did not contain a 'properties' key. "
                    "this is acceptable, but I cannot test thhat the reported "
                    "properties are correct."
                )
    except KIMTestDriver.NonKIMModelError:
        warnings.warn(
            "WARNING Your Test Driver is unable to run with non-KIM calculators "
            "because it requests the kim_model_name property. This is acceptable, "
            "but only if strictly necessary (e.g. needing a KIM model name to run a "
            "LAMMPS calculation)"
        )
