KIM Validation and Verification
===============================

|Testing| |PyPI|

.. |Testing| image:: https://github.com/openkim/kimvv/actions/workflows/test.yml/badge.svg
   :target: https://github.com/openkim/kimvv/actions/workflows/test.yml
.. |PyPI| image:: https://img.shields.io/pypi/v/kimvv.svg
   :target: https://pypi.org/project/kimvv/

This package allows the user to run any Test Drivers written using the `kim-tools <https://kim-tools.readthedocs.io>`_ package locally. Currently, all Test Drivers require the AFLOW software to be installed and in your PATH. See https://kim-tools.readthedocs.io/en/stable/#doc-standalone-installation for installation info.

Basic usage example:
--------------------

.. code-block:: python

    from kimvv import EquilibriumCrystalStructure
    from ase.build import bulk
    from json import dumps

    # If a string is passed when instantiating the class, it is assumed to be a KIM model name
    relax = EquilibriumCrystalStructure('LennardJones_Ar')

    # Pass an Atoms object
    relax(bulk('Ar','fcc',5.0))

    # Access the results dictionary
    print(dumps(relax.property_instances,indent=2))

    # You can also use a generic ASE calculator (as long as the Test Driver doesn't use external simulation codes)
    # In this case you don't even need kimpy or the KIM API installed.

    from ase.calculators.lj import LennardJones
    relax = EquilibriumCrystalStructure(LennardJones(sigma=3.4,epsilon=0.0104,rc=8.15))
    relax(bulk('Ar','fcc',5.0), optimize=True)


Usage example 2
---------------
Querying for all DFT-relaxed structures for a given combination of elements in OpenKIM and relaxing them with your potential

.. code-block:: python

    from kimvv import EquilibriumCrystalStructure
    from kim_tools import (
      query_crystal_structures,
      get_deduplicated_property_instances
    )
    from json import dumps
    from ase.calculators.lj import LennardJones

    # Query for all relaxed Argon reference data in OpenKIM
    raw_structs = query_crystal_structures(stoichiometric_species=["Ar"])

    # Deduplicate them
    unique_structs = get_deduplicated_property_instances(raw_structs, allow_rotation=True)

    # Instantiate the Driver with your model
    relax = EquilibriumCrystalStructure(LennardJones(sigma=3.4,epsilon=0.0104,rc=8.15))

    # Run the Driver with each structure. As this is run, the driver internally accumulates
    # Property Instances
    for struct in unique_structs:
      relax(struct)

    # Access the results as a dictionary. For each structure, there are 3 properties
    # (structure, binding energy, density)
    print(dumps(relax.property_instances,indent=2))
