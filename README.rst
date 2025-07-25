KIM Validation and Verification
===============================

|Testing| |PyPI|

.. |Testing| image:: https://github.com/openkim/kimvv/actions/workflows/test.yml/badge.svg
   :target: https://github.com/openkim/kimvv/actions/workflows/test.yml
.. |PyPI| image:: https://img.shields.io/pypi/v/kimvv.svg
   :target: https://pypi.org/project/kimvv/

This package allows the user to run any `OpenKIM <https://openkim.org/>`_ Test Drivers written using the `kim-tools <https://kim-tools.readthedocs.io>`_ package locally. A "Test Driver" is
a computational protocol that reports one or more material properties using the `KIM Properties Framework <https://openkim.org/doc/schema/properties-framework/>`_

List of included Test Drivers:

  * EquilibriumCrystalStructure
  * ElasticConstantsCrystal

Currently, all Test Drivers require the AFLOW software to be installed and in your PATH. See https://kim-tools.readthedocs.io/en/stable/#doc-standalone-installation for installation info.

Basic usage example:
--------------------

.. code-block:: python

    from kimvv import EquilibriumCrystalStructure, ElasticConstantsCrystal
    from ase.build import bulk
    from json import dumps

    # If a string is passed when instantiating the class, it is assumed to be a KIM model name
    relax = EquilibriumCrystalStructure('LennardJones_Ar')

    # Every Test Driver is able to take an Atoms object
    relax(bulk('Ar','fcc',5.0))

    # Access the list of dictionaries containing the material properties reported by the Test Driver
    print(dumps(relax.property_instances,indent=2))

    # All Test Drivers besides EquilibriumCrystalStructure expect to be
    # passed a relaxed structure. This can be either a relaxed Atoms
    # object, or a results dictionary from an EquilibriumCrystalStructure
    # run. Any element of the list returned by EquilibriumCrystalStructure
    # will do, as they all contain a description of the crystal structure
    elastic = ElasticConstantsCrystal('LennardJones_Ar')
    elastic(relax.property_instances[0])
    print(dumps(elastic.property_instances,indent=2))

    # You can also use a generic ASE calculator (as long as the Test Driver only uses ASE for calculations,
    # i.e. this will not work for Test Drivers that do MD using LAMMPS)
    # In this case you don't even need kimpy or the KIM API installed.
    from ase.calculators.lj import LennardJones
    relax = EquilibriumCrystalStructure(LennardJones(sigma=3.4,epsilon=0.0104,rc=8.15))
    relax(bulk('Ar','fcc',5.0))


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
