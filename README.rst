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
  * CrystalStructureAndEnergyVsPressure
  * GroundStateCrystalStructure
  * VacancyFormationEnergyRelaxationVolumeCrystal

Installation
------------
``kimvv`` is installable with ``pip install kimvv``, but it requires some non-Python rquirements to be installed first.
These prerequisites are decribed in installation info for ``kim-tools``, the backend for ``kimvv``, here: https://kim-tools.readthedocs.io/en/stable/#doc-standalone-installation.

Basic usage example:
--------------------
Computing elastic constants for FCC argon using an example KIM potential

.. code-block:: python

  from kimvv import ElasticConstantsCrystal
  from ase.build import bulk
  from json import dumps

  # The Test Driver must be instantiated with an ASE Calculator object
  # or a string indicating a KIM model name
  elast = ElasticConstantsCrystal('LennardJones_Ar')

  # To perform the computation, call the Test Driver object. The first argument
  # to most Test Drivers is the crystal structure to perform the computation on.
  # To see an explanation of the calculation and a description of the
  # additonal arguments, use .printdoc()
  elast.printdoc()

  # For the sake of speed, let's compute the elastic constants with the
  # "stress-condensed" method, instead of the default robust computation loop.
  # The crystal structure can be specified as an Atoms object. Any dependencies
  # (e.g. relaxing the crystal structure with EquilibriumCrystalStructure) are
  # automatically run.
  atoms = bulk("Ar", "fcc", 5.0)
  results = elast(atoms, method="stress-condensed")

  # Each Test Driver computes a list of one or more dictionaries, each defining
  # a material property in the format specified by the KIM Properties Framework.
  # The name of the property is in the "property-id" key. See
  # https://openkim.org/properties for the definition of each property.
  print(dumps(results, indent=2))


Usage example 2
---------------
Getting the anisotropic pressure-volume curve of HCP Ag using a non-KIM ASE Calculator and saving
the output files

.. code-block:: python

  from kimvv import CrystalStructureAndEnergyVsPressure
  from ase.build import bulk
  from ase.calculators.emt import EMT

  # The Test Driver must be instantiated with an ASE Calculator object
  # or a string indicating a KIM model name
  scan = CrystalStructureAndEnergyVsPressure(EMT())

  # To perform the computation, call the Test Driver object. The first argument
  # to most Test Drivers is the crystal structure to perform the compuation on.
  # To see the additonal arguments, use .printdoc() to print the docstring
  scan.printdoc()

  # The default volume range of 0.25-4.0 will take a long time to scan. Let's
  # do a much smaller range
  results = scan(
      bulk("Ag", "hcp", 2.92), min_fractional_volume=0.98, max_fractional_volume=1.02
  )

  # In addition to accessing the results as a Python dictionary, you can save them to
  # a file in .edn format. This is especially useful if the Test Driver produces
  # auxiliary files, like the pressure scan does. All auxiliary files will be written
  # to the parent directory of the path you specified.
  scan.write_property_instances_to_file("scan_output/results.edn")

Usage example 3
---------------
This example is functionally identical to the previous example, except the crystal is specified by
passing a dictionary specifying the symmetry-reduced description of the crystal

.. code-block:: python

  from kimvv import CrystalStructureAndEnergyVsPressure
  from ase.calculators.emt import EMT

  scan = CrystalStructureAndEnergyVsPressure(EMT())

  # Specify the material using a symmetry-reduced dictionary.
  # Internally, all kimvv Test Drivers use this representation,
  # so this allows more direct control, as an Atoms object will
  # be converted to this regardless. This allows you to specify
  # the crystal in a specific orientation that will be maintained
  # Notionally, this should be an instance of the
  # `crystal-structure-npt` OpenKIM Property, but the exact schema
  # is not enforced. As long as the following fields are present,
  # it will work: "prototype-label.source-value",
  # "stoichiometric-species.source-value", "a.source-value",
  # and, if the crystal has any free parameters,
  # "parameter-values.source-value".
  # For an exact definition of these fields, see
  # https://openkim.org/properties/show/crystal-structure-npt
  # For more info about the AFLOW Prototype Designation,
  # see section B here: https://arxiv.org/pdf/2401.06875
  material = {
      "prototype-label": {"source-value": "A_hP2_194_c"},
      "stoichiometric-species": {"source-value": ["Ag"]},
      "a": {
          "source-value": 2.933,
          "source-unit": "angstrom",
      },
      "parameter-names": {"source-value": ["c/a"]},
      "parameter-values": {"source-value": [1.6373338]},
  }


  results = scan(material, min_fractional_volume=0.98, max_fractional_volume=1.02)
  scan.write_property_instances_to_file("scan_output/results.edn")

Usage example 4
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
  # You can narrow the query further by specifying more information
  # about the crystal, see
  # https://kim-tools.readthedocs.io/en/stable/kim_tools.test_driver.html#kim_tools.test_driver.core.query_crystal_structures
  raw_structs = query_crystal_structures(stoichiometric_species=["Ar"])

  # Deduplicate them
  unique_structs = get_deduplicated_property_instances(raw_structs, allow_rotation=True)

  # Instantiate the Driver with your model
  relax = EquilibriumCrystalStructure(LennardJones(sigma=3.4,epsilon=0.0104,rc=8.15))

  # Run the Driver with each structure. As this is run, the driver internally accumulates
  # Property Instances
  for struct in unique_structs:
    relax(struct)

  # In addition to returning the Property Instances for the current run, Test Drivers
  # accumulate all computed Property Instances. They can be accessed like this:
  print(dumps(relax.property_instances, indent=2))
