Simulating Photodynamic Inactivation of a Cocci bacterium 
------------------------------------------------------------------------

|PyPI version| |Actions Status| |Downloads| |License|

.. |PyPI version| image:: https://img.shields.io/pypi/v/pdipy.svg?logo=PyPI&logoColor=brightgreen
   :target: https://pypi.org/project/pdipy/
   :alt: PyPI version

.. |Actions Status| image:: https://github.com/freiburgermsu/pdipy/workflows/Test%20PDIpy/badge.svg
   :target: https://github.com/freiburgermsu/pdipy/actions
   :alt: Actions Status

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

.. |Downloads| image:: https://pepy.tech/badge/pdipy
   :target: https://pepy.tech/project/pdipy
   :alt: Downloads

Antibiotic resistance is developing medical crisis that is projected to surpass cancer in annual deaths by mid-21st century. Photodynamic Inactivation (PDI) is a promising treatment method that escapes resistance evolution and may be an essential technology to hamper the growing threat of resistant pathogens; however, the requisite rate of research to mitigate the somber projections of resistant pathogens requires computational tools that can improve experimental design and efficiency in developing PDI treatments.

`PDIpy <https://pypi.org/project/pdipy/>`_ is offered to fulfill this role as the first API that simulates PDI biochemistry. PDIpy accepts user inputs of an experimental system and calculates parameters that are implemented into a `Tellurium <https://tellurium.readthedocs.io/en/latest/walkthrough.html>`_ kinetic system and visually expressed. Post-processing of the simulation data is supported with a built-in function, and the simulation contents can be exported for external post-processing at the discretion of the user. The directory of the `PDIpy GitHub <https://github.com/freiburgermsu/pdipy>`_ articulates a few case studies of replicating experimental data in Jupyter Notebooks. Users and developers are encouraged to critique and improve PDIpy, as an open-source library, through `GitHub issues <https://github.com/freiburgermsu/pdipy/issues>`_. 

____________


PDIpy API
--------------

++++++++++++++++++++++
Installation
++++++++++++++++++++++

pdipy is installed in a command prompt, Powershell, Terminal, or Anaconda Command Prompt via ``pip``::

 pip install pdipy

++++++++++++++++++++++
__init__
++++++++++++++++++++++

The simulation environment is defined:

.. code-block:: python

 import pdipy
 pdi = pdipy.PDI(total_time, solution_dimensions = {}, surface_system = False, well_count = 24, timestep = 3, verbose = False, jupyter = False)

- *total_time* ``float``: specifies the total simulated time.
- *simulation* ``dict``: defines the physical dimensions of the simulated solution, which are used to calculate photonic density and photosensitizer volume proportion.
- *surface_system* ``bool``: specifies whether a photodynamic system with a surface-bound, cross-linked, photosensitizer will be simulated.
- *well_count* ``int``: specifies the petri dish well count that will be simulated, which begets default dimensions of the simulated solution.
- *timestep* ``int``: specifies the timestep value of the simulation, which subtly affects the log-reduction predictions at the end of the simulated time.  
- *verbose* ``bool``: specifies whether simulation details and calculated values will be printed. This is valuable for troubleshooting.
- *jupyter* ``bool``: specifies whether the simulation is being conducted in a Jupyter Notebook, which allows ``display()`` to illustrate data tables and figures.

++++++++++++++++++++++
define_bacterium()
++++++++++++++++++++++

The characteristics of the simulated bacterium are defined from either a parameter file in the ``pdipy/parameters/bacteria`` directory or in a dictionary argument:

.. code-block:: python

 pdi.define_bacterium(bacterial_specie = None, bacterial_characteristics = {}, bacterial_cfu_ml = 1e6, biofilm = False)

- *bacterial_specie* ``str``: specifies one of the bacteria in the ``parameters`` directory to simulate, where *S. aureus* is imported by default. The imported parameters from this selection are overwritten by the ``bacterial_characteristics``, and thus the two arguments are not exclusive.
- *bacterial_characteristics* ``dict``: provides custom characteristics of the simulated bacterium, which supplants characteristics from the ``bacterial_specie`` argument. The expected dictionary keys ``shape``, ``membrane_thickness_nm``, ``cell_mass_pg``, ``"cell_volume_fL``, ``eps_oxidation_rate_constant``, and ``cellular_dry_mass_proportion_biofilm`` are all themselves dictionaries that follow a simple structure:

.. code-block:: json

 {
    	  "value": 0.268,
    	  "reference": "A. G. O’DONNELL, M. R. NAHAIE, M. GOODFELLOW, D. E. MINNIKINI, and V. HAJEK. Numerical Analysis of Fatty Acid Profiles in the Identification of Staphylococci. Journal of General Microbiology (1989). 131, 2023-2033. https://doi.org/10.1099/00221287-131-8-2023",
    	  "notes": "All saturated SCFAs were summed from Table 2 for all S. aureus entries."
 }

The ``reference`` and ``notes`` keys are optional, yet may be important for provenance and reproducibility of the simulation results. The final key of ``bacterial_characteristics`` is ``membrane_chemicals``, which contains as values each chemical group that constitutes the respective membrane, such as ``BCFA`` for branch-chain fatty acids and ``SCFA`` for straight-chain fatty acids. The sub-structure of these values are provided exemplified by the following content for the ``"SCFA"`` entry for *S. aureus*:

.. code-block:: json

 {
     "density_gL": {
          "value": 0.94,
          "reference": ["https://pubchem.ncbi.nlm.nih.gov/compound/Stearic-acid#section=Density"],
          "notes": "The density for all saturated fatty acids is estimated as stearic acid."
          },
     "formula": [
          "C20_H38_O2",
          "C18_H34_O2",
          "C16_H30_O2"
          ],
	 "proportion": {
	      "value": 0.268,
	      "reference": "A. G. O’DONNELL, M. R. NAHAIE, M. GOODFELLOW, D. E. MINNIKINI, and V. HAJEK. Numerical Analysis of Fatty Acid Profiles in the Identification of Staphylococci. Journal of General Microbiology (1989). 131, 2023-2033. https://doi.org/10.1099/00221287-131-8-2023",
	      "notes": "All saturated SCFAs were summed from Table 2 for all S. aureus entries."
          }
  }


- *bacterial_cfu_ml* ``float``: specifies the bacterial concentration for simulations of solution-based photosensitizers. 
- *biofilm* ``bool``: specifies whether a biofilm will be simulated.

+++++++++++++++++++++++++++++++
define_photosensitizer()
+++++++++++++++++++++++++++++++

Defines the simulated photosensitizer:

.. code-block:: python

 pdi.define_photosensitizer(photosensitizer = 'A3B_4Zn', photosensitizer_characteristics = {}, photosensitizer_molar = None, photosensitizer_g = 90e-9, cross_linked_sqr_m = 0.0191134)

- *photosensitizer* ``str``: specifies which photosensitizer from the predefined options in the ``pdipy/parameters/photosensitizers.json`` parameter file will be simulated.
- *photosensitizer_characteristics* ``dict``: defines characteristics of the simulation photosensitizer, which can be used to refine the parameterized photosensitizer. The expected structure of the dictionary are keys with dictionary substructure according to the following example:

.. code-block:: json

 {
		"e_quantum_yield": {
			"value": 0.6,
			"reference": "Singlet Oxygen Yields and Radical Contributions in the Dye-Sensitised Photo-oxidation in methanol of esters of polyunsaturated fatty acids _oleic, linoleic, linolenic, and arachidonic) Chacon et al., 1988"
		},
		"so_specificity": {
			"value": 0.8,
			"reference": null
		},
		"formula": {
			"value": "C76_H48_N16_F12_Zn",
			"reference": null
		},
		"soret_nm": {
			"value": [ 400, 430 ],
			"reference": null
		},
		"q_nm": {
			"value": [ 530, 625 ],
			"reference": null
		},
		"charge": 4,
		"photobleaching_constant (cm^2/J)": {
			"value": 1.74e-7,
			"reference": "“Photobleaching kinetics, photoproduct formation, and dose estimation during ALA induced PpIX PDT of MLL cells under well oxygenated and hypoxic conditions” by Dysart et al., 2005",
			"notes": "The 0.015 value from literature is divided by 8.64e4 -- the quantity of seconds in a day -- to yield a sensible value. A similar value is discovered from “PHOTOBLEACHING OF PORPHYRINS USED IN PHOTODYNAMIC THERAPY AND  IMPLICATIONS FOR THERAPY” by Mang et al., 1987"
			},
		"dimensions": {
			"shape": "disc",
			"length_A": 32.8,
			"width_A": 32.8,
			"depth_A": 1.5,
			"notes": "The depth is atomic thickness, as quantified by this paper https://www.nature.com/articles/ncomms1291."
		} 
 }

The ``value`` sub-key in the dictionary substructures, where it is present, is the only necessary sub-key for each parameter.

- *photosensitizer_molar* ``float``: specifies the photosensitizer molar concentration for solution simulations.
- *photosensitizer_g* ``float``: specifies the mass of photosensitizer that is surface-bound in cross-linked simulations.
- *cross_linked_sqr_m* ``float``: defines the square-meters area that is coated with the bound photosensitizer from the ``photosensitizer_g`` parameter, for cross-linked  simulations.
- *parameterized_ph_charge* ``bool``: specifies whether the pH will be charged balance, where ``True`` prevents the parameterization of alkalinity in the feed solution. 



++++++++++++++++++++++
define_light()
++++++++++++++++++++++

This function is used to parse and execute pre-existing input file:

.. code-block:: python

 pdipy.define_light(measurement, light_source = None, light_characteristics = {})

- *measurement* ``dict``: provides the unit and quantity of the photonic intensity measurement of the light source in a key-value pair. The supported unit options are: ``irradiance`` in :math:`\frac{mW}{cm^2}`, ``exposure`` in :math:`\frac{J}{cm^2}`, ``lux`` in :math:`\frac{lumen}{m^2}`, and ``lumens`` in :math:`lumens`.
- *light_source* ``str``: specifies a light source from the predefined options in the ``pdipy/parameters/light_source.json`` parameter file will be simulated. 
- *light_characteristics* ``dict``: specifies custom characteristics of the light source, which overwrite characteristics that are specified from the ``light_source`` option. The expected structure of the dictionary are keys with dictionary substructure according to the following example:

.. code-block:: json

 {
    "visible_proportion": {
      "value": 0.1,
      "reference": "Macisaac et al., 1999"
    },
    "lumens_per_watt": {
      "value": 3,
      "reference": "Michael F. Hordeski. Dictionary Of Energy Efficiency Technologies. Fairmont Press. ISBN: 9780824748104"
    }
  }

where the ``value`` sub-key in the dictionary substructures is the only necessary sub-key for each parameter.


++++++++++++++++++++++
simulate()
++++++++++++++++++++++

The aforementioned system specifications are refined into chemical parameters and are executed in a ``Tellurium`` kinetic model:

.. code-block:: python

 pdi.simulate(figure_title = None, y_label = 'log10', exposure_axis = False, display_fa_oxidation = False, display_ps_excitation = False)

- *figure_title* & *y_label* ``str``: specify the title and y-axis label of the simulation figure, respectively. The value of ``None`` defaults to **Cytoplasmic oxidation and inactivation of < bacterial genera_specie > via PDI**.
- *exposure_axis* ``bool``: specifies whether the x-axis of the simulation figure will be defined with cumulative exposure :math:`\frac{J}{cm^2}` over the simulation or in minutes of simulation time, where the latter is default.
- *display_fa_oxidation* & *display_ps_excitation* ``bool``: determine whether the fatty acid oxidation or the photosensitizer excitation proportions, respectively, will be plotted with the reduction data.


++++++++++++++++++++++
export()
++++++++++++++++++++++

The simulation contents, including the regression plot and information, are exported to the desired location:

.. code-block:: python

 pdi.export(self, export_name = None, export_directory = None)

- *export_name* & *export_directory* ``str``: specify the name and directory, respectively, to which the simulation contents will be saved, where ``None`` defaults to a folder name with simulation parameters **PDIpy-<photosensitizer_selection>-<bacterial_specie>-<count>** within the current workign directory.


++++++++++++++++++++++
data_parsing()
++++++++++++++++++++++

The processed data can be automatically processed through this function, as a convenient form of post-processing within the ``PDI`` object environment:

.. code-block:: python

 pdi.data_parsing(log_reduction = None, target_time = None)

- *log_reduction* ``float``: inquires at what time the specified log-reduction is achieved 
- *target_time* ``float``: inquires what log-reduction is achieved that the specified time

____________


Accessible content
----------------------

Numerous entities are stored within the ``PDI`` object, and can be subsequently used in a workflow. The complete list of content within the ``PDI`` object can be identified and printed through the built-in ``dir()`` function in the following example sequence:

.. code-block:: python

 # conduct a pdipy simulation
 from pdipy import PDI
 pdi = PDI(total_time, solution_dimensions = {}, surface_system = False, well_count = 24, timestep = 3, verbose = False, jupyter = False)
 pdi.define_bacterium(bacterial_specie = None, bacterial_characteristics = {}, bacterial_cfu_ml = 1e6, biofilm = False)
 pdi.define_photosensitizer(photosensitizer = 'A3B_4Zn', photosensitizer_characteristics = {}, photosensitizer_molar = None, photosensitizer_g = 90e-9, cross_linked_sqr_m = 0.0191134)
 pdipy.define_light(measurement, light_source = None, light_characteristics = {})
 pdi.simulate(figure_title = None, y_label = 'log10', exposure_axis = False, display_fa_oxidation = False, display_ps_excitation = False)
 pdi.export(self, export_name = None, export_directory = None)

 # evaluate the PDI object contents
 print(dir(pdi))

The following list highlights stored content in the ``PDI`` object after a simulation:

- *raw_data* & *processed_data* ``Pandas.DataFrame``: `Pandas DataFrames <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that contain the raw and processed simulation data, respectively. This files are also exported through the export function.
- *model* & *phrasedml_str* ``str``: The kinetic model and its corresponding `SED-ML <https://sed-ml.org/>`_ plot, respectively, composed in a string that can be read by Tellurium and converted into the standard XML formats of these languages.
- *bacterium*, *photosensitizer*, & *light* ``dict``: Dictionaries of the simulation parameters for the bacterium, photosensitizer, and light, respectively.
- *parameters*, *variables*, & *results* ``dict``: Dictionaries that possess the input parameters, calculation variables, and simulation results, respectively.
- *figure* & *ax* ``MatplotLib.Pyplot.subplots``: The `MatPlotLib objects <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplot.html#matplotlib.pyplot.subplot>`_ of the simulation figure, which allows the user to externally manipulate the figure without recreating a new figure from the raw or processed data.
- *chem_mw* ``ChemMW``: The ``ChemMW`` object from the `ChemW module <https://pypi.org/project/ChemW/>`_, which allows users to calculate the molecular weight from a string of any chemical formula. The formatting specifications are detailed in the README of the ChemW module. 
- *hf* ``HillFit``: The `HillFit object <https://pypi.org/project/hillfit/>`_ is stored, from which the Hill-equation regrssion parameters, equation string, and R\ :sup:`2`\ of the fitted equation can be programmatically accessed, in addition to being exported with the ``PDIpy`` content through the ``export()`` function.
- *bacteria* ``list``: A list of all the predefined bacteria parameter files, from which a user can easily simulate via the ``PDI`` object.
- *light_parameters*, *photosensitizers*, & *solution* ``dict``: Dictionaries of the predefined options and parameters for the light sources, photosensitizers, and solution dimensions, respectively.