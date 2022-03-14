PDIpy API
--------------

++++++++++++++++++++++
PDI()
++++++++++++++++++++++

The simulation environment is defined:

.. code-block:: python

 import pdipy
 pdi = pdipy.PDI(total_time, solution_dimensions = {}, surface_system = False, 
                 well_count = 24, timestep = 3, verbose = False, jupyter = False)

- *total_time* ``float``: specifies the total simulated time in minutes.
- *simulation* ``dict``: defines the physical dimensions of the simulated solution, which are used to calculate photonic density and photosensitizer volume proportion.
- *surface_system* ``bool``: specifies whether a photodynamic system with a surface-bound, cross-linked, photosensitizer will be simulated.
- *well_count* ``int``: specifies the petri dish well count that will be simulated, which begets default dimensions of the simulated solution.
- *timestep* ``int``: specifies the timestep value in minutes of the simulation, which subtly affects the log-reduction predictions at the end of the simulated time.  
- *verbose* ``bool``: specifies whether simulation details and calculated values will be printed. This is valuable for troubleshooting.
- *jupyter* ``bool``: specifies whether the simulation is being conducted in a Jupyter Notebook, which allows ``display()`` to illustrate data tables and figures.


define_conditions()
++++++++++++++++++++++

The characteristics of the simulated system are concisely defined in a single function:

.. code-block:: python

 pdi.define_conditions(bacterial_specie = None, bacterial_characteristics = {}, 
       bacterial_cfu_ml = 1e6, biofilm = False, photosensitizer = 'A3B_4Zn', 
       photosensitizer_characteristics = {}, photosensitizer_molar = None, photosensitizer_g = 90e-9, 
       cross_linked_sqr_m = 0.0191134, light_source = None, light_characteristics = {}, measurement = {})

- *bacterial_specie* ``str``: specifies one of the bacteria in the ``pdipy/parameters/bacteria`` directory, where *S. aureus* is imported by default. These imported parameters are overwritten by the ``bacterial_characteristics``, and can be used in conjunction.
- *bacterial_characteristics* ``dict``: provides custom characteristics of the simulated bacterium that supplant imported characteristics from the ``bacterial_specie`` argument. The expected dictionary keys ``shape``, ``membrane_thickness_nm``, ``cell_mass_pg``, ``"cell_volume_fL``, ``eps_oxidation_rate_constant``, and ``cellular_dry_mass_proportion_biofilm`` are all themselves dictionaries that follow a simple structure:

.. code-block:: json

 {
    	  "value": 0.268,
    	  "reference": "A. G. O’DONNELL, M. R. NAHAIE, M. GOODFELLOW, D. E. MINNIKINI, and V. HAJEK. Numerical Analysis of Fatty Acid Profiles in the Identification of Staphylococci. Journal of General Microbiology (1989). 131, 2023-2033. https://doi.org/10.1099/00221287-131-8-2023",
    	  "notes": "All saturated SCFAs were summed from Table 2 for all S. aureus entries."
 }

The ``reference`` and ``notes`` keys are optional, yet may be important for provenance and reproducibility of the simulation results. The final key of ``bacterial_characteristics`` is ``membrane_chemicals``, which characterizes the chemicals that constitutes the bacterial membrane such as ``BCFA`` for branch-chain fatty acids and ``SCFA`` for straight-chain fatty acids. The sub-structure of these values are exemplified by the following content for the ``"SCFA"`` entry for *S. aureus*:

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
- *photosensitizer* ``str``: specifies which photosensitizer from the predefined options in the ``pdipy/parameters/photosensitizers.json`` parameter file will be simulated.
- *photosensitizer_characteristics* ``dict``: defines characteristics of the simulation photosensitizer, which supplant the characteristics from the ``photosensitizer`` parameter. The expected structure of the dictionary are keys with dictionary substructure according to the following example:

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

- *photosensitizer_molar* ``float``: specifies the photosensitizer molar concentration for simulations of a solution-based photosensitizer.
- *photosensitizer_g* ``float``: specifies the mass of photosensitizer that is surface-bound in cross-linked simulations.
- *cross_linked_sqr_m* ``float``: defines the square-meters area that is coated with the bound photosensitizer from the ``photosensitizer_g`` parameter, for cross-linked  simulations.
- *parameterized_ph_charge* ``bool``: specifies whether the pH will be charged balance, where ``True`` prevents the parameterization of alkalinity in the feed solution. 
- *measurement* ``dict``: provides the unit and quantity of the photonic intensity measurement of the light source in a key-value pair. The supported unit options are: ``irradiance`` in mW/cm\ :sup:`2`\, ``exposure`` in J/cm\ :sup:`2`\, ``lux`` in lumen/m\ :sup:`2`\, and ``lumens`` in lumens.
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


simulate()
++++++++++++++++++++++

The aforementioned system specifications are refined into chemical parameters and are executed in a ``Tellurium`` kinetic model:

.. code-block:: python

 pdi.simulate(export_name = None, export_directory = None, figure_title = None, 
            y_label = 'log10', exposure_axis = False, display_fa_oxidation = False, 
            display_ps_excitation = False, export_content = True)

- *export_name* & *export_directory* ``str``: specify the name and directory, respectively, to which the simulation contents will be saved, where ``None`` defaults to a folder name with simulation parameters **PDIpy-<photosensitizer_selection>-<bacterial_specie>-<count>** within the current workign directory.
- *figure_title* & *y_label* ``str``: specify the title and y-axis label of the simulation figure, respectively. The y-axis label is vague to support generalization to plots where the fatty acid oxidation and photosensitizer excitation content is overlaid, and thus would be not appropriately described by more descriptive labels. The value of ``None`` defaults to **Cytoplasmic oxidation and inactivation of < bacterial genera_specie > via PDI**. 
- *exposure_axis* ``bool``: specifies whether the x-axis of the simulation figure will be defined with cumulative exposure J/cm\ :sup:`2`\ over the simulation or in minutes of simulation time, where the latter is default.
- *display_fa_oxidation* & *display_ps_excitation* ``bool``: determine whether the fatty acid oxidation or the photosensitizer excitation proportions, respectively, will be plotted with the reduction data.
- *export_content* ``bool``: specifies whether the simulation content will be exported.

parse_data()
++++++++++++++++++++++

The processed data can be automatically processed through this function, as a convenient form of post-processing within the ``PDI`` object environment:

.. code-block:: python

 value, unit = pdi.data_parsing(log_reduction = None, target_hours = None)

- *log_reduction* ``float``: inquires at what time the specified log-reduction is achieved 
- *target_hours* ``float``: inquires what log-reduction is achieved that the specified time

**Returns** *value* ``float``: The value of the search inquiry, reported in the respective units.

**Returns** *unit* ``str``: The units of the search inquiry result, being either log-reduction or hours.