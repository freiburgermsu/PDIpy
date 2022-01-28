ROSSpy
_______

-----------
Motivation
-----------

Desalination is an unavoidable technology for meeting the 6th UN Sustainable Development Goal of providing potable for all people. Reverse Osmosis (RO) is the leading desalination technology, although, it can be further improved in energy efficiency and economic practicality by mitigating membrane fouling like mineral scaling. The geochemistry of mineral scaling is generally inaccessible to physical experimentation, and existing software programs to simulate scaling geochemistry -- e.g. French Creek -- are esoteric and/or financially expensive. 

We therefore developed a `Reverse Osmosis Scaling Software in Python (ROSSpy) <https://pypi.org/project/ROSSpy/>`_ as an open-source software that simulates scaling geochemistry through `PHREEQC <https://www.usgs.gov/software/phreeqc-version-3>`_. ROSSpy essentially translates user specifications of an RO system into a 1D reactive transport model of the RO membrane-solution interface in the feed channel that is then translated into PHREEQ code and executed via `PHREEQpy <https://pypi.org/project/phreeqpy/>`_. Examples abound in the `ROSSpy GitHub <https://github.com/freiburgermsu/ROSSpy>`_, and demonstrate the breadth and accuracy of ROSSpy for research. The underlying calculations and logic are detailed in our manuscript of ROSSpy [under-review]. We encourage users and developers to share critiques and suggestions for improving ROSSpy as an open-source community resource that may expedite RO research and resolving water insecurities.

++++++++++++++++
Installation
++++++++++++++++

ROSSpy is installed in a command prompt, Powershell, Terminal, or Anaconda Command Prompt (to use the program in Anaconda environments like Jupyter Notebook and Spyder) via the PyPI command::

 pip install rosspy

-----------
Concept
-----------

The ROSSpy framework represents RO desalination as a 1D reactive transport model of the membrane-solution interface in the feed channel of an RO module. The feed solution can be represented either by a single, homogeneous, solution domain or by a dual solution domain -- where the bulk solution is simulated separately from the concentration polarization (CP) that forms adjacent to the filtration membrane in the feed channel. Experimental evidence supports that the dual-domain is more fundamentally accurate, however, representing this model in PHREEQC has proven elusive, thus ROSSpy currently only supports the single-domain. The accuracy of the single-domain is nevertheless verified in the *examples/scaling/scaling_validation* directory of our `ROSSpy GitHub <https://github.com/freiburgermsu/ROSSpy>`_. The inlet boundary is defined by the Dirichlet condition, where the feed is assumed to be an infinite reservoir, and the outlet boundary is defined by the Cachy condition, where the effluent is assumed to be dependent upon the reactive transport processes of the RO module. 


----------------------
Functions
----------------------

ROSSpy is organized into Python functions, within the ``ROSSPkg`` class object, that serve specified purposes or types of calculations. Each of these functions are detailed in the following sub-sections.


+++++++++++
__init__
+++++++++++

The simulation environment is defined::

 import rosspy
 ross = rosspy.ROSSPkg(operating_system = 'windows', verbose = False, jupyter = False)

- *operating_system* ``str``: specifies whether the user is using a ``windows`` or ``unix`` system, which directs subtle differences in importing the PHREEQpy package and incorporating comments to the ``PQI`` PHREEQ input files.
- *verbose* ``bool``: specifies whether simulation details and calculated values will be printed. This is valuable for trobuleshooting.
- *jupyter* ``bool``: specifies whether the simulation is being conducted in a Jupyter Notebook, which allows ``display()`` of data tables and figures.


++++++++++++++++
define_general
++++++++++++++++

The general conditions of the simulation are defined::

 ross.define_general(database_selection, simulation = 'scaling', domain_phase = None, 
 quantity_of_modules = 1, simulation_type = 'transport', simulation_title = None)

- *database_selection* ``str``: specifies which PHREEQ database file -- ``Amm``, ``ColdChem``, ``core10``, ``frezchem``, ``iso``, ``llnl``, ``minteq``, ``minteq.v4``, ``phreeqc``, ``pitzer``, ``sit``, ``Tipping_Hurley``, or ``wateq4f`` -- will be imported and used to execute the simulation.
- *simulation* ``str``: specifies whether the ``scaling`` or ``brine`` of the simulation will be evaluated.
- *domain_phase* ``str``: specifies whether the ``mobile`` (i.e. bulk solution) or the ``immobile`` (i.e. the CP solution layer) will be evaluated for dual domain simulations. Parameterizing an argument other than ``None`` implicitly signifies that simulation of the dual-domain model, as opposed to the default single-domain model.  
- *quantity_of_modules* ``int``: specifies the number of in-series RO modules that will be simulated.
- *simulation_type* ``str``: specifies whether RO reactive transport ``transport``, or the geochemistry of ``evaporation``, will be simulated with the parameterized feed solution.
- *simulation_title* ``str``: specifies the title of the simulation, which is only observed in the PHREEQC ``PQI`` input file.


+++++++++++
transport
+++++++++++

Spatiotemporal conditions for reactive transport simulations are defined::

 ross.transport(simulation_time, simulation_perspective = None, 
 module_characteristics = {}, timestep = None, cells_per_module = 12, 
 kinematic_flow_velocity = None, exchange_factor = 1e5)

- *simulation_time* ``float``: specifies the total simulated time in seconds.
- *simulation_perspective* ``str``: specifies whether the simulation data is parsed to view the end of the module over the entire simulated time "all_time" or to view the entire module distance at the final timestep ``all_distance``. The ``None`` parameter defaults to ``all_time`` for brine simulations and ``all_distance`` for scaling simulations.
- *module_characteristics* ``dict``: specifies custom RO specifications that diverge from those of the default DOW FILMTEC BW30-400 RO module. The expected ``keys`` of the dictionary are 

 + 'module_diameter_mm'
 + 'permeate_tube_diameter_mm'
 + 'module_length_m'
 + 'permeate_flow_m3_per_day' 
 + 'max_feed_flow_m3_per_hour'
 + 'membrane_thickness_mm' 
 + 'feed_thickness_mm'
 + 'active_m2'
 + 'permeate_thickness_mm'
 + 'polysulfonic_layer_thickness_mm'
 + 'support_layer_thickness_mm'. 

 The ``values`` for the dictionary are all floats in the units that are listed at the end of the corresponding key.
 
- *timestep* ``float``: specifies the simulation timestep in seconds. The ``None`` parameter defaults to the maximum timestep that still adheres to the Courant Condition of maintaining simulated resolution.
- *cells_per_module* ``int``: specifies the quantity of cells into an RO module is discretized. This primarily controls the resolution of data and plots over the distance of the module, and thus is only consequential for ``simulation_perspective = "all_distance"``.
- *kinematic_flow_velocity* ``float``: specifies the kinetic flow velocity of the feed solution. The ``None`` parameter defaults to 9.33E-7 (m^2/sec).
- *exchange_factor* ``float``: specifies the kinetic rate of exchange (1/sec) between the mobile and immobile phases of a dual domain simulation.


+++++++++++
reaction
+++++++++++

The permeate flux gradient in reactive transport simulations, or the rate of evaporation in evaporation simulations, is calculated::

 ross.reaction(final_cf = None, permeate_efficiency = 1, head_loss = 0.89, evaporation_steps = 15)

- *final_cf* ``float``: specifies the effluent CF of the last module in the simulated RO system. The ``None`` parameter indicates that the ``linear_permeate`` permeate flux method will be used, while any numerical value indicates that a ``linear_cf`` permeate flux method will be used. 
- *permeate_efficiency* ``float``: specifies 0<=PE<=1 proportion of calculated permeate flux that actually filters from the feed solution. This is useful for distinguishing fresh RO modules from aged and partly compromised RO modules in ROSSpy simulation.
- *head_loss* ``float``: specifies the 0<=HL<=1 proportion of effluent pressure relative to the influent. The default value of 0.89 (“Reverse osmosis desalination: Modeling and experiment” by Fraidenraich et al., 2009) corresponds to an 11% pressure drop.


+++++++++++
solutions
+++++++++++

The geochemistry of the feed solution is parameterized, either through specifying a complete parameter file that is imported from the *rosspy/water_bodies* directory of the ROSSpy package, or through passing a dictionary of same content as organization as an argument::

 ross.solutions(water_selection = '', water_characteristics = {}, 
 solution_description = '', parameterized_ph_charge = True)

- *water_selection* ``str``: specifies which feed water from the *rosspy/water_bodies* directory will be simulated. ROSSpy offers by default  parameter files for natural waters -- i.e. the ``red_sea`` and the ``mediterranean_sea`` -- and produced waters from fracking oil wells -- i.e. ``bakken_formation``, ``marcellus_appalachian_basin``, ``michigan_basin``, ``north_german_basin``, ``palo_duro_basin``, or ``western_pennsylvania_basin``. Other parameter files can be created and called in simulations by emulating the syntax of these default files and storing the created parameter files in the aforementioned directory of these files.
- *water_characteristics* ``dict``: defines the geochemistry and conditions that will simulate the feed solution. The expected ``keys`` are 

 + 'element'
 + 'temperature'
 + 'pe'
 + 'Alkalinity' 
 + 'pH'
 
Each of the ``values`` of these keys is itself a dictionary, with the keys of "value" that denotes the numerical value of the entry and "reference" that denotes the experimental citation for the numerical value. The "element" key deviates slightly from this organization, where another layer of dictionaries is introduced for each ion in the feed. Each ion dictionary possesses the "concentration (ppm)" key to specify the ppm concentration of the designated ion and the "form" key to signify the mineral form of the parameterized ion, in addition to the aforementioned "reference" key. The following dictionary illustrates this organization, which is also exemplified in the default water body parameter files.

{
    "element": {
        "Mn": {
            "concentration (ppm)": 3000,

            "reference": "Haluszczak, Rose, and Kump, 2013 [estimated from another Marcellus publication]"
			
        }, 

        "Si": {
            "concentration (ppm)": 95,

            "reference": "Haluszczak, Rose, and Kump, 2013 [reported average from another Marcellus publication]",

            "form": "SiO2"
			
        }
		
    },

    "temperature": {
        "value": 24,

        "reference": "Dresel and Rose, 2010"
		
    }
	
}

- *solution_description* ``str``: briefly describes the solution in the name of the simulation folder.
- *parameterized_ph_charge* ``bool``: specifies whether the pH will be charged balance, which consequently prevents the parameterization of alkalinity in the feed solution. 


+++++++++++++++++++++
equilibrium_phases
+++++++++++++++++++++

The minerals, and pre-existing equilibria conditions, that will be explored in scaling simulations are defined::

 ross.equilibrium_phases(block_comment = '', ignored_minerals = [], 
 existing_parameters = {})

- *block_comment* ``str``: describes the minerals or scaling phenomena only in the ``PQI`` file.
- *ignored_minerals* ``list``: defines the minerals that will be excluded from the determined set of minerals that can potentially precipitate from the parameterized feed ions.
- *existing_parameters* ``dict``: specifies pre-existing equilibria conditions that influence the geochemical calculations of PHREEQ. The expected ``keys`` are the referenced mineral names and the respective ``values`` are

 + 'saturation'
 + 'initial_moles'
 
which correspond to the pre-existing saturation index and the initial moles of the respective mineral in the simulated system.


++++++++++++++++
selected_output
++++++++++++++++

The simulation content that will be incorporated to the output file is defined::

 ross.selected_output(output_filename = None)

- *output_filename* ``str``: specifies the name of a simulation output file that will be created when the input file is executed.


+++++++++++
export
+++++++++++

The simulation parameters, raw and processed data, figures, and the input file are exported into a designated labeled folder for the simulation::

 ross.export(simulation_name = None, input_path = None, 
 output_path = None, external_file = False)

- *simulation_name* ``str``: specifies the name of the simulation folder to which simulation content will be exported. The ``None`` parameter assigns a default name for the simulation folder, which follows the format of **date-ROSSpy-water_selection-simulation_type-database_selection-simulation-simulation_perspective-#**. 
- *input_path* ``str``: specifies the directory path to where the input file will be exported. The ``None`` parameter exports the input file as "input.pqi" to the curent working directory. 
- *output_path* ``str``: specifies the directory path to where the output file will be exported. The ``None`` parameter exports the output file as "selected_output.csv" to the curent working directory.
- *external_file* ``str``: specifies whether the input file of the simulation was imported and parsed from a pre-existing ``PQI`` file, and thus was not created through the aforementioned ROSSpy functions.


++++++++++++++++
parse_input
++++++++++++++++

A pre-existing input file is parsed and interpreted for simulation::

 ross.parse_input(input_file_path, simulation, water_selection = None, 
 simulation_name = None, active_m2= None)

- *input_file_path* ``str``: specifies the path of the input file that will be imported, parsed, and simulated. 
- *simulation* ``str``: defines parsing the simulation data for either ``scaling`` or ``brine``. 
- *water_selection* ``str``: describes the feed water of the input file that will be simulated. 
- *simulation_name* ``str``: specifies the name of the simulation folder to which all of the simulation files will be exported, identical to this parameter for the aforementioned ``export()`` function of ROSSpy. The ``None`` parameter likewise defaults to the aforementioned naming scheme. 
- *active_m2* ``float``: defines the area of active filtration in the simulated RO module. The ``None`` parameter defaults to 37 from the FILMTEC BW30-400 module. 


+++++++++++
execute
+++++++++++

The input file is executed through PHREEQ::

 ross.execute(simulated_to_real_time = 9.29)

- *simulated_to_real_time* ``float``: specifies the ratio of simulated time to real computational time when executing ROSSpy simulations. This is used to approximate the time that is required for the simulation to complete. The default ``9.29`` ratio represented simulations of multiple days or weeks, while shorter simulations of minutes/hours have a higher ratio, perhaps around ``20``.

- Note: The raw simulation data is returned by this function as a ``pandas.DataFrame`` object, which can be arbitrarily manipulated by the user through `pandas operations <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_.


++++++++++++++++++++++++++
process_selected_output
++++++++++++++++++++++++++

The simulation data is processed into figures and corresponding data tables::

 ross.process_selected_output(selected_output_path = None, scale_ions = True, 
 plot_title = None, title_font = 'xx-large', label_font = 'x-large', 
 x_label_number = 6, export_name = None, export_format = 'svg', individual_plots = None)

- *selected_output_path* ``str``: specifies the path of a simulation output file that will be processed into data tables and figures. This imported file can be independent of executing ROSSpy, and thus can be used to process old data. The ``None`` parameter constrains the function to only process data that resides in the ROSSpy object from a recent execution.
- *scale_ions* ``bool``: specifies whether the scale, from ``scaling`` simulations, will be collectivized and refined into quantities of individual ions that constitute the mineral scale. The ionic quantities are exported in a JSON file to the simulation folder, with the other simulation content. The default value is ``True``.
- *plot_title* ``str``: specifies the title of the simulation figure. The ``None`` parameter defaults to titles that are customized with the simulation (``scaling`` or ``brine``), the water body, and the total simulation time.
- *title_font* & *label_font* ``str``: these specify the fonts of the figure title and axis labels, respectively, in terms of MatPlotLib font identifications: ``xx-small``, ``x-small``, ``small``, ``medium``, ``large``, ``x-large``, or ``xx-large``. 
- *x_label_number* ``int``: quantifies the ticks along the x-axis in the simulation figure.
- *export_name* ``str``: specifies the export name of the simulation figure. The default names are ``brine`` for ``brine`` simulations, and ``all_minerals`` or an individual mineral name (e.g. ``Gypsum``) for ``scaling`` simulations, depending upon a ``False`` or ``True`` value of the *individual_plots* argument, respectively.
- *export_format* ``str``: specifies the format of the exported simulation figure, from the MatPlotLib options: ``svg``, ``pdf``, ``png``, ``jpeg``, ``jpg``, or ``eps``.
- *individual_plots* ``bool``: specifies whether each mineral of ``scaling`` simulations are plotted and exported in individual figures, or whether all precipitated minerals are plotted and exported together in a combined single figure. The ``None`` parameter defaults to "True" for the "all_time" *simulation_perspective* or "False" otherwise.
- Note: The processed simulation data, which predicated the figures, is returned by this function as a ``pandas.DataFrame`` object, which can be arbitrarily manipulated by the user through `pandas operations <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_.


----------------------
Execution
----------------------

ROSSpy is executed through a deliberate sequence of the aforementioned functions::
 
 import rosspy
 ross = rosspy.ROSSPkg()
 ross.define_general(database_selection, simulation)
 ross.transport(simulation_time, simulation_perspective, )
 ross.reaction(final_cf)
 ross.solutions(water_selection, water_characteristics)
 ross.equilibrium_phases()
 ross.selected_output()
 ross.export()
 raw_data = ross.execute()
 processed_data = ross.process_selected_output()

ROSSpy can be tested with a simple built-in ``test()`` function, which can be executed through these three lines::

 import rosspy
 ross = rosspy.ROSSPkg(operating_system = 'windows', verbose = False, jupyter = False)
 ross.test()

The ``Test()`` function executes a predefined sample simulation, which should exemplify the use of the ROSSpy API and the analysis of the exported outputs.