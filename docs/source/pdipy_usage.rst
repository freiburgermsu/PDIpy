Usage 
---------

.. code-block:: python

 from pdipy import PDI

 # define the simulation conditions
 pdi = PDI(total_time = 360)
 pdi.define_conditions(bacterial_specie = 'S_aureus', bacterial_cfu_ml = 1e7, 
 photosensitizer = 'A3B_4Zn', photosensitizer_molar = 18e-9, measurement = {'irradiance': 8}, light_source = 'LED')
 
 # execute and export the simulation
 pdi.simulate()

 # parse the data and evaluate the PDI object contents
 value, unit = pdi.parse_data(log_reduction = 5)
 print(dir(pdi))
 
Accessible content
++++++++++++++++++++++

Numerous entities are stored within the ``PDI`` object, and can be subsequently used in a workflow. The complete list of content within the ``PDI`` object can be identified and printed through the built-in ``dir()`` function, while the following list highlights stored content in the ``PDI`` object after a simulation:

- *raw_data* & *processed_data* ``Pandas.DataFrame``: `Pandas DataFrames <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that contain the raw and processed simulation data, respectively. This files are also exported through the export function.
- *model* & *phrasedml_str* ``str``: The kinetic model and its corresponding `SED-ML <https://sed-ml.org/>`_ plot, respectively, composed in a string that can be read by Tellurium and converted into the standard XML formats of these languages.
- *bacterium*, *photosensitizer*, & *light* ``dict``: Dictionaries of the simulation parameters for the bacterium, photosensitizer, and light, respectively.
- *parameters*, *variables*, & *results* ``dict``: Dictionaries that possess the input parameters, calculation variables, and simulation results, respectively.
- *figure* & *ax* ``MatplotLib.Pyplot.subplots``: The `MatPlotLib objects <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplot.html#matplotlib.pyplot.subplot>`_ of the simulation figure, which allows the user to externally manipulate the figure without recreating a new figure from the raw or processed data.
- *chem_mw* ``chemw.ChemMW``: The ``ChemMW`` object from the `ChemW module <https://pypi.org/project/ChemW/>`_, which allows users to calculate the molecular weight from a string of any chemical formula. The formatting specifications are detailed in the README of the ChemW module. 
- *hf* ``hillfit.HillFit``: The ``HillFit`` object from the `Hillfit module <https://pypi.org/project/hillfit/>`_ is stored, from which the Hill-equation regrssion parameters, equation string, and R\ :sup:`2`\ of the fitted equation can be programmatically accessed, in addition to being exported with the ``PDIpy`` content through the ``export()`` function.
- *bacteria* ``list``: A list of all the predefined bacteria parameter files, from which a user can easily simulate via the ``PDI`` object.
- *light_parameters*, *photosensitizers*, & *solution* ``dict``: Dictionaries of the predefined options and parameters for the light sources, photosensitizers, and solution dimensions, respectively.