Usage 
---------

.. code-block:: python

 from pdipy import PDI

 # define the simulation conditions
 pdi = PDI(        
     solution_dimensions = {
         'area_sqr_cm':pi/4,
         'depth_cm': 46,
        },
        printing = False
        )
 pdi.define_conditions(
     bacterial_specie = 'S_aureus', 
     bacterial_cfu_ml = 3e7, 
     photosensitizer = 'protoporphyrin IX', 
     photosensitizer_characteristics = {
        'formula': {
            'value': 'C34_H34_N4_O4'   
        },
        'dimensions':{
            'length (A)': length,
            'width (A)': width,
            'depth (A)': 4,
            'shape': 'rect',
        }
     },
     photosensitizer_molar = 18e-9, 
     measurement = {'irradiance': 8}, 
     )
 
 # execute and export the simulation
 pdi.simulate(
    export_name = 'Bozja_et_al, 1E4 Lux',
    figure_title = '1E4 Lux, Protoporphyrin IX',
    display_ps_excitation = True,
 )

 # parse the data and evaluate the PDI object contents
 value, unit = pdi.parse_data(log_reduction = 5)
 print(dir(pdi))
 
Accessible content
++++++++++++++++++++++

Numerous entities are stored within the ``PDI`` object, and can be subsequently used in a workflow. The complete list of content within the ``PDI`` object can be identified and printed through the Python function ``dir()``, while the following list highlights a few notable contents within the ``PDI`` object:

- *raw_data* & *processed_data* ``Pandas.DataFrame``: `Pandas DataFrames <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that contain the raw and processed simulation data, respectively. These files are also exported through the export function.
- *model* & *phrasedml_str* ``str``: The kinetic model and its corresponding `SED-ML <https://sed-ml.org/>`_ plot, respectively, composed in a string that is read by Tellurium.
- *bacterium*, *photosensitizer*, & *light* ``dict``: Simulation parameters for the bacterium, photosensitizer, and light, respectively, in a for that is read by the code.
- *parameters*, *variables*, & *results* ``dict``: Input parameters, calculation variables, and simulation results, respectively, that are subsequently processed into DataFrames and exported as CSVs.
- *figure* & *ax* ``MatplotLib.Pyplot.subplots``: The `MatPlotLib objects <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplot.html#matplotlib.pyplot.subplot>`_ of the simulation figure, which allows the user to externally manipulate the figure without recreating a new figure from the raw or processed data.
- *chem_mw* ``chemw.ChemMW``: The ``ChemMW`` object from the `ChemW module <https://pypi.org/project/ChemW/>`_, which allows users to calculate the molecular weight from a string of any chemical formula or common name. 
- *bacteria* ``list``: A list of all the predefined bacteria parameter files, from which a user can easily simulate via the ``PDI`` object.
- *light_parameters*, *photosensitizers*, & *solution* ``dict``: Dictionaries of the predefined options and parameters for the light sources, photosensitizers, and solution dimensions, respectively.