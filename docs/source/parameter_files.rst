PDIpy API
--------------

Parameters may be more succinctly provided to the code through JSON files that are automatically searched by the code. These files complement hard-coded parmeters through the package functions. Each parameter files pertains to a distinct aspect of the calculations, which are individually detailed in the following sections.


photosensitizers
++++++++++++++++++++++

The simulation environment is defined:

.. code-block:: json



- *total_time* ``float``: specifies the total simulated time in minutes.
- *simulation* ``dict``: defines the physical dimensions of the simulated solution, which are used to calculate photonic density and photosensitizer volume proportion.
- *surface_system* ``bool``: specifies whether a photodynamic system with a surface-bound, cross-linked, photosensitizer will be simulated.
- *well_count* ``int``: specifies the petri dish well count that will be simulated, which begets default dimensions of the simulated solution.
- *timestep* ``int``: specifies the timestep value in minutes of the simulation, which subtly affects the log-reduction predictions at the end of the simulated time.  
- *verbose* ``bool``: specifies whether simulation details and calculated values will be printed. This is valuable for troubleshooting.
- *jupyter* ``bool``: specifies whether the simulation is being conducted in a Jupyter Notebook, which allows ``display()`` to illustrate data tables and figures.


wells
+++++++++




Bacteria
++++++++++



light
+++++++