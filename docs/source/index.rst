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

Antibiotic resistance is developing medical crisis that is projected to surpass cancer in annual deaths by mid-21st century. Photodynamic Inactivation (PDI) is a promising treatment method that escapes resistance evolution and may be an essential technology to hamper the growing threat of resistant pathogens. The requisite rate of research to mitigate these somber projections requires computational tools that can improve and expedite experimental research in developing PDI treatments.

`PDIpy <https://pypi.org/project/pdipy/>`_ is offered as the first comprehensive software of PDI to fulfill this by simulating PDI biochemistry from a chemical kinetics model. PDIpy accepts user inputs of an experimental system, executes a `Tellurium <https://tellurium.readthedocs.io/en/latest/walkthrough.html>`_ kinetic system, and expresses and exports the results through CSV spreadsheets and SVG images. Post-processing of the simulation data is further supported with a function that parses the based upon user parameters. The `examples directory <https://github.com/freiburgermsu/pdipy/examples>`_ of the PDIpy GitHub exemplifies PDIpy through replicating experimental observations. Users and developers are encouraged to critique and improve PDIpy, as an open-source library, through `GitHub issues <https://github.com/freiburgermsu/pdipy/issues>`_. 


.. note::

   This project is under active development.


++++++++++++++++++++++
Installation
++++++++++++++++++++++

pdipy is installed in a command prompt, Powershell, Terminal, or Anaconda Command Prompt via ``pip``::

 pip install pdipy
   
   

Contents
--------

.. toctree::

    pdipy_api
    parameter_files