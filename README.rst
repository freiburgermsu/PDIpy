Simulating Photodynamic Inactivation (PDI) from chemical kinetics
------------------------------------------------------------------------

|PyPI version| |Actions Status| |docs| |Downloads| |License|

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
   
.. |docs| image:: https://readthedocs.org/projects/pdipy/badge/?version=latest
   :target: https://pdipy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Antibiotic resistance is a developing medical crisis that is projected to surpass cancer in annual deaths by mid-21st century. Photodynamic Inactivation (PDI) escapes resistance evolution and may therefore be an essential antibiotic method to hamper the growing threat of resistant pathogens. The requisite rate of research to mitigate these somber projections requires computational tools that can improve and expedite experimental research in developing PDI treatments.

`PDIpy <https://pypi.org/project/pdipy/>`_ is the first comprehensive software of PDI that simulates PDI biochemistry from a chemical kinetics model. PDIpy accepts user inputs of the simulated system, constructs and executes a `Tellurium <https://tellurium.readthedocs.io/en/latest/walkthrough.html>`_ kinetic system, and processes and exports the results in spreadsheets and SVG figures. Investigation of the simulation data is further supported with a PDIpy function that parses the data based upon either a specified inactivation or time. The `examples directory <https://github.com/freiburgermsu/pdipy/examples>`_ of the PDIpy GitHub exemplifies PDIpy through replicating experimental observations and conducting multiple sensitivity analyses. Users and developers are encouraged to critique, improve, and join the open-source project through `GitHub issues <https://github.com/freiburgermsu/pdipy/issues>`_ or emailing afreiburger@uvic.ca, respectively. 


.. note::

   This project is under active development, and may therefore be subject to incompatible changes in the API. We seek to minimize these types of upgrades.


++++++++++++++++++++++
Installation
++++++++++++++++++++++

pdipy is installed in a command prompt, Powershell, Terminal, or Anaconda Command Prompt via ``pip``::

 pip install pdipy


The full documentation is provided by `ReadTheDocs <https://pdipy.readthedocs.io/en/latest/index.html>`_.