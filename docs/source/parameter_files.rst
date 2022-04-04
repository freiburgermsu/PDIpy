PDIpy parameter files
-----------------------

Parameters may be more succinctly provided through comprehensive JSON files, which are automatically imported by the code, than through function arguments; although, these parameter sources are complementary. Each JSON file pertains to a distinct category of simulation parameters, which are individually detailed in the following sections.


photosensitizers
++++++++++++++++++++++

The chemical attributes of the simulated photosensitizer are articulated in the following format:

.. code-block:: json

 {
	"A3B_4Zn": {
		"e_quantum_yield": {
			"value": 0.6,
			"reference": "Singlet Oxygen Yields and Radical Contributions in the Dye-Sensitised Photo-oxidation in methanol of esters of polyunsaturated fatty acids _oleic, linoleic, linolenic, and arachidonic) Chacon et al., 1988"
		},
		"so_specificity": {
			"value": 0.8,
			"reference": null
		},
		"so_quantum_yield":{
			"value": 0.48,
			"reference": null
		},
		"formula": {
			"value": "C76_H48_N16_F12_Zn",
			"reference": null
		},
		"excitation_nm": {
			"value": [
				[400, 430],
				[530, 625]
				 ],
			"reference": null
		},
		"ps_rise (fs)": {
    		        "value": 50,
    		        "reference": "Anderssonet al., 1999, 'Photoinduced ... State' ; Gurzadyan et al., 1998, 'Time-resolved ... Zn-tetraphenylporphyrin'" 
		},
		"ps_decay (ns)": {
    		        "value": 1.5,
    		        "reference": "Akimoto et al., 1999, 'Ultrafast ... Porphyrins'"
		},
		"ps_charge_transfer (ns)": {
    		        "value": 100,
    		        "reference": "Kupper et al., 2002, 'Kinetics ... Oxygen'" 
		},
		"photobleaching_constant (cm2/(J*M))": {
			"value": 600,
			"reference": "Dysart et al., 2005, 'Calculation of Singlet Oxygen Dose ... and Photobleaching During mTHPC Photodynamic Therapy of MLL Cells'"
		},
		"dimensions": {
			"shape": "disc",
			"length_A": 32.8,
			"width_A": 32.8,
			"depth_A": 1.5,
			"notes": "The depth is atomic thickness, as quantified by this paper https://www.nature.com/articles/ncomms1291."
		}
	}


- *e_quantum_yield*, *so_specificity*, & *so_quantum_yield* ``float``: quantum yields of the photosensitizer (PS). The ``e_quantum_yield`` and ``so_specificity`` quantum yields correspond to the excitation of a PS after collision with a photon and the generation of singlet oxygen (SO) by that excited PS, respectively. The ``so_quantum_yield`` quatum yield is the collective probability, and thus the product, of the two aforementioned quantum yields, and is disaggregated into estimates of the first two quantum yields by simply taking its square-root.
- *formula* ``str``: defines the chemical formula of the PS, where underscores may optionally separate elements to improve readability.
- *excitation_nm* ``2D-array``: specifies the excitation range(s) of the photosensitizer in nanometers. 
- *ps_rise (fs)* & *ps_decay (ns)* ``float``: specifies the excitation and fluorescence times in units of femtoseconds and nanoseconds, respectively, for the PS.
- *ps_charge_transfer (ns)* ``float``: specifies the average time after a photon collision when the excited PS relays excitation energy to another molecule, ideally tripley oxygen.  
- *photobleaching_constant (cm2/(J*M))* ``float``: specifies the rate constant of the oxygen-dependent photobleaching reaction for the PS.
- *dimensions* ``dict``: specifies the physical dimensions of the PS molecule in angstroms and the approximate shape of the simulated PS, which are used in the absence of absorption data for the PS at the simulated concentration to calculate the molecular volume and then estimate the absorption properties of the chemical in the system.

The ``reference`` and ``notes`` keys are optional, and other keys of metadata may be added as well at the user's discretion.


wells
+++++++++

The dimensions of the solution in which the PDI simulation is conducted may be described through the ``wells.json`` file, in addition to the function argument. The default entries of this file reflect standard dimensions of individual wells in the 6-, 12-, 24-, 48-, and 96-well bioassay plates; however, any solution can be defined that provides the ``area_sqr_cm``, ``depth_cm``, and ``extinction_coefficient (1/m)`` parameters. 

.. code-block:: json

 {
	"12": {
		"area_sqr_cm": 3.85,
		"depth_cm": 1.766,
		"extinction_coefficient (1/m)":0.013,
		"dimensions_reference":"https://ca.vwr.com/assetsvc/asset/en_CA/id/25423331/contents/vwr-essential-products-for-tissue-culture.pdf",
                "coefficient_reference":"The effect of photocarrier generating light on light scattering in the Sea. Lorenzen, 1972"
	}
 }

Other key:value pairs may be defined to specify references or other notes about the system. A sample entry in the ``wells.json`` file is provided below:

- *area_sqr_cm* ``float``: specifies the area on the bottom of the simulated solution, in units of cm\ :sup:`2`\.
- *depth_cm* ``float``: specifies the depth (height) of the simulated solution, in units of cm.
- *extinction_coefficient (1/m)* ``float``: specifies the rate constant for the scattering of light through the solution, as a function of depth, via the light attenuation equation: remaining fraction = e\ :sup:`(-k*z)`\.

bacteria
++++++++++

This folder contains a different JSON file for each bacterial specie. The only bacterium that is specified by default is *Staphylococcus aureus* (``S_aureus.json``), however, other organisms can be defined by replicating the structure of the default parameter file. The values of each sub-dictionary are stored in the ``value`` key.

.. code-block:: json

 {
    "membrane_chemicals": {
      "BC_SFA": {
        "density_gL": {
          "value": 0.9,
          "reference": ["https://pubchem.ncbi.nlm.nih.gov/compound/Stearic-acid#section=Density", "https://pubchem.ncbi.nlm.nih.gov/compound/445639"],
	  "notes": "The density is estimated to be between stearic acid and oleic acid"
        },
        "formula": ["C18_H34_O2","C16_H30_O2"],
	"proportion": {
	  "value": 0.662,
	  "reference": "A. G . O’DONNELL, M. R . NAHAIE, M. GOODFELLOW, D. E. MINNIKINI, and V . HAJEK. Numerical Analysis of Fatty Acid Profiles in the Identification of Staphylococci. Journal of General Microbiology (1989). 131, 2023-2033. https://doi.org/10.1099/00221287-131-8-2023",
	  "notes": "All BCFAs were summed from Table 2 for all S. aureus entries."
	}
      },
      "SC_SFA": {
        "density_gL": {
          "value": 0.94,
          "reference": ["https://pubchem.ncbi.nlm.nih.gov/compound/Stearic-acid#section=Density"],
	  "notes": "The density for all saturated fatty acids is estimated as stearic acid."
        },
        "formula": ["C20_H38_O2","C18_H34_O2","C16_H30_O2"],
	"proportion": {
	  "value": 0.268,
	  "reference": "A. G. O’DONNELL, M. R. NAHAIE, M. GOODFELLOW, D. E. MINNIKINI, and V. HAJEK. Numerical Analysis of Fatty Acid Profiles in the Identification of Staphylococci. Journal of General Microbiology (1989). 131, 2023-2033. https://doi.org/10.1099/00221287-131-8-2023",
	  "notes": "All saturated SCFAs were summed from Table 2 for all S. aureus entries."
    	}
      },
      "SC_UFA": {
        "density_gL": {
          "value": 0.94,
          "reference": ["https://pubchem.ncbi.nlm.nih.gov/compound/Stearic-acid#section=Density"],
	  "notes": "The density for all saturated fatty acids is estimated as stearic acid."
        },
        "formula": ["C20_H38_O2","C18_H34_O2","C16_H30_O2"],
	"proportion": {
	  "value": 0.07,
	  "reference": "A. G. O’DONNELL, M. R. NAHAIE, M. GOODFELLOW, D. E. MINNIKINI, and V. HAJEK. Numerical Analysis of Fatty Acid Profiles in the Identification of Staphylococci. Journal of General Microbiology (1989). 131, 2023-2033. https://doi.org/10.1099/00221287-131-8-2023",
	  "notes": "All saturated SCFAs were summed from Table 2 for all S. aureus entries."
	}
      }
    },
    "membrane_thickness_nm": {
      "value": 4,
      "reference": "W.Rawicz, K.C.Olbrich, T.McIntosh, D.Needham, E.Evans (2000). Effect of Chain Length and Unsaturation on Elasticity of Lipid Bilayers. Biophysical Journal, 79(1), 328–339. https://doi.org/10.1016/S0006-3495(00)76295-3  ; “The electrical capacity of suspensions with special reference to blood” by Fricke, 1925"
    },
    "cell_mass_pg": {
      "value": 1.048,
      "reference": "Lewis, C. L., Craig, C. C., & Senecal, A. G. (2014). Mass and density measurements of live and dead Gram-negative and Gram-positive bacterial populations. Applied and environmental microbiology, 80(12), 3622–3631. https://doi.org/10.1128/AEM.00117-14"
    },
    "cell_volume_fL": {
      "value": 0.9357,
      "reference":  "Lewis, C. L., Craig, C. C., & Senecal, A. G. (2014). Mass and density measurements of live and dead Gram-negative and Gram-positive bacterial populations. Applied and environmental microbiology, 80(12), 3622–3631. https://doi.org/10.1128/AEM.00117-14"
  },
    "eps_oxidation_rate_constant":{
      "value": 37.75,
      "reference": null,
      "notes": "This rate constant was empirically determined after calibrating the predictions with the Beirao et al., 2014 paper that constituted one of our examples"
  },
    "cellular_dry_mass_proportion_biofilm":{
	  "value": 0.1,
	  "reference": "The biofilm matrix; Flemming et al.; 2010"
  },
    "doubling_rate_constant":{
	  "value": 0.00038441,
	  "reference": "Baines et al.; mBio; 2015"
  },
    "biofilm_oxidation_fraction_lysis":{
	  "value": 0.0014125,
	  "note": "empirically derived through training with Beirao et al."
  }
 }

- *membrane_chemicals* ``dict``: specifies the fatty acid constituent the phospholipids of the bacterial cytoplasmic membrane. Each fatty acid (FA) entry is defined with sub-dictionaries of its chemical ``formula``, its ``density_gL`` density, and its ``proportion`` of total FAs in the cytoplasmic membrane. 
- *membrane_thickness_nm* ``dict``: specifies the thickness of the cytoplasmic membrane in nanometers.
- *cell_mass_pg* & *cell_volume_fL* ``dict``: specifies the mass and volume of the bacterial cell in picograms and femtoliters, respectively.
- *eps_oxidation_rate_constant* ``dict``: defines the rate constant of oxidizing the extracellular polymeric substance for biofilm simulations of this organism.
- *cellular_dry_mass_proportion_biofilm* ``dict``: defines the ratio of biofilm mass that is comprised of cellular dry mass.  
- *doubling_rate_constant* ``dict``: specifies the rate constant at which the bacteria doubles
- *biofilm_oxidation_fraction_lysis* ``dict``: specifies the threshold of membrane oxidation that manifests in cellular lysis.

The ``reference`` and ``notes`` keys are optional in each sub-dictionary, as are other keys of metadata at the user's discretion.

light
+++++++
 
The emission spectrum and visual intensity per energy input of the simulated light source are defined through the ``light.json`` file, complementarily to the function argument. The default light sources are ``incandescent``, ``LED``, and ``fluorescent``; however, other light sources can be defined by emulating the same syntactic structure. The values of each sub-dictionary are stored in the ``value`` key.
      
.. code-block:: json
      
 {
  "incandescent": {
    "visible_proportion": {
      "value": 0.1,
      "reference": "Macisaac et al., 1999"
    },
    "lumens_per_watt": {
      "value": 3,
      "reference": "Michael F. Hordeski. Dictionary Of Energy Efficiency Technologies. Fairmont Press. ISBN: 9780824748104"
    }
  }
 }
	  
- *visible_proportion* ``dict``: specifies proportion of the emission spectrum that resides within the visible region.
- *lumens_per_watt* ``dict``: specifies the visual lumens that are emitted per watt of energy. This is used to convert between parameterized light intensity in units of lux or lumens into watts, which is the necessary unit for subsequent PDIpy calculations.