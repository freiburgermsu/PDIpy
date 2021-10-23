define 
""" Parameterize the model 
            Arguments (type, units):
                bacterial_species (string, \capitalGenus\fullLowerSpecies) = the selected bacterial species for study 
                porphyrin_selection (string, grantCode) = the selected porphyrin for study
                porphyrin_conc (float, mg/L) = the concentration of the selected porphyrin in the aqueous environment
                light_source (string, lightName) = the light source for the study
                wattage (float, joules/second) = the wattage of the light source 
"""

time_to_threshold 
''' Execute the model for this simulation
            Arguments (type, units):
                end_time (float, minutes) = the conclusion time for the simulation in minutes
                timestep (float, minutes) = the loop time for the simulation in minutes
                kinetic_constant (float, ___) = the kinetic constant for the oxidation of fatty acids via singlet oxygen
'''

kinetic_calculation 
""" Execute the kinetic calculations in Tellurium
            Used:    simulate()
"""       

singlet_oxygen_calculations
""" Calculate the intermediates and values that yield the [singlet_oxygen] from PDI
            Used:    simulate()
"""

define_bacterium
""" Define the model bacteria PDI and is biological parameters
            Used:    define()
"""    

define_photosensitizer
""" Define the porphyrin PDI 
            Used:     define()
""" 

define_light
""" Define the PDI light source
            Used:    define()
""" 