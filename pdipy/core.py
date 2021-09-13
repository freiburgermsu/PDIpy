# import libraries
from scipy.constants import femto, pico, nano, micro, milli, kilo, liter, N_A, h, c, minute
from math import pow, pi, cos, radians
from chemicals.elements import periodic_table
from datetime import date
import json
import sys
import re
import os

class PDIBacterialPkg():
    def __init__(self):
        """ Define the initial PDI and biological parameters
        """        
        # define base organizations
        self.parameters = {}
        self.variables = {}
        self.results = {}
        
        # initial parameters
        self.parameters['time'] = 0
        self.parameters['singlet_oxygen_diffusion_distance'] = 80 * nano       # Moan1990 
        self.parameters['oxidation_angle'] = 5           # degrees
        self.parameters['cwd'] = re.sub('(?!\\\\)(\w+\.py)', '', __file__)

    def define_bacterium(self, bacterial_specie):
        """ Define the model bacteria PDI and is biological parameters
            Used:    define()
        """        
        self.results['cellular_vitality'] = True
        self.results['biofilm_growth'] = True
        
        # load the bacterial parameters
        self.parameters['bacterial_specie'] = bacterial_specie
        bacterium = json.load(open('{}/parameters/{}.json'.format(self.parameters['cwd'], bacterial_specie)))[bacterial_specie]

        self.membrane_chemicals = bacterium['membrane_chemicals']
        self.cell_shape = bacterium['shape']['value']
        self.parameters['stop_biofilm_threshold'] = bacterium['stop_biomass_threshold (%)']['value']
        self.parameters['oxidized_death_threshold'] = bacterium['death_threshold (%)']['value']
        self.membrane_thickness = bacterium['membrane_thickness (nm)']['value'] * nano                    
        self.cell_mass = bacterium['cell_mass (pg)']['value'] * pico
        self.cell_volume = bacterium['cell_volume (fL)']['value'] * femto  / liter


    def define_porphyrin(self, porphyrin_conc, porphyrin_selection = 'A3B_4Zn'):
        """ Define the porphyrin PDI 
            Used:     define()
        """ 
        # load the photosensitizer parameters
        self.parameters['porphyrin_selection'] = porphyrin_selection
        self.photosensitizer = json.load(open('{}/parameters/photosensitizers.json'.format(self.parameters['cwd'])))[porphyrin_selection]
        
        self.parameters['soret'] = {'upper': self.photosensitizer['soret (nm)']['value'][1] * nano, 'lower': self.photosensitizer['soret (nm)']['value'][0] * nano}
        self.parameters['q'] = {'upper': self.photosensitizer['q (nm)']['value'][1] * nano, 'lower': self.photosensitizer['q (nm)']['value'][0] * nano}

        self.variables['porphyrin_conc'] = float(porphyrin_conc)
        

    def define_light(self, light_source, irradiance = None, lux = None, lumens = None, wattage = None, distance = None, reflection = False):
        """ Define the PDI light source
            Used:    define()
        """ 
        self.parameters['visible'] = {'upper': 780 * nano, 'lower': 390 * nano}

        light_parameters = json.load(open('{}/parameters/light_source.json'.format(self.parameters['cwd'])))
        self.parameters['visible_proportion'] = light_parameters[light_source]['visible_proportion']['value']
        
        if irradiance is not None:
            self.parameters['irradiance'] = irradiance
            return 'irradiance'
        
        # define the available light
        if lux is not None:
            # convert lux into watts at the plate
            self.parameters['lux'] = lux
            return 'lux'
        
        if lumens is not None:
            # convert umens into watts at the plate
            self.parameters['lumens = lumens']
            return 'lumens'
            
        else:
            # calculate the fraction of photons that strike the plate
            self.parameters['light_watts'] = wattage


    def singlet_oxygen_calculations(self, timestep, kinetic_constant, photon_collision_proportion, healing_kinetics = 5):
        """ Calculate the intermediates and values that yield the [singlet_oxygen] from PDI
            Used:    simulate()
        """
        try:
            if self.parameters['defined']:
                pass
        except:
            import sys
            sys.exit('ERROR: The model must first be defined.')
            
        self.parameters['timestep'] = timestep
        seconds_per_timestep = self.parameters['timestep'] * minute
            
        effective_visible_light_watts = self.parameters['light_watts'] * self.parameters['visible_proportion']
        visible_region = self.parameters['visible']['upper'] - self.parameters['visible']['lower']
        excitation_visible_proportion = ((self.parameters['q']['upper'] - self.parameters['q']['lower']) + (self.parameters['soret']['upper'] - self.parameters['soret']['lower'])) / visible_region
        effective_excitation_light_watts = excitation_visible_proportion * effective_visible_light_watts  # homogeneous light intesity throughout the visible spectrum is assumed
        
        # photonic calculations
        average_excitation_wavelength = (self.parameters['q']['upper'] + self.parameters['soret']['lower']) / 2
        joules_per_photon = (h * c) / average_excitation_wavelength
        photons_per_second = effective_excitation_light_watts / joules_per_photon
        self.variables['photons_per_timestep'] = photons_per_second * seconds_per_timestep

        # singlet oxygen calculations
        '''so_from_light = photons_per_second * molecules_dissolved_oxygen * excitation_constant'''
        mw_molecular_oxygen = periodic_table.O.MW * 2 * kilo          #mg / mole
        dissolved_oxygen_concentration = 9              # mg / L, ambient water quality criteria for DO, EPA         # this must be adjusted to only consider oxygen in the water in the vacinity of the photosensitizer, which is geometrically limited to the material surface. The continuum assumption of the aqueous solution may be implemented such that only a oxygen within fractional volume of the total solution volume.       
        self.variables['molecules_dissolved_oxygen'] = dissolved_oxygen_concentration / mw_molecular_oxygen
        
        # define the kinetic parameters
        self.variables['quantum_yield'] = self.photosensitizer['quantum_yield']['value'] * self.photosensitizer['so_specificity']['value']
        self.variables['photon_collisions'] = photon_collision_proportion
        self.variables['healing'] = healing_kinetics
        self.variables['k'] = kinetic_constant
        
    def geometric_oxidation(self):
        if self.cell_shape == "sphere":
            # define calculation functions
            def shell_volume(radi_1, radi_2, coeff = 1):
                volume = (4*pi/3) * coeff * (radi_1**3 - radi_2**3)
                return volume
            def cap_volume(r, h):
                volume = (pi/3) * h**2 * (3*r - h)
                return volume
            
            # calculate the cellular dimensions
            cell_radius = pow((self.cell_volume * 3) / (4 * pi), 1/3)
            membrane_inner_radius = cell_radius - self.membrane_thickness

            # calculate the cellular and oxidation volumes
            membrane_volume = shell_volume(cell_radius, membrane_inner_radius) # M^3
            outer_h = cell_radius * (1 - cos(radians(self.parameters['oxidation_angle'])))
            cap_volume = cap_volume(cell_radius, outer_h)
            oxidized_membrane_volume_ratio = cap_volume / membrane_volume
#             oxidized_membrane_volume_ratio = self.parameters['oxidation_angle'] / 360
            print('oxidized volume proportion: ', oxidized_membrane_volume_ratio)

            #calculate the cellular and oxidation areas
            oxidized_cap_area = 2 * pi * cell_radius * outer_h
            cell_area = 4 * pi * cell_radius ** 2
            oxidized_area_ratio = oxidized_cap_area / cell_area
            print('oxidized area proportion: ',oxidized_area_ratio, '\n\n')

            # fatty acid oxidation
            singlet_oxygen_interaction_radius = cell_radius + self.parameters['singlet_oxygen_diffusion_distance']
#             membrane_solution_interface_volume = shell_volume(singlet_oxygen_interaction_radius, cell_radius)

            c17_oxidized_volume = membrane_volume * self.membrane_chemicals['anteiso_C17']['proportion']['value'] * oxidized_membrane_volume_ratio  # M^3 
            c17_oxidized_ppm = c17_oxidized_volume * self.membrane_chemicals['anteiso_C17']['density (g/L)']['value'] * milli / liter
            c15_oxidized_volume = membrane_volume * self.membrane_chemicals['anteiso_C15']['proportion']['value'] * oxidized_membrane_volume_ratio  # M^3
            c15_oxidized_ppm = c15_oxidized_volume * self.membrane_chemicals['anteiso_C15']['density (g/L)']['value'] * milli / liter

            self.variables['fa17_conc'] = c17_oxidized_ppm  / self.membrane_chemicals['anteiso_C17']['mw']
            self.variables['fa15_conc'] = c15_oxidized_ppm / self.membrane_chemicals['anteiso_C15']['mw']


    def kinetic_calculation(self, end_time, omex_file_path, omex_file_name = None):
        """ Execute the kinetic calculations in Tellurium
            Used:    simulate()
        """
        import tellurium
        
        # define the first equation
        k_so = self.variables['quantum_yield'] * self.variables['photon_collisions'] * self.variables['photons_per_timestep'] * self.variables['porphyrin_conc']
        mo = self.variables['molecules_dissolved_oxygen']
        
        # define the second equation
        k = self.variables['k']
        fa17 = self.variables['fa17_conc']
        fa15 = self.variables['fa15_conc']
        
        # define constants
        healing_kinetics = self.variables['healing']
        biofilm_threshold = self.parameters['stop_biofilm_threshold']
        death_threshold = self.parameters['oxidized_death_threshold']

        # define the SBML model
        model = (f'''
          model pdi_oxidation
            # expressions
            o -> so;  {k_so}*o
            so + fa17 -> ofa; {k}*so*fa17 - {healing_kinetics}*ofa  # time      # the aggregated photons / second must be programmatically inserted into the rate expression         
            so + fa15 -> ofa; {k}*so*fa15 - {healing_kinetics}*ofa

            # define the first expression 
            o = {mo};
            
            # define the second expression
            so = 0;
            fa17 = {fa17};
            fa15 = {fa15};
            
            # define constants
            biofilm = {biofilm_threshold};
            vitality = {death_threshold};
            oxidation := ofa / (ofa + fa15 + fa17);
            
          end
        ''')
        tellurium_model = tellurium.loada(model)
        print(tellurium_model.getCurrentAntimony())
        print('\nCurrent integrator:', '\n', tellurium_model.integrator)
        
        
        # define the SEDML plot
        total_points = end_time / self.parameters['timestep'] 
        phrasedml_str = f'''
          model1 = model "pdi_oxidation"
          sim1 = simulate uniform(0, {end_time}, {total_points})
          task1 = run sim1 on model1
          plot "Figure 1" time vs biofilm, vitality, oxidation, so, fa17, fa15, ofa
        '''

        # create, execute, and export an OMEX file
        inline_omex = '\n'.join([model, phrasedml_str])               
        tellurium.executeInlineOmex(inline_omex)
        
        if omex_file_name is None:
            count = 0
            omex_file_name = '_'.join([str(date.today()), self.parameters['porphyrin_selection'], self.parameters['bacterial_specie'], str(count)])
            while os.path.exists(omex_file_name):
                count += 1
                omex_file_name = '_'.join([str(date.today()), self.parameters['porphyrin_selection'], self.parameters['bacterial_specie'], str(count)])
            omex_file_name += '.omex'
        tellurium.exportInlineOmex(inline_omex, os.path.join(omex_file_path, omex_file_name))
              

    def define(self, bacterial_species, porphyrin_selection, porphyrin_conc, light_source, wattage):
        """ Parameterize the model 
            Arguments (type, units):
                bacterial_species (string, \capitalGenus\fullLowerSpecies) = the selected bacterial species for study 
                porphyrin_selection (string, grantCode) = the selected porphyrin for study
                porphyrin_conc (float, mg/L) = the concentration of the selected porphyrin in the aqueous environment
                light_source (string, lightName) = the light source for the study
                wattage (float, joules/second) = the wattage of the light source 
        """
        # parameterize the simulation
        self.define_bacterium(bacterial_species)
        self.define_porphyrin(porphyrin_conc, porphyrin_selection)
        light = self.define_light(light_source, wattage = wattage)
        
        self.parameters['defined'] = True
        
        return light
        
    def simulate(self, light, timestep, kinetc_constant, photon_collision_proportion, omex_file_path, end_time):
        ''' Execute the model for this simulation
            Arguments (type, units):
                end_time (float, minutes) = the conclusion time for the simulation in minutes
                timestep (float, minutes) = the loop time for the simulation in minutes
                kinetic_constant (float, ___) = the kinetic constant for the oxidation of fatty acids via singlet oxygen
        '''
        # execute the simulation
        if light == 'irradiance':
            pass
        else:
            self.singlet_oxygen_calculations(timestep, kinetc_constant, photon_collision_proportion)
            self.geometric_oxidation()
            self.kinetic_calculation(end_time, omex_file_path)