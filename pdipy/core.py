# import libraries
from scipy.constants import femto, pico, angstrom, nano, micro, milli, kilo, liter, N_A, h, c, minute
from math import pow, pi, sin, cos, radians
from chemicals.elements import periodic_table
from to_precision import sci_notation
from datetime import date
import json
import sys
import re
import os


def average(num_1, num_2 = None):
    if num_2 is not None:
        numbers = [num_1, num_2]
        average = sum(numbers) / len(numbers)
        return average
    else:
        return num_1

# chemical dimensions in Angstroms (as the averages from https://en.wikipedia.org/wiki/Bond_length) and degrees
chemical_dimensions = {
    'bond':{
        'c-c':average(1.2,1.54),
        'c-h':average(1.06,1.12),
        'c-n':average(1.47,2.1),
        'c-f':average(1.34),
        'n=n':average(1.23) # https://doi.org/10.1016/B978-0-08-101033-4.00003-6
    },
    'angle':{
        'sp3':109.5,
        'sp2':120
    }
}
sigfigs = 2


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
        

    def define_light(self, light_source, irradiance = None, well_surface_area = 1, exposure = None, simulation_time = None, lux = None, lumens = None, wattage = None, distance = None, reflection = False):
        """ Define the PDI light source
            Used:    define()
        """ 
        self.parameters['visible'] = {'upper': 780 * nano, 'lower': 390 * nano}

        light_parameters = json.load(open('{}/parameters/light_source.json'.format(self.parameters['cwd'])))
        self.parameters['visible_proportion'] = light_parameters[light_source]['visible_proportion']['value']
        
        # define the available light
        if irradiance is not None: # mW / cm^2
            self.parameters['watts'] = irradiance * milli * well_surface_area
            return 'watts'
        
        if exposure is not None:  # J / cm^2
            if simulation_time is not None: # minutes
                simulation_time = simulation_time * minute
                self.parameters['watts'] = exposure / simulation_time * well_surface_area
                return 'watts'
            
#             else:
#                 print('--> ERROR: The simulation_time must be defined to calculation singlet oxygen generation')
#                 self.parameters['joules'] = exposure
#                 return 'joules'
        
#         if lux is not None:
#             # convert lux into watts at the plate
#             self.parameters['lux'] = lux
#             return 'lux'
        
#         if lumens is not None:
#             # convert umens into watts at the plate
#             self.parameters['lumens = lumens']
#             return 'lumens'
            
        else:
            print('--> ERROR: The light source has an unrecognized dimension.')
#             # calculate the fraction of photons that strike the plate
#             self.parameters['light_watts'] = wattage


    def define_photosensitizer_volume(self, molecular_proportion = None, verbose = False):
        if molecular_proportion is not None:
            self.parameters['molecular_proportion'] = molecular_proportion
        else:
            # determine the individual molecular components
            if self.parameters['porphyrin_selection'] == 'A3B_4Zn':
                center_porphyrin_length = sci_notation(2*(chemical_dimensions['bond']['c-c']*(2*cos(radians(chemical_dimensions['angle']['sp2']-90))+cos(radians(180-chemical_dimensions['angle']['sp2'])))), sigfigs)

                sp2_extension = sci_notation(chemical_dimensions['bond']['c-c'] * (2 + cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-n'] * cos(radians(180-chemical_dimensions['angle']['sp2'])) + chemical_dimensions['bond']['c-n'], sigfigs)

                sp3_diazirine = sci_notation(chemical_dimensions['bond']['c-c']*cos(radians(chemical_dimensions['angle']['sp3']-90)) + 2*chemical_dimensions['bond']['c-c']*(1+cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-c']*sin(radians(chemical_dimensions['angle']['sp3']-90)) + chemical_dimensions['bond']['c-f'], sigfigs)

                diazirine_length = float(sp2_extension) + float(sp3_diazirine)

                # determine the molecular dimensions
                length = sci_notation(float(center_porphyrin_length) + 2*diazirine_length, sigfigs)
                min_thickness = 2*chemical_dimensions['bond']['c-f']*sin(radians(chemical_dimensions['angle']['sp3']))
                max_thickness = 2*float(sp3_diazirine)
                thickness = sci_notation(average(min_thickness, max_thickness), sigfigs)

                # calculate the proportion of well volume that is constituted by the molecular volume, in meters
                molecular_volume = sci_notation((float(length) * angstrom)**2 * (float(thickness) * angstrom), sigfigs)
                self.parameters['well_solution_volume'] = 0.75 * micro * liter
                num_photosensitizers_well = (milli*N_A)*(liter*milli)/96 # 1 millimolar, with one 1 mL, per one of 96 wells
                self.variables['photon_collision_proportion'] = num_photosensitizers_well*float(molecular_volume) / self.parameters['well_solution_volume']

                if verbose:
                    print(f'The center porphyrin object is {center_porphyrin_length} angstroms')
                    print(f'The benzyl extension is {sp2_extension} angstroms')
                    print(f'The diazirine is {sp3_diazirine} angstroms')
                    print(f'The molecular length is {length} angstroms')
                    print(f'The molecular thickness is {thickness} angstroms')
                    print(f'The molecular volume is {molecular_volume} cubic meters')
                    print('The photosensitizer volume proportion is {proportion}'.format(self.variables['photon_collision_proportion']))
            else:
                print('--> ERROR: The {} porphyrin selection is not defined.'.format(self.parameters['porphyrin_selection']))


    def singlet_oxygen_calculations(self, timestep, total_time, kinetic_constant, healing_kinetics = 5, verbose = True):
        """ Calculate the intermediates and values that yield the [singlet_oxygen] from PDI
            Used:    simulate()
        """
        try:
            if self.parameters['defined']:
                pass
        except:
            import sys
            sys.exit('ERROR: The model must first be defined.')
            
        self.parameters['total_time'] = total_time * minute
        self.parameters['timestep'] = timestep
        seconds_per_timestep = self.parameters['timestep'] * minute            
            
        if 'watts' in self.parameters:                   
            # define the light watts
            effective_visible_light_watts = self.parameters['watts'] * self.parameters['visible_proportion']
            visible_region = self.parameters['visible']['upper'] - self.parameters['visible']['lower']
            excitation_visible_proportion = ((self.parameters['soret']['upper'] - self.parameters['soret']['lower'])) / visible_region
            effective_excitation_light_watts = excitation_visible_proportion * effective_visible_light_watts  # homogeneous light intesity throughout the visible spectrum is assumed
#             print(effective_excitation_light_watts)
            
            # photonic calculations
            average_excitation_wavelength = (self.parameters['q']['upper'] + self.parameters['soret']['lower']) / 2
            joules_per_photon = (h * c) / average_excitation_wavelength
            photons_per_second = effective_excitation_light_watts / joules_per_photon
            self.variables['photons_per_timestep'] = photons_per_second * seconds_per_timestep
            print('photons per timestep: ', self.variables['photons_per_timestep'])
            
            # singlet oxygen calculations
            '''so_from_light = photons_per_second * molecules_dissolved_oxygen * excitation_constant'''
            mw_molecular_oxygen = periodic_table.O.MW*2*kilo          #mg / mole
            dissolved_oxygen_concentration = 9              # mg / L, ambient water quality criteria for DO, EPA         # this must be adjusted for the material surface system to only consider oxygen in the water in the vacinity of the surface material photosensitizer. The continuum assumption of the aqueous solution may be implemented such that only a oxygen within fractional volume of the total solution volume.       
            
            self.variables['molecules_dissolved_oxygen'] = dissolved_oxygen_concentration / mw_molecular_oxygen * self.parameters['well_solution_volume'] * N_A
            estimated_excited_photosensitizers = photons_per_second * self.parameters['total_time'] * self.variables['photon_collision_proportion']
            
            if verbose:
                print('molecular oxygen molecules: ', sci_notation(self.variables['molecules_dissolved_oxygen'], sigfigs))
                print('excited photosensitizer molecules: ', sci_notation(estimated_excited_photosensitizers, sigfigs))

            # define the kinetic parameters
            self.variables['quantum_yield'] = self.photosensitizer['quantum_yield']['value'] * self.photosensitizer['so_specificity']['value']
            self.variables['healing'] = healing_kinetics
            self.variables['k'] = kinetic_constant

        else:
            print('--> ERROR: The singlet oxygen generation cannot be calculated from the light source')
            
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


    def kinetic_calculation(self, omex_file_path, omex_file_name = None):
        """ Execute the kinetic calculations in Tellurium
            Used:    simulate()
        """
        import tellurium
        
        # define the first equation
        k_so = self.variables['quantum_yield'] * self.variables['photons_per_timestep'] * self.variables['porphyrin_conc'] * self.variables['photon_collision_proportion']
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
        total_points = self.parameters['total_time'] / self.parameters['timestep'] 
        phrasedml_str = '''
          model1 = model "pdi_oxidation"
          sim1 = simulate uniform(0, {}, {})
          task1 = run sim1 on model1
          plot "Figure 1" time vs biofilm, vitality, oxidation, so, fa17, fa15, ofa
        '''.format(self.parameters['total_time'], total_points)

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
              

    def define(self, bacterial_species, porphyrin_selection, porphyrin_conc, light_source, irradiance = None, well_surface_area = 1, exposure = None, simulation_time = None, molecular_proportion = None, verbose = False):
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
        light = self.define_light(light_source, irradiance, well_surface_area, exposure, simulation_time)
        self.define_photosensitizer_volume(molecular_proportion = None, verbose = False)
        
        self.parameters['defined'] = True
        
        return light
        
    def simulate(self, light, timestep, kinetc_constant, omex_file_path, total_time):
        ''' Execute the model for this simulation
            Arguments (type, units):
                end_time (float, minutes) = the conclusion time for the simulation in minutes
                timestep (float, minutes) = the loop time for the simulation in minutes
                kinetic_constant (float, ___) = the kinetic constant for the oxidation of fatty acids via singlet oxygen
        '''
        # execute the simulation
        if light == 'watts':
            self.singlet_oxygen_calculations(timestep, total_time, kinetc_constant)
            self.geometric_oxidation()
            self.kinetic_calculation(omex_file_path)
        else:
            print(f'--> ERROR: The {light} light source is not supported by the code')