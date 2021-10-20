# import libraries
from scipy.constants import femto, pico, angstrom, nano, micro, milli, centi, kilo, liter, N_A, h, c, minute
from math import pow, pi, sin, cos, radians
from chemicals.elements import periodic_table
from to_precision import sci_notation
from matplotlib import pyplot
from datetime import date
import tellurium
import pandas
import json
import sys
import re
import os

pandas.set_option('max_colwidth', None)


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
    def __init__(self, solution_volume = 0.75, verbose = False):
        """ Define the initial PDI and biological parameters
        """        
        # define base organizations
        self.parameters = {}
        self.variables = {}
        self.results = {}
        
        # initial parameters
        self.parameters['singlet_oxygen_diffusion_distance'] = 80 * nano       # Moan1990 
        self.parameters['oxidation_angle'] = 5           # degrees
        self.parameters['root_path'] = re.sub('(\w+\.py)', '', __file__)
        self.parameters['well_solution_volume'] = solution_volume * micro * liter
        self.verbose = verbose

    def define_bacterium(self, bacterial_specie):
        """ Define the model bacteria PDI and is biological parameters
            Used:    define()
        """        
        self.results['cellular_vitality'] = True
        self.results['biofilm_growth'] = True
        
        # load the bacterial parameters
        self.parameters['bacterial_specie'] = bacterial_specie
        bacterium = json.load(open('{}/parameters/{}.json'.format(self.parameters['root_path'], bacterial_specie)))[bacterial_specie]

        self.membrane_chemicals = bacterium['membrane_chemicals']
        self.cell_shape = bacterium['shape']['value']
        self.parameters['stop_biofilm_threshold'] = bacterium['stop_biomass_threshold (%)']['value']
        self.parameters['oxidized_death_threshold'] = bacterium['death_threshold (%)']['value']
        self.membrane_thickness = bacterium['membrane_thickness (nm)']['value'] * nano        # meters                 
        self.cell_mass = bacterium['cell_mass (pg)']['value'] * pico                          # grams
        self.cell_volume = bacterium['cell_volume (fL)']['value'] * femto  / liter            #cubic meters


    def define_porphyrin(self, porphyrin_conc, porphyrin_selection = 'A3B_4Zn'):
        """ Define the porphyrin PDI 
            Used:     define()
        """ 
        # load the photosensitizer parameters
        self.parameters['porphyrin_selection'] = porphyrin_selection
        self.photosensitizer = json.load(open('{}/parameters/photosensitizers.json'.format(self.parameters['root_path'])))[porphyrin_selection]
        
        self.parameters['soret'] = {'upper': self.photosensitizer['soret (nm)']['value'][1] * nano, 'lower': self.photosensitizer['soret (nm)']['value'][0] * nano}
        self.parameters['q'] = {'upper': self.photosensitizer['q (nm)']['value'][1] * nano, 'lower': self.photosensitizer['q (nm)']['value'][0] * nano}

        self.variables['porphyrin_ppm'] = porphyrin_conc
        self.variables['porphyrin_moles'] = porphyrin_conc / self.photosensitizer['mw']['value'] * self.parameters['well_solution_volume']/liter

    def define_light(self, light_source, irradiance = None, well_surface_area = 1, exposure = None, simulation_time = None, lux = None, lumens = None, wattage = None, distance = None, reflection = False):
        """ Define the PDI light source
            Used:    define()
        """ 
        self.parameters['visible'] = {'upper': 780 * nano, 'lower': 390 * nano}

        light_parameters = json.load(open('{}/parameters/light_source.json'.format(self.parameters['root_path'])))
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


    def define_photosensitizer_volume(self, molecular_proportion = None, porphyrin_per_square_cm = None):
        if molecular_proportion is not None:
            self.parameters['molecular_proportion'] = molecular_proportion
        else:
            # determine the individual molecular components
            if self.parameters['porphyrin_selection'] == 'A3B_4Zn':
                # calculate the volume of photosensitizing region in the molecule
                center_porphyrin_length = sci_notation(2*(chemical_dimensions['bond']['c-c']*(2*cos(radians(chemical_dimensions['angle']['sp2']-90))+cos(radians(180-chemical_dimensions['angle']['sp2'])))), sigfigs)
                sp2_extension = sci_notation(chemical_dimensions['bond']['c-c'] * (2 + cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-n'] * cos(radians(180-chemical_dimensions['angle']['sp2'])) + chemical_dimensions['bond']['c-n'], sigfigs)
                sp3_diazirine = 0
#                 sp3_diazirine = sci_notation(chemical_dimensions['bond']['c-c']*cos(radians(chemical_dimensions['angle']['sp3']-90)) + 2*chemical_dimensions['bond']['c-c']*(1+cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-c']*sin(radians(chemical_dimensions['angle']['sp3']-90)) + chemical_dimensions['bond']['c-f'], sigfigs)
#                 min_thickness = 2*chemical_dimensions['bond']['c-f']*sin(radians(chemical_dimensions['angle']['sp3']))
#                 max_thickness = 2*float(sp3_diazirine)
#                 print('min', min_thickness, 'max', max_thickness)
#                 thickness = sci_notation(average(min_thickness, max_thickness), sigfigs)

                diazirine_length = float(sp2_extension) + float(sp3_diazirine)
                length = sci_notation(float(center_porphyrin_length) + 2*diazirine_length, sigfigs)
                thickness = 0.5 * angstrom
                
                # calculate the proportion of well volume that is constituted by the molecular volume, in meters
                molecular_volume = sci_notation((float(length) * angstrom)**2 * thickness, sigfigs)
                num_photosensitizers_well = (milli*N_A)*(liter*milli)/96                  # 1 millimolar, with one 1 mL, per one of 96 wells
                self.variables['photon_collision_proportion'] = num_photosensitizers_well*float(molecular_volume) / self.parameters['well_solution_volume']
                
                # calculate the layer distribution for a slice 
                if porphyrin_per_square_cm:
                    self.parameters['porphyrin_per_square_cm'] = porphyrin_per_square_cm
                    porphyrin_layer_area = porphyrin_layer_area * centi**2
                else:
                    solution_depth = 1*centi
                    layers = solution_depth / thickness
                    porphyrin_per_layer = self.variables['porphyrin_ppm'] / layers
                    solution_incident_area = self.parameters['well_solution_volume'] / solution_depth
                    orthogonal_area = pi*(18*angstrom)**2
                    parallel_area = (18*angstrom) * thickness
                    average_area = average(orthogonal_area, parallel_area)
                    porphyrin_layer_area = porphyrin_per_layer * average_area
                    
                # total porphyrin area proportion
                self.variables['area_proportion'] = porphyrin_layer_area*self.variables['porphyrin_moles']*N_A/solution_incident_area
                
                if self.verbose:
                    print(f'The center porphyrin object is {center_porphyrin_length} angstroms')
                    print(f'The benzyl extension is {sp2_extension} angstroms')
                    print(f'The diazirine is {sp3_diazirine} angstroms')
                    print(f'The molecular length is {length} angstroms')
                    print(f'The molecular thickness is {thickness} angstroms')
                    print(f'The molecular volume is {molecular_volume} cubic meters')
                    print('The photosensitizer volume proportion is {}'.format(self.variables['photon_collision_proportion']))
                    print('The photosensitizer area proportion is {}'.format(self.variables['area_proportion']))
            else:
                print('--> ERROR: The {} porphyrin selection is not defined.'.format(self.parameters['porphyrin_selection']))
                
        self.defined_model = True


    def singlet_oxygen_calculations(self, timestep, total_time, kinetic_constant, initial_time = 0, healing_kinetics = 5):
        """ Calculate the intermediates and values that yield the [singlet_oxygen] from PDI
            Used:    simulate()
        """
        if not self.defined_model:
            import sys
            sys.exit('ERROR: The model must first be defined.')
            
        self.parameters['initial_time'] = initial_time
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
            
            self.variables['dissolved_mo_moles'] = dissolved_oxygen_concentration / mw_molecular_oxygen * self.parameters['well_solution_volume'] * N_A
            estimated_excited_photosensitizers = photons_per_second * self.parameters['total_time'] * self.variables['photon_collision_proportion']
            
            if self.verbose:
                print('molecular oxygen molecules: ', sci_notation(self.variables['dissolved_mo_moles'], sigfigs))
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

            if 'bcfa' in self.membrane_chemicals:
                self.variables['bcfa_conc'] = self.membrane_chemicals['bcfa']['concentration']
            else:
                for chemical in self.membrane_chemicals:
                    if re.search('anteiso', chemical):
                        oxidized_volume = membrane_volume * self.membrane_chemicals[chemical]['proportion']['value'] * oxidized_membrane_volume_ratio  # M^3 
                        oxidized_g = oxidized_volume/liter * self.membrane_chemicals[chemical]['density (g/L)']['value'] 
                        self.variables['bcfa_moles'] = oxidized_g / self.membrane_chemicals[chemical]['mw']
                        self.variables['bcfa_conc'] = self.membrane_chemicals[chemical]['density (g/L)']['value'] / self.membrane_chemicals[chemical]['mw']


    def kinetic_calculation(self):
        """ Execute the kinetic calculations in Tellurium
            Used:    simulate()
        """        
        # define the first equation
        k_so = self.variables['quantum_yield'] * self.variables['photons_per_timestep'] * self.variables['porphyrin_moles'] * self.variables['photon_collision_proportion'] # 0.0000777859848166666
        mo = self.variables['dissolved_mo_moles']
        
        # define the second equation
        k = self.variables['k']
        bcfa = self.variables['bcfa_conc']
        
        # define constants
        healing_kinetics = self.variables['healing']
        biofilm_threshold = self.parameters['stop_biofilm_threshold']
        death_threshold = self.parameters['oxidized_death_threshold']

        # define the SBML model
        self.model = (f'''
          model pdi_oxidation
            # expressions, chemical quantities in lieu of concentrations to accommodate different compartment volumes 
            mo -> so;  {k_so}*mo
            so + bcfa -> ofa; {k}*so*bcfa - {healing_kinetics}*ofa    # the aggregated photons / second must be programmatically inserted into the rate expression         

            # define the first expression 
            mo = {mo}
            so = 0
            
            # define the second expression
            ofa = 0
            bcfa = {bcfa}
            
            # define constants
            biofilm = {biofilm_threshold};
            vitality = {death_threshold};
            oxidation := ofa / (ofa + bcfa);
            
          end
        ''')       
        # define the SEDML plot
        total_points = self.parameters['total_time'] / self.parameters['timestep'] 
        self.phrasedml_str = '''
          model1 = model "pdi_oxidation"
          sim1 = simulate uniform({}, {}, {})
          task1 = run sim1 on model1
          plot "Oxidation proportion of prokaryotic membrane fatty acids" time vs biofilm, vitality, oxidation
        '''.format(self.parameters['initial_time'], self.parameters['total_time'], total_points)
        
        # execute the figure data
        tellurium_model = tellurium.loada(self.model)            
        result = tellurium_model.simulate(self.parameters['initial_time'], int(self.parameters['total_time']), int(total_points))
        self.result_df = pandas.DataFrame(result)
        self.result_df.index = self.result_df[0]
        del self.result_df[0]
        self.result_df.index.name = 'Time (s)'
        self.result_df.columns = ['[o]', '[so]', '[bcfa]', '[ofa]']
        
        if self.verbose:
            print(tellurium_model.getCurrentAntimony())
            print('\nCurrent integrator:', '\n', tellurium_model.integrator)
            print(self.result_df)
            
        return self.result_df
        
    def export(self, simulation_path = None, figure_title = 'Oxidation proportion of prokaryotic membrane fatty acids', x_label = 'Time (s)', y_label = 'oxidation proportion', biofilm_threshold = 0.05, vitality_threshold = 0.1, ):
        # parse the simulation results
        x_values = []
        y_values = []
        for index, point in self.result_df.iterrows():
            x_values.append(index)
            oxidation_proportion = point['[ofa]'] / (point['[ofa]']+point['[bcfa]']) 
            y_values.append(oxidation_proportion)
        biofilm_threshold = [biofilm_threshold for x in range(len(x_values))]
        vitality_threshold = [vitality_threshold for x in range(len(x_values))]
        
        # define the simulation_path
        self.simulation_path = simulation_path
        if self.simulation_path is None:
            count = 0
            self.simulation_path = '_'.join([str(date.today()), self.parameters['porphyrin_selection'], self.parameters['bacterial_specie'], str(count)])
            while os.path.exists(self.simulation_path):
                count += 1
                self.simulation_path = '_'.join([str(date.today()), self.parameters['porphyrin_selection'], self.parameters['bacterial_specie'], str(count)])
            self.simulation_path = os.path.join(os.getcwd(), self.simulation_path)        
        os.mkdir(self.simulation_path)
        
        # create and export the OMEX file
        inline_omex = '\n'.join([self.model, self.phrasedml_str])  
        omex_file_name = 'input.omex'
        tellurium.exportInlineOmex(inline_omex, os.path.join(self.simulation_path, omex_file_name))
        self.result_df.to_csv(os.path.join(self.simulation_path, 'raw_data.csv'))
        
        # export the figure
        figure_path = os.path.join(self.simulation_path, 'output.svg')
        pyplot.rcParams['figure.figsize'] = (15, 9)
        pyplot.rcParams['figure.dpi'] = 150
        figure, ax = pyplot.subplots()
        ax.plot(x_values, y_values, label = 'oxidation_proportion')
        ax.plot(x_values, biofilm_threshold, label = 'biofilm_thershold')
        ax.plot(x_values, vitality_threshold, label = 'vitality_threshold')
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        ax.set_title(figure_title)
        ax.legend()
        figure.savefig(figure_path)
        figure.show()
        
        # define a table of parameters
        parameters = {'parameter':[], 'value':[]}
        parameters['parameter'].append('simulation_path')
        parameters['value'].append(self.simulation_path)
        for parameter in self.parameters:
            parameters['parameter'].append(parameter)
            parameters['value'].append(self.parameters[parameter])
        parameters_table = pandas.DataFrame(parameters)
        if self.verbose:
            print(parameters_table)
        
        parameters_path = os.path.join(self.simulation_path, 'parameters.csv')
        parameters_table.to_csv(parameters_path)
        
        # define a table of variables
        variables = {'variable':[], 'value':[]}
        variables['variable'].append('simulation_path')
        variables['value'].append(self.simulation_path)
        for variable in self.variables:
            variables['variable'].append(variable)
            variables['value'].append(self.variables[variable])
        variables_table = pandas.DataFrame(variables)
        if self.verbose:
            print(variables_table)
        
        variables_path = os.path.join(self.simulation_path, 'variables.csv')
        variables_table.to_csv(variables_path)
             

    def define(self, bacterial_species, porphyrin_selection, porphyrin_conc, light_source, irradiance = None, well_surface_area = 1, exposure = None, simulation_time = None, molecular_proportion = None):
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
        self.define_photosensitizer_volume(molecular_proportion = None)
        
        self.defined_model = True
        
        return light
        
    def simulate(self, timestep, total_time, kinetc_constant, simulation_path = None, figure_title = 'Oxidation proportion of prokaryotic membrane fatty acids', x_label = 'time', y_label = 'oxidation proportion'):
        ''' Execute the model for this simulation
            Arguments (type, units):
                end_time (float, minutes) = the conclusion time for the simulation in minutes
                timestep (float, minutes) = the loop time for the simulation in minutes
                kinetic_constant (float, ___) = the kinetic constant for the oxidation of fatty acids via singlet oxygen
        '''
        # execute the simulation
        if 'watts' in self.parameters:
            self.singlet_oxygen_calculations(timestep, total_time, kinetc_constant)
            self.geometric_oxidation()
            raw_data = self.kinetic_calculation()
            self.export(simulation_path, figure_title, x_label, y_label, )
        else:
            print(f'--> ERROR: The {light} light source is not supported by the code')
            
        return raw_data