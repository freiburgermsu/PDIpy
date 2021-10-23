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
    def __init__(self, solution_volume = 0.75, verbose = False, jupyter = False):
        """ Define the initial PDI and biological parameters
        """        
        # define base organizations
        self.parameters = {}
        self.variables = {}
        self.results = {}
        
        # initial parameters
        self.parameters['singlet_oxygen_diffusion_distance'] = 80 * nano       # Moan1990 
        self.parameters['oxidation_angle'] = 5           # degrees
        self.parameters['root_path'] = re.sub('(\\\\.+\.py)', '', __file__)
        self.parameters['well_solution_volume'] = solution_volume * micro * liter
        self.verbose = verbose

    def define_bacterium(self, bacterial_specie):    
        self.results['cellular_vitality'] = True
        self.results['biofilm_growth'] = True
        
        # load the bacterial parameters
        self.parameters['bacterial_specie'] = bacterial_specie
        self.bacterium = json.load(open('{}/parameters/{}.json'.format(self.parameters['root_path'], bacterial_specie)))[bacterial_specie]

        self.membrane_chemicals = self.bacterium['membrane_chemicals']
        self.parameters['stop_biofilm_threshold'] = self.bacterium['stop_biomass_threshold (%)']['value']
        self.parameters['oxidized_death_threshold'] = self.bacterium['death_threshold (%)']['value']

    def define_photosensitizer(self, photosensitizer_molar, photosensitizer = 'A3B_4Zn'):
        # load the photosensitizer parameters
        self.parameters['photosensitizer_selection'] = photosensitizer
        self.photosensitizer = json.load(open('{}/parameters/photosensitizers.json'.format(self.parameters['root_path'])))[photosensitizer]
        
        self.parameters['soret'] = {'upper': self.photosensitizer['soret (nm)']['value'][1] * nano, 'lower': self.photosensitizer['soret (nm)']['value'][0] * nano}
        self.parameters['q'] = {'upper': self.photosensitizer['q (nm)']['value'][1] * nano, 'lower': self.photosensitizer['q (nm)']['value'][0] * nano}

        self.parameters['photosensitizer_molar'] = photosensitizer_molar
        self.variables['photosensitizers'] = (self.parameters['photosensitizer_molar']*N_A) * (self.parameters['well_solution_volume']/liter)

    def define_light(self, light_source, irradiance = None, well_surface_area = 1, exposure = None, simulation_time = None, lux = None, lumens = None, wattage = None, distance = None, reflection = False):
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
            # calculate the fraction of photons that strike the plate
#             self.parameters['light_watts'] = wattage


    def define_photosensitizer_volume(self, molecular_proportion = None, photosensitizer_moles_per_square_cm = None, solution_depth = 1):
        if molecular_proportion is not None:
            self.parameters['molecular_proportion'] = molecular_proportion
        else:
            # determine the individual molecular components
            if self.parameters['photosensitizer_selection'] == 'A3B_4Zn':
                # calculate the volume of photosensitizing region in the molecule
                center_porphyrin_length = sci_notation(2*(chemical_dimensions['bond']['c-c']*(2*cos(radians(chemical_dimensions['angle']['sp2']-90))+cos(radians(180-chemical_dimensions['angle']['sp2'])))), sigfigs) * angstrom
                sp2_extension = sci_notation(chemical_dimensions['bond']['c-c'] * (2 + cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-n'] * cos(radians(180-chemical_dimensions['angle']['sp2'])) + chemical_dimensions['bond']['c-n'], sigfigs) * angstrom
                sp3_diazirine = 0 * angstrom
#                 sp3_diazirine = sci_notation(chemical_dimensions['bond']['c-c']*cos(radians(chemical_dimensions['angle']['sp3']-90)) + 2*chemical_dimensions['bond']['c-c']*(1+cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-c']*sin(radians(chemical_dimensions['angle']['sp3']-90)) + chemical_dimensions['bond']['c-f'], sigfigs)
#                 min_thickness = 2*chemical_dimensions['bond']['c-f']*sin(radians(chemical_dimensions['angle']['sp3']))
#                 max_thickness = 2*float(sp3_diazirine)
#                 print('min', min_thickness, 'max', max_thickness)
#                 thickness = sci_notation(average(min_thickness, max_thickness), sigfigs)

                diazirine_length = float(sp2_extension) + float(sp3_diazirine)
                length = float(sci_notation(float(center_porphyrin_length) + 2*diazirine_length, sigfigs))
                atomic_thickness = 1.5 * angstrom   # https://www.nature.com/articles/ncomms1291
                
                # calculate the proportion of well volume that is constituted by the molecular volume, in meters
                self.variables['molecular_volume'] = float(sci_notation(length**2 * atomic_thickness, sigfigs))
                self.variables['volume_proportion'] = (self.variables['photosensitizers'] * self.variables['molecular_volume']) / self.parameters['well_solution_volume']
                
                # calculate the layer distribution for a slice 
                if photosensitizer_moles_per_square_cm:
                    self.parameters['photosensitizers_per_square_cm'] = photosensitizer_moles_per_square_cm * N_A
                    porphyrin_layer_area = self.parameters['photosensitizers_per_square_cm'] * centi**2
                else:
                    solution_depth = solution_depth*centi
                    tilted_height = length*float(sin(radians(45)))
                    layers = solution_depth / tilted_height
                    photosensitizers_per_layer = self.variables['photosensitizers'] / layers
                    solution_incident_area = self.parameters['well_solution_volume'] / solution_depth
                    orthogonal_area = pi*(length)**2
                    parallel_area = (length)*atomic_thickness
                    average_area = average(orthogonal_area, parallel_area)
                    photosensitizers_layer_area = photosensitizers_per_layer*average_area
                    
                # total porphyrin area proportion
                self.variables['area_proportion'] = photosensitizers_layer_area/solution_incident_area
                
                if self.verbose:
                    print(f'The {solution_depth} m deep solution was divided into {layers} layers')
                    print(f'The center porphyrin object is {center_porphyrin_length} meters')
                    print(f'The benzyl extension is {sp2_extension} meters')
                    print(f'The diazirine is {sp3_diazirine} meters')
                    print(f'The molecular length is {length} meters')
                    print(f'The atomic_thickness is {atomic_thickness} meters')
                    print('The molecular volume is {} cubic meters').format(self.variables['molecular_volume'])
                    print('The photosensitizer volume proportion is {}'.format(self.variables['volume_proportion']))
                    print('The photosensitizer area proportion is {}'.format(self.variables['area_proportion']))
            else:
                print('--> ERROR: The {} porphyrin selection is not defined.'.format(self.parameters['porphyrin_selection']))
                
        self.defined_model = True

    def singlet_oxygen_calculations(self, timestep, total_time, initial_time = 0, healing_kinetics = 5):
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
            self.variables['photon_moles_per_timestep'] = photons_per_second * seconds_per_timestep / N_A
            
            print('photons per timestep: ', self.variables['photons_per_timestep'])
            
            # singlet oxygen calculations
            '''so_from_light = photons_per_second * molecules_dissolved_oxygen * excitation_constant'''
            mw_molecular_oxygen = periodic_table.O.MW*2*kilo          #mg / mole
            dissolved_oxygen_concentration = 9              # mg / L, ambient water quality criteria for DO, EPA         # this must be adjusted for the material surface system to only consider oxygen in the water in the vacinity of the surface material photosensitizer. The continuum assumption of the aqueous solution may be implemented such that only a oxygen within fractional volume of the total solution volume.       
            self.variables['dissolved_mo_molar'] = dissolved_oxygen_concentration*milli / mw_molecular_oxygen # * self.parameters['well_solution_volume'] * N_A
            estimated_excited_photosensitizers = photons_per_second * self.parameters['total_time'] * self.variables['volume_proportion']

            # define the kinetic parameters
            self.variables['quantum_yield'] = self.photosensitizer['quantum_yield']['value'] * self.photosensitizer['so_specificity']['value']
            self.variables['healing'] = healing_kinetics
            
            if self.verbose:
                print('molecular oxygen molecules: ', sci_notation(self.variables['dissolved_mo_moles'], sigfigs))
                print('excited photosensitizer molecules: ', sci_notation(estimated_excited_photosensitizers, sigfigs))
        else:
            print('--> ERROR: The singlet oxygen generation cannot be calculated from the light source')
            
    def geometric_oxidation(self):
        if self.bacterium['shape']['value'] == "sphere":
            # define calculation functions
            def shell_volume(radi_1, radi_2, coeff = 1):
                volume = (4*pi/3) * coeff * (radi_1**3 - radi_2**3)
                return volume
            def cap_volume(r, h):
                volume = (pi/3) * h**2 * (3*r - h)
                return volume
            
            # calculate the cellular dimensions
            self.variables['cell_radius'] = pow((bacterium['cell_volume (pL)']['value'] * pico  / liter * 3) / (4 * pi), 1/3)
            membrane_inner_radius = cell_radius - bacterium['membrane_thickness (nm)']['value'] * nano

            # calculate the cellular and oxidation volumes
            membrane_volume = shell_volume(cell_radius, membrane_inner_radius) # M^3
            outer_h = cell_radius * (1 - cos(radians(self.parameters['oxidation_angle'])))
            cap_volume = cap_volume(cell_radius, outer_h)
            self.variables['oxidized_membrane_volume_ratio'] = cap_volume / membrane_volume
#             oxidized_membrane_volume_ratio = self.parameters['oxidation_angle'] / 360

            #calculate the cellular and oxidation areas
            oxidized_cap_area = 2 * pi * cell_radius * outer_h
            cell_area = 4 * pi * cell_radius ** 2
            self.variables['oxidized_area_ratio'] = oxidized_cap_area / cell_area

            # fatty acid oxidation
            singlet_oxygen_interaction_radius = cell_radius + self.parameters['singlet_oxygen_diffusion_distance']
#             membrane_solution_interface_volume = shell_volume(singlet_oxygen_interaction_radius, cell_radius)

            if 'BCFA' in self.membrane_chemicals:
                self.variables['fa_g/L_conc'] = self.membrane_chemicals['BCFA']['concentration']['value'] 
                if 'SCFA' in self.membrane_chemicals:
                    self.variables['fa_g/L_conc'] += self.membrane_chemicals['SCFA']['concentration']['value']
            else:
                fa_g = self.variables['fa_moles'] = self.variables['fa_g/L_conc'] = self.variables['fa_conc'] = 0
                oxidized_volume = membrane_volume * self.variables['oxidized_membrane_volume_ratio']
                for chemical in self.membrane_chemicals:
                    if re.search('FA', chemical):
                        fa_g = oxidized_volume/liter * self.membrane_chemicals[chemical]['density (g/L)']['value'] * self.membrane_chemicals[chemical]['proportion']['value']
                        self.variables['fa_g/L_conc'] += self.membrane_chemicals[chemical]['density (g/L)']['value'] * self.membrane_chemicals[chemical]['proportion']['value']
                        self.variables['fa_moles'] += fa_g / self.membrane_chemicals[chemical]['mw']
                        
            self.variables['k'] = 2.7E2 * self.variables['fa_g/L_conc'] # https://www.jstage.jst.go.jp/article/jos/68/1/68_ess18179/_pdf/-char/ja
                        
            if self.verbose:
                print('oxidized volume proportion: ', self.variables['oxidized_membrane_volume_ratio'])
                print('oxidized area proportion: ',self.variables['oxidized_area_ratio'], '\n\n')

    def kinetic_calculation(self): 
        # define the first equation
        k_so = self.variables['quantum_yield'] * (self.variables['photon_moles_per_timestep'] * self.variables['volume_proportion']) * self.variables['photosensitizer_molar'] # * micro**2*milli/10 # 0.0000777859848166666
        mo = self.variables['dissolved_mo_molar'] # self.variables['dissolved_mo_moles']
        
        # define the second equation
        k = self.variables['k'] # 1/sec
        fa = self.variables['fa_g/L_conc']
        
        # define constants
        healing_kinetics = self.variables['healing']
        biofilm_threshold = self.parameters['stop_biofilm_threshold']
        death_threshold = self.parameters['oxidized_death_threshold']

        # define the SBML model
        self.model = (f'''
          model pdi_oxidation
            # expressions, chemical quantities in lieu of concentrations to accommodate different compartment volumes 
            mo -> so;  {k_so}*mo
            so + fa -> ofa; {k}*so*fa - {healing_kinetics}*ofa    # the aggregated photons / second must be programmatically inserted into the rate expression         

            # define the first expression 
            mo = {mo}
            so = 0
            
            # define the second expression
            ofa = 0
            fa = {fa}
            
            # define constants
            biofilm = {biofilm_threshold};
            vitality = {death_threshold};
            oxidation := ofa / (ofa + fa);
            
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
            if jupyter:
                display(self.result_df)
            else:
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
            self.simulation_path = '_'.join([str(date.today()), self.parameters['photosensitizer_selection'], self.parameters['bacterial_specie'], str(count)])
            while os.path.exists(self.simulation_path):
                count += 1
                self.simulation_path = '_'.join([str(date.today()), self.parameters['photosensitizer_selection'], self.parameters['bacterial_specie'], str(count)])
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
            if jupyter:
                display(parameters_table)
            else:
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
            if jupyter:
                display(variables_table)
            else:
                print(variables_table)
        
        variables_path = os.path.join(self.simulation_path, 'variables.csv')
        variables_table.to_csv(variables_path)
             

    def define(self, bacterial_species, photosensitizer, photosensitizer_conc, light_source, irradiance = None, well_surface_area = 1, exposure = None, simulation_time = None, molecular_proportion = None):
        # parameterize the simulation
        self.define_bacterium(bacterial_species)
        self.define_photosensitizer(photosensitizer_conc, photosensitizer)
        light = self.define_light(light_source, irradiance, well_surface_area, exposure, simulation_time)
        self.define_photosensitizer_volume(molecular_proportion = None)
        
        self.defined_model = True
        
        return light
        
    def time_to_threshold(self, timestep, total_time, kinetc_constant, simulation_path = None, figure_title = 'Oxidation proportion of prokaryotic membrane fatty acids', x_label = 'time', y_label = 'oxidation proportion'):
        # execute the simulation
        if 'watts' in self.parameters:
            self.singlet_oxygen_calculations(timestep, total_time, kinetc_constant)
            self.geometric_oxidation()
            raw_data = self.kinetic_calculation()
            self.export(simulation_path, figure_title, x_label, y_label, )
        else:
            print(f'--> ERROR: The {light} light source is not supported by the code')
            
        return raw_data