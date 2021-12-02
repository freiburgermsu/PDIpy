# import libraries
from scipy.constants import femto, pico, angstrom, nano, micro, milli, centi, kilo, liter, N_A, h, c, minute
from math import pow, pi, sin, cos, radians, ceil
from chemicals.elements import periodic_table
from matplotlib import pyplot
from datetime import date
from sigfig import round
import tellurium
import graphviz
import pandas
import json
import sys
import re
import os

pandas.set_option('max_colwidth', None)

# define useful functions
def sigfigs_conversion(num, sigfigs_in = 2):
    return round(num, sigfigs=sigfigs_in, notation = 'sci')

def average(num_1, num_2 = None):
    if type(num_1) is list:
        return sum(num_1) / len(num_1)
    else:
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


class PDIBacterialPkg():
    def __init__(self, verbose = False, jupyter = False):
        """ Define the initial PDI and biological parameters
        """        
        # define base organizations
        self.parameters = {}
        self.variables = {}
        self.results = {}
        self.verbose = verbose
        self.jupyter = jupyter
        self.messages = []
        
        # initial parameters
        self.parameters['singlet_oxygen_diffusion_distance (nm)'] = 80 * nano       # Moan1990 
        self.parameters['oxidation_angle'] = 5           # degrees
        self.parameters['root_path'] = os.path.dirname(__file__)
        
    def define_system(self, surface_area = None, solution_volume = None, solution_depth = None, photosensitizer_mg_per_sqr_cm = None, surface_system = False, medium = 'water'):
        # parameterize the photosensitizer system
        self.surface_system = surface_system
        self.parameters['medium'] = medium
        if self.surface_system:
            if photosensitizer_mg_per_sqr_cm is None:
                self.parameters['photosensitizer_mg_per_disc'] = 0.3
                self.parameters['cm2_per_disc'] = 1.5**2*pi
                self.parameters['surface_area (m^2)'] = self.parameters['cm2_per_disc']*centi**2
                photosensitizer_mg_per_sqr_cm = self.parameters['photosensitizer_mg_per_disc']/self.parameters['cm2_per_disc']
            self.parameters['photosensitizer_mg_per_sqr_cm'] = photosensitizer_mg_per_sqr_cm
            
        elif surface_area is None:
            self.parameters['solution_volume (m^3)'] = solution_volume*micro*liter
            self.parameters['solution_depth (m)'] = solution_depth*centi
            self.parameters['surface_area (m^2)'] = self.parameters['solution_volume (m^3)'] / (solution_depth*centi)
            
        elif solution_volume is None:
            self.parameters['solution_depth (m)'] = solution_depth*centi
            self.parameters['surface_area (m^2)'] = surface_area*centi**2
            self.parameters['solution_volume (m^3)'] = self.parameters['solution_depth (m)'] * self.parameters['surface_area (m^2)']
            
        elif solution_depth is None:
            self.parameters['solution_volume (m^3)'] = solution_volume*centi**3
            self.parameters['surface_area (m^2)'] = surface_area*centi**2
            self.parameters['solution_depth (m)'] = self.parameters['solution_volume (m^3)'] / self.parameters['surface_area (m^2)']
                                        
    def define_bacterium(self, bacterial_specie):    
        self.results['cellular_vitality'] = True
        self.results['biofilm_growth'] = True
        
        # load the bacterial parameters
        self.parameters['bacterial_specie'] = bacterial_specie
        self.bacterium = json.load(open('{}/parameters/{}.json'.format(self.parameters['root_path'], bacterial_specie)))[bacterial_specie]
        self.membrane_chemicals = self.bacterium['membrane_chemicals']

    def define_photosensitizer(self, photosensitizer_molar = None, photosensitizer = 'A3B_4Zn', molecular_proportion = None, photosensitizer_moles_per_square_cm = None):
        # load the photosensitizer parameters
        self.parameters['photosensitizer_selection'] = photosensitizer
        self.photosensitizer = json.load(open('{}/parameters/photosensitizers.json'.format(self.parameters['root_path'])))[photosensitizer]
        self.parameters['soret (m)'] = {'upper': self.photosensitizer['soret (nm)']['value'][1] * nano, 'lower': self.photosensitizer['soret (nm)']['value'][0] * nano}
        self.parameters['q (m)'] = {'upper': self.photosensitizer['q (nm)']['value'][1] * nano, 'lower': self.photosensitizer['q (nm)']['value'][0] * nano}
        
        self.parameters['atomic_thickness'] = 1.5 * angstrom   # https://www.nature.com/articles/ncomms1291
        if self.parameters['photosensitizer_selection'] == 'A3B_4Zn':
            # calculate the volume of conjugated region of the photosensitizer
            self.variables['center_porphyrin_length'] = (2*(chemical_dimensions['bond']['c-c']*(2*cos(radians(chemical_dimensions['angle']['sp2']-90))+cos(radians(180-chemical_dimensions['angle']['sp2']))))) * angstrom
            self.variables['sp2_extension'] = (chemical_dimensions['bond']['c-c'] * (2 + cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-n'] * cos(radians(180-chemical_dimensions['angle']['sp2'])) + chemical_dimensions['bond']['c-n']) * angstrom
            self.variables['sp3_diazirine'] = (chemical_dimensions['bond']['c-c']*cos(radians(chemical_dimensions['angle']['sp3']-90)) + 2*chemical_dimensions['bond']['c-c']*(1+cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-c']*sin(radians(chemical_dimensions['angle']['sp3']-90)) + chemical_dimensions['bond']['c-f']) * angstrom
            self.variables['linked_sp3_diazirine'] = (chemical_dimensions['bond']['c-c']*cos(radians(chemical_dimensions['angle']['sp3']-90)) + 2*chemical_dimensions['bond']['c-c']*(1+cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-c']*sin(radians(chemical_dimensions['angle']['sp3']-90))) * angstrom
            
            total_length = self.variables['center_porphyrin_length'] + 2*(self.variables['sp2_extension'] + self.variables['sp3_diazirine'])
            linked_total_length = self.variables['center_porphyrin_length'] + 2*self.variables['sp2_extension'] + self.variables['sp3_diazirine'] + self.variables['linked_sp3_diazirine']

        if not self.surface_system:
            self.parameters['photosensitizer_molar'] = photosensitizer_molar
            self.variables['photosensitizers'] = (self.parameters['photosensitizer_molar']*N_A) * (self.parameters['solution_volume (m^3)']/liter)
        else:
            self.variables['photosensitizers'] = self.parameters['photosensitizer_mg_per_sqr_cm']*milli/centi**2 / self.photosensitizer['mw']['value']*N_A * self.parameters['surface_area (m^2)']
            self.parameters['solution_volume (m^3)'] = self.parameters['surface_area (m^2)']*linked_total_length
            self.parameters['photosensitizer_molar'] = self.variables['photosensitizers']/N_A / (self.parameters['solution_volume (m^3)']/liter)
            
        # determine the individual molecular components
        if self.parameters['photosensitizer_selection'] == 'A3B_4Zn':
            conjugated_length = self.variables['center_porphyrin_length'] + 2*self.variables['sp2_extension']

        # calculate the proportion of well volume that is constituted by the molecular volume, in meters
        self.variables['molecular_volume (m^3)'] = conjugated_length**2 * self.parameters['atomic_thickness']
        molecules_volume = (self.variables['photosensitizers'] * self.variables['molecular_volume (m^3)'])
        self.variables['volume_proportion'] = molecules_volume / self.parameters['solution_volume (m^3)']

        # calculate the layer distribution for a slice 
        if photosensitizer_moles_per_square_cm:
            photosensitizers_layer_area = self.parameters['photosensitizers_per_square_cm'] * centi**2
            self.parameters['photosensitizers'] = photosensitizer_moles_per_square_cm * N_A * photosensitizers_layer_area
        else:
            tilted_height = conjugated_length*float(sin(radians(45)))
            layers = 1
            if not self.surface_system:
                layers = self.parameters['solution_depth (m)'] / tilted_height
            photosensitizers_per_layer = self.variables['photosensitizers'] / layers
            orthogonal_area = pi*(conjugated_length)**2
            parallel_area = (conjugated_length)*self.parameters['atomic_thickness']
            photosensitizers_layer_area = photosensitizers_per_layer*average(orthogonal_area, parallel_area)
            
        # calculate the volume proportion in the cross-linked layer
        self.parameters['layer_volume (m^3)'] = self.parameters['solution_volume (m^3)'] / layers
        self.variables['volume_proportion'] = molecules_volume / self.parameters['layer_volume (m^3)']

        # total porphyrin area proportion
        self.variables['area_proportion'] = photosensitizers_layer_area/self.parameters['surface_area (m^2)']               
        self.defined_model = True
        if self.verbose:
            message1 = 'The center porphyrin object is {} meters'.format(sigfigs_conversion(self.variables['center_porphyrin_length']))
            message2 = 'The benzyl extension is {} meters'.format(sigfigs_conversion(self.variables['sp2_extension']))
            message3 = 'The diazirine is {} meters'.format(sigfigs_conversion(self.variables['sp3_diazirine']))
            message4 = 'Solution depth is not calculated for surface systems'
            if not self.surface_system:
                message4 = 'The {} m deep solution was divided into {} layers'.format(self.parameters['solution_depth (m)'], ceil(layers))
            message5 = f'The molecular length is {sigfigs_conversion(conjugated_length)} meters'
            message6 = 'The molecular volume is {} cubic meters'.format(sigfigs_conversion(self.variables['molecular_volume (m^3)']))
            message7 = 'The photosensitizer volume proportion is {}'.format(sigfigs_conversion(self.variables['volume_proportion']))
            message8 = 'The photosensitizer area proportion is {}'.format(sigfigs_conversion(self.variables['area_proportion']))
            self.messages.extend([message1, message2, message3, message4, message5, message6, message7, message8])
            
            print(message1)
            print(message2)
            print(message3)
            print(message4)
            print(message5)            
            print(message6)
            print(message7)
            print(message8)

    def define_light(self, light_source, irradiance = None, exposure = None, simulation_time = None, lux = None, lumens = None, wattage = None, distance = None, reflection = False):
        self.parameters['visible (nm)'] = {'upper': 780 * nano, 'lower': 390 * nano}

        light_parameters = json.load(open(os.path.join(self.parameters['root_path'],'parameters','light_source.json')))
        self.parameters['visible_proportion'] = light_parameters[light_source]['visible_proportion']['value']
        
        # define the available light
        if irradiance is not None: # mW / cm^2
            self.parameters['watts'] = (irradiance*milli/centi**2) * self.parameters['surface_area (m^2)']
        elif exposure is not None:  # J / cm^2
            simulation_time = simulation_time * minute
            self.parameters['watts'] = (exposure/centi**2) / simulation_time * self.parameters['surface_area (m^2)']
        elif lux is not None: # lumen / cm^2
            irradiance =  (lux/centi**2) / light_parameters[light_source]['lumens_per_watt']['value']
            self.parameters['watts'] = irradiance * self.parameters['surface_area (m^2)']
        elif lumens is not None: # lumen
            self.parameters['watts'] = lumens / light_parameters[light_source]['lumens_per_watt']['value']
        else:
            error = '--> ERROR: The light source has an unrecognized dimension.'
            self.messages.append(error)
            print(error)

    def singlet_oxygen_calculations(self, timestep, total_time, initial_time = 0, healing_kinetics = 5, bacterial_cfu_ml = 1E8):
        if not self.defined_model:
            import sys
            error = 'ERROR: The model must first be defined.'
            self.messages.append(error)
            sys.exit(error)
            
        # timestep conditions
        self.parameters['initial_time'] = initial_time
        self.parameters['total_time (s)'] = total_time * minute
        self.parameters['timestep (s)'] = timestep * minute
        
        # SO lifetime determination
        cfu_lifetime_slope = (40-10) / (1E8-1E4)               # “The role of singlet oxygen and oxygen concentration in photodynamic inactivation of bacteria” by Maisch et al., 2007
        initial_lifetime = 10
        self.variables['so_decay_time (s)'] = (cfu_lifetime_slope*bacterial_cfu_ml+initial_lifetime) * micro 
        self.variables['so_rise_time (s)'] = 2 * micro        # “Time-Resolved Investigations of Singlet Oxygen Luminescence in Water, in Phosphatidylcholine, and in Aqueous Suspensions of Phosphatidylcholine or HT29 Cells” by Baier et al., 2005
        
        # excited photosensitizer lifetime
        self.variables['e_ps_charge_transfer'] = 500 * nano   # “Kinetics and efficiency of excitation energy transfer from chlorophylls, their heavy metal-substituted derivatives, and pheophytins to singlet oxygen” by Küpper et al., 2002  & “The role of singlet oxygen and oxygen concentration in photodynamic inactivation of bacteria” by Maisch et al., 2007     
        
        # photobleaching
        self.variables['hv_photobleaching'] = 0.015 * (self.parameters['watts']/self.parameters['surface_area (m^2)'])  # “Photobleaching kinetics, photoproduct formation, and dose estimation during ALA induced PpIX PDT of MLL cells under well oxygenated and hypoxic conditions” by Dysart et al., 2005
            
        if 'watts' in self.parameters:                   
            # define the light watts
            effective_visible_light_watts = self.parameters['watts'] * self.parameters['visible_proportion']
            visible_region = self.parameters['visible (nm)']['upper'] - self.parameters['visible (nm)']['lower']
            excitation_visible_proportion = ((self.parameters['soret (m)']['upper'] - self.parameters['soret (m)']['lower'])) / visible_region
            effective_excitation_watts = excitation_visible_proportion * effective_visible_light_watts  # homogeneous light intesity throughout the visible spectrum is assumed
            
            # photonic calculations
            average_excitation_wavelength = (self.parameters['q (m)']['upper'] + self.parameters['soret (m)']['lower']) / 2
            joules_per_photon = (h * c) / average_excitation_wavelength
            photons_per_second = effective_excitation_watts / joules_per_photon
            self.variables['photon_moles_per_timestep'] = photons_per_second * self.parameters['timestep (s)'] / N_A
            
            # singlet oxygen calculations
            '''so_from_light = photons_per_second * molecules_dissolved_oxygen * excitation_constant'''
            mw_molecular_oxygen = periodic_table.O.MW*2*kilo          #mg / mole
            dissolved_oxygen_concentration = 9              # mg / L, ambient water quality criteria for DO, EPA         # this must be adjusted for the material surface system to only consider oxygen in the water in the vacinity of the surface material photosensitizer. The continuum assumption of the aqueous solution may be implemented such that only a oxygen within fractional volume of the total solution volume.    
#             if self.parameters['medium'] == 'air':
#                 self.variables['so_relaxation'] = 1e4
#             if self.parameters['medium'] == 'water':
#                 self.variables['so_relaxation'] = 1e10
#             self.variables['ps_relaxation'] = 1e4
            self.variables['dissolved_mo_molar'] = dissolved_oxygen_concentration*milli / mw_molecular_oxygen # * self.parameters['well_solution_volume'] * N_A
            estimated_excited_photosensitizers = photons_per_second * self.parameters['total_time (s)'] * self.variables['volume_proportion']

            # define the kinetic parameters
            self.variables['quantum_yield'] = self.photosensitizer['quantum_yield']['value'] * self.photosensitizer['so_specificity']['value']
            self.variables['healing'] = healing_kinetics
            
            if self.verbose:
                message1 = 'photons per timestep: ', self.variables['photon_moles_per_timestep']
                message2 = 'molecular oxygen molecules: ', sigfigs_conversion(self.variables['dissolved_mo_molar'])
                message3 = 'excited photosensitizer molecules: ', sigfigs_conversion(estimated_excited_photosensitizers)
                message4 = 'effective excitation watts: ', sigfigs_conversion(effective_excitation_watts)
                self.messages.extend([message1, message2, message3, message4])

                print(message1)
                print(message2)
                print(message3)
                print(message4)
        else:
            error = '--> ERROR: The singlet oxygen generation cannot be calculated from the light source'
            self.messages.append(error)
            print(error)
            
    def geometric_oxidation(self):
        if self.bacterium['shape']['value'] == "sphere":
            # define calculation functions
            def shell_volume(radi_1, radi_2, coeff = 1):
                volume = (4*pi/3) * coeff * (radi_1**3 - radi_2**3)
                return volume        
            def sector_volume(r, h):
                volume = (2*pi/3)*r**2*h
                return volume       
            
            # calculate the cellular dimensions
            self.variables['cell_radius (m)'] = pow((self.bacterium['cell_volume (pL)']['value']*pico/liter*3) / (4*pi), 1/3)
            membrane_inner_radius = self.variables['cell_radius (m)'] - self.bacterium['membrane_thickness (nm)']['value'] * nano

            # calculate the cellular and oxidation volumes
            membrane_volume = shell_volume(self.variables['cell_radius (m)'], membrane_inner_radius) # M^3
            outer_h = self.variables['cell_radius (m)'] * (1 - cos(radians(self.parameters['oxidation_angle'])))
            inner_h = membrane_inner_radius * (1 - cos(radians(self.parameters['oxidation_angle'])))
            shell_sector_volume = sector_volume(self.variables['cell_radius (m)'], outer_h) - sector_volume(membrane_inner_radius, inner_h)
            self.variables['oxidized_membrane_volume_ratio'] = shell_sector_volume / membrane_volume

            #calculate the cellular and oxidation areas
            oxidized_cap_area = 2 * pi * self.variables['cell_radius (m)'] * outer_h
            cell_area = 4 * pi * self.variables['cell_radius (m)'] ** 2
            self.variables['oxidized_area_ratio'] = oxidized_cap_area / cell_area

            # fatty acid concentrations in the phospholipid membrane
            self.variables['fa_g/L_conc'] = total_proportion = 0
            oxidized_volume = membrane_volume * self.variables['oxidized_membrane_volume_ratio']
            for chemical in self.membrane_chemicals:
                if re.search('FA', chemical):
                    self.variables['fa_molar'] = self.membrane_chemicals[chemical]['density (g/L)']['value'] / average(self.membrane_chemicals[chemical]['mw']) * self.membrane_chemicals[chemical]['proportion']['value'] # * oxidized_volume/liter 
                    self.variables['fa_g/L_conc'] += self.membrane_chemicals[chemical]['density (g/L)']['value'] * self.membrane_chemicals[chemical]['proportion']['value'] / milli
                    total_proportion += self.membrane_chemicals[chemical]['proportion']['value']
                        
            self.variables['fa_g/L_conc'] /= total_proportion
            self.variables['fa_molar'] /= total_proportion
            
        self.variables['k_fa'] = 2.7E2 * self.variables['fa_g/L_conc'] # https://www.jstage.jst.go.jp/article/jos/68/1/68_ess18179/_pdf/-char/ja
        if self.verbose:
            message1 = 'oxidized volume proportion: ', self.variables['oxidized_membrane_volume_ratio']
            message2 = 'volume:area consistency', round(self.variables['oxidized_area_ratio'],7) == round(self.variables['oxidized_membrane_volume_ratio'],7)
            self.messages.extend([message1, message2])

            print(message1)
            print(message2)

    def kinetic_calculation(self): 
        # define photosensitizer excitation
        ps = self.parameters['photosensitizer_molar']
        hv = self.variables['photon_moles_per_timestep']
        qy = self.variables['quantum_yield']
        k_ps = self.variables['k_ps'] = 1e10
        
        # define photosensitizer photobleaching
        k_hv_ps = self.variables['hv_photobleaching']
        
#         self.parameters['so_rise_time (s)']
        
        # define singlet oxygen generation
        mo = self.variables['dissolved_mo_molar'] # self.variables['dissolved_mo_moles']
        k_so = 1/self.variables['e_ps_charge_transfer']
        k_rlx_so = 1/self.variables['so_decay_time (s)']
        
        # define the third equation
        k_fa = self.variables['k_fa'] # 1/sec
        fa = self.variables['fa_molar']
        
        # define constants
#         so_relaxation = self.variables['so_relaxation']
#         ps_relaxation = self.variables['ps_relaxation'] # - ({ps_relaxation}*e_ps)

        # define the SBML model
        self.model = (f'''
          model pdipy_oxidation
            # kinetic expressions
            ps -> e_ps; {k_ps}*{hv}*ps 
            ps => b_ps ; {k_hv_ps}*ps
            e_ps + mo -> so + ps;  {qy}*{k_so}*e_ps*mo - ({k_rlx_so}*so)
            so + fa => ofa; {k_fa}*so*fa

            # define concentrations
            ps = {ps}
            e_ps = 0
            b_ps = 0
            mo = {mo}
            so = 0
            ofa = 0
            fa = {fa}
            
            # calculate the oxidation proportion
            oxidation := ofa / (ofa + fa);
            
          end
        ''')       
        # define the SEDML plot
        total_points = self.parameters['total_time (s)'] / self.parameters['timestep (s)'] 
        self.phrasedml_str = '''
          model1 = model "pdipy_oxidation"
          sim1 = simulate uniform({}, {}, {})
          task1 = run sim1 on model1
          plot "Oxidation proportion of prokaryotic membrane fatty acids" time vs oxidation
        '''.format(self.parameters['initial_time'], self.parameters['total_time (s)'], total_points)
        
        # execute the model
        tellurium_model = tellurium.loada(self.model)            
#         tellurium_model.draw()
#         print(tellurium_model.draw())
        result = tellurium_model.simulate(self.parameters['initial_time'], int(self.parameters['total_time (s)']), int(total_points))
        
        # process the data
        self.result_df = pandas.DataFrame(result)
        self.result_df.index = self.result_df[0]
        del self.result_df[0]
        self.result_df.index.name = 'Time (s)'
        self.result_df.columns = ['[ps]', '[e_ps]', '[b_ps]', '[mo]', '[so]', '[fa]', '[ofa]']
        
        if self.verbose:
            message1 = tellurium_model.getCurrentAntimony()
            message2 = '\nCurrent integrator:', '\n', tellurium_model.integrator
            self.messages.extend([message1, message2])
            print(message1)
            print(message2)
            
            if self.jupyter:
                display(self.result_df)
            else:
                print(self.result_df)
            
        return self.result_df
        
    def export(self, simulation_path = None, figure_title = 'Oxidation proportion of prokaryotic membrane fatty acids', x_label = 'Time (min)', y_label = 'proportion of total'):
        # parse the simulation results
        x_values = []
        oxidation_y_values = []
        excitation_y_values = []
        for index, point in self.result_df.iterrows():
            x_values.append(index/minute) # sigfigs_conversion(index))
            oxidation_proportion = point['[ofa]'] / (point['[ofa]']+point['[fa]']) 
            oxidation_y_values.append(oxidation_proportion)
            excitation_proportion = point['[ps]'] / (point['[ps]']+point['[e_ps]']) 
            excitation_y_values.append(excitation_proportion)
        
        # define the simulation_path
        self.simulation_path = simulation_path
        if self.simulation_path is None:
            count = 0
            self.simulation_path = '-'.join([re.sub(' ', '_', str(x)) for x in [date.today(), 'PDIpy', self.parameters['photosensitizer_selection'], self.parameters['bacterial_specie']]])
            while os.path.exists(self.simulation_path):
                count += 1
                self.simulation_path = '-'.join([re.sub(' ', '_', str(x)) for x in [date.today(), 'PDIpy', self.parameters['photosensitizer_selection'], self.parameters['bacterial_specie'], str(count)]])
            self.simulation_path = os.path.join(os.getcwd(), self.simulation_path)        
        os.mkdir(self.simulation_path)
        
        # create and export the OMEX file
        inline_omex = '\n'.join([self.model, self.phrasedml_str])  
        omex_file_name = 'input.omex'
        tellurium.exportInlineOmex(inline_omex, os.path.join(self.simulation_path, omex_file_name))
        self.result_df.to_csv(os.path.join(self.simulation_path, 'raw_data.csv'))
        
        # export the figure
        figure_path = os.path.join(self.simulation_path, 'output.svg')
        pyplot.rcParams['figure.figsize'] = (11, 7)
        pyplot.rcParams['figure.dpi'] = 150
        figure, ax = pyplot.subplots()
        ax.plot(x_values, oxidation_y_values, label = 'oxidation_proportion')
        ax.plot(x_values, excitation_y_values, label = 'excitation_proportion')
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        ax.set_title(figure_title)
        ax.legend()
        figure.savefig(figure_path)
        if self.verbose:
            if self.jupyter:
                display(figure)
            else:
                figure.show()
        
        # define a table of parameters
        parameters = {'parameter':[], 'value':[]}
        parameters['parameter'].append('simulation_path')
        parameters['value'].append(self.simulation_path)
        for parameter in self.parameters:
            parameters['parameter'].append(parameter)
            parameters['value'].append(self.parameters[parameter])
        parameters_table = pandas.DataFrame(parameters)
        parameters_path = os.path.join(self.simulation_path, 'parameters.csv')
        parameters_table.to_csv(parameters_path)
        if self.verbose:
            if self.jupyter:
                display(parameters_table)
            else:
                print(parameters_table)
        
        # define a table of variables
        variables = {'variable':[], 'value':[]}
        variables['variable'].append('simulation_path')
        variables['value'].append(self.simulation_path)
        for variable in self.variables:
            variables['variable'].append(variable)
            variables['value'].append(self.variables[variable])
        variables_table = pandas.DataFrame(variables)        
        variables_path = os.path.join(self.simulation_path, 'variables.csv')
        variables_table.to_csv(variables_path)
        if self.verbose:
            if self.jupyter:
                display(variables_table)
            else:
                print(variables_table)             

    def define(self, bacterial_species, photosensitizer, photosensitizer_conc, light_source, irradiance = None, well_surface_area = 1, exposure = None, simulation_time = None, molecular_proportion = None):
        # parameterize the simulation
        self.define_bacterium(bacterial_species)
        self.define_photosensitizer(photosensitizer_conc, photosensitizer)
        self.define_light(light_source, irradiance, well_surface_area, exposure, simulation_time)
        self.photosensitizer_proportion(molecular_proportion = None)
        
    def execute(self, timestep, total_time, kinetc_constant, simulation_path = None, figure_title = 'Oxidation proportion of prokaryotic membrane fatty acids', x_label = 'time', y_label = 'oxidation proportion'):
        self.singlet_oxygen_calculations(timestep, total_time, kinetc_constant)
        self.geometric_oxidation()
        raw_data = self.kinetic_calculation()
        self.export(simulation_path, figure_title, x_label, y_label, )
            
        return raw_data