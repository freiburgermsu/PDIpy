# import libraries
from scipy.constants import femto, pico, angstrom, nano, micro, milli, centi, kilo, liter, N_A, h, c, minute, hour
from math import pow, pi, sin, cos, radians, ceil
from chemicals.elements import periodic_table
from matplotlib import pyplot
from datetime import date
from sigfig import round
from glob import glob
import tellurium
import graphviz
import pandas
import json, sys, re, os

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
    
def isnumber(num):
    try:
        float(num)
        return True
    except:
        try:
            int(num)
            return True
        except:
            return False

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
        self.parameters['singlet_oxygen_diffusion_distance (m)'] = 80 * nano       # Moan1990 
        self.parameters['oxidation_angle'] = 5           # degrees
        self.parameters['root_path'] = os.path.dirname(__file__)
        
        # identify options
        self.bacteria = [re.search('(?<=bacteria\\\\)(.+)(?=\.json)', x).group() for x in glob(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'*.json'))]
        self.light_parameters = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'light_source.json')))
        
    def define_system(self, timestep, total_time, initial_time = 0, healing_kinetics = 5, bacterial_cfu_ml = 1E6, solution_dimensions = {}, photosensitizer_mg_per_sqr_cm = None, surface_system = False, biofilm = False, individual_cell = False, medium = 'water'):
        # parameterize the photosensitizer system
        self.surface_system = surface_system
        self.biofilm = biofilm
        self.individual_cell = individual_cell
        self.parameters['medium'] = medium
        
        # timestep conditions
        self.bacterial_cfu_ml = bacterial_cfu_ml
        self.parameters['initial_time'] = initial_time
        self.parameters['total_time (s)'] = total_time * minute
        self.parameters['timestep (s)'] = timestep * minute
                
        if self.surface_system:
            if photosensitizer_mg_per_sqr_cm is None:
                self.parameters['photosensitizer_mg_per_disc'] = 0.09
                self.parameters['cm2_per_disc'] = (1.56/2)**2*pi
                self.parameters['surface_area (m^2)'] = self.parameters['cm2_per_disc']*centi**2
                photosensitizer_mg_per_sqr_cm = self.parameters['photosensitizer_mg_per_disc']/self.parameters['cm2_per_disc']
            self.parameters['photosensitizer (g/m^2)'] = photosensitizer_mg_per_sqr_cm*milli/centi**2
            
        else:
            if solution_dimensions != {}:
                solution_volume = solution_dimensions['solution_volume (m^3)']
                surface_area = solution_dimensions['surface_area (m^2)']
                solution_depth = solution_dimensions['solution_depth (m)']
            else:
                print('-> ERROR: The solution_dimensions dictionary must be populated for solution systems.')
                
            if solution_volume is None:
                self.parameters['solution_depth (m)'] = solution_depth
                self.parameters['surface_area (m^2)'] = surface_area
                self.parameters['solution_volume (m^3)'] = self.parameters['solution_depth (m)']*self.parameters['surface_area (m^2)']
            elif surface_area is None:
                self.parameters['solution_volume (m^3)'] = solution_volume
                self.parameters['solution_depth (m)'] = solution_depth
                self.parameters['surface_area (m^2)'] = self.parameters['solution_volume (m^3)']/solution_depth
            elif solution_depth is None:
                self.parameters['solution_volume (m^3)'] = solution_volume
                self.parameters['surface_area (m^2)'] = surface_area
                self.parameters['solution_depth (m)'] = self.parameters['solution_volume (m^3)']/self.parameters['surface_area (m^2)']
            
                                        
    def define_bacterium(self, bacterial_specie):    
        self.results['cellular_vitality'] = True
        self.results['biofilm_growth'] = True
        
        # load the bacterial parameters
        if bacterial_specie in self.bacteria:
            self.parameters['bacterial_specie'] = bacterial_specie
            self.bacterium = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'{bacterial_specie}.json')))[bacterial_specie]
            self.membrane_chemicals = self.bacterium['membrane_chemicals']
        elif type(bacterial_species) is dict:
            pass
        else:
            error = f'--> ERROR: The {bacterial_species} bacterial specie is neither a predefined option nor is a customized dictionary.'
            self.messages.append(error)
            print(error)
        
        
        # fatty acid concentrations in the phospholipid membrane
        self.variables['fa_g/L_conc'] = total_proportion = 0
        for chemical in self.membrane_chemicals:
            if re.search('FA', chemical):
                fa_density_proportion = self.membrane_chemicals[chemical]['density (g/L)']['value'] * self.membrane_chemicals[chemical]['proportion']['value']
                self.variables['fa_molar'] = fa_density_proportion / average(self.membrane_chemicals[chemical]['mw'])
                self.variables['fa_g/L_conc'] += fa_density_proportion
                total_proportion += self.membrane_chemicals[chemical]['proportion']['value']
        

    def define_photosensitizer(self, photosensitizer_molar = None, photosensitizer = 'A3B_4Zn', molecular_proportion = None):
        def cylinder_volume(radius, height):
            return radius**2*pi*height
        
        # load the photosensitizer parameters
        self.parameters['photosensitizer_selection'] = photosensitizer
        self.photosensitizer = json.load(open('{}/parameters/photosensitizers.json'.format(self.parameters['root_path']), encoding='utf-8'))[photosensitizer]
        self.parameters['soret (m)'] = {'upper': self.photosensitizer['soret (nm)']['value'][1]*nano, 'lower': self.photosensitizer['soret (nm)']['value'][0]*nano}
        self.parameters['q (m)'] = {'upper': self.photosensitizer['q (nm)']['value'][1]*nano, 'lower': self.photosensitizer['q (nm)']['value'][0]*nano}
        
        photosensitizer_length = self.photosensitizer['dimensions']['length (A)']*angstrom
        photosensitizer_width = self.photosensitizer['dimensions']['width (A)']*angstrom
        photosensitizer_depth = self.photosensitizer['dimensions']['depth (A)']*angstrom
        photosensitizer_shape = self.photosensitizer['dimensions']['shape']

        if not self.surface_system:
            self.parameters['photosensitizer_molar'] = photosensitizer_molar
            if self.parameters['photosensitizer_molar'] is None:
                print('-> ERROR: A molar photosensitizer concentration must be defined for systems of dissolved photosensitizers.') 
            self.variables['photosensitizers'] = (self.parameters['photosensitizer_molar']*N_A) * (self.parameters['solution_volume (m^3)']/liter)
        else:
            self.variables['photosensitizers'] = self.parameters['photosensitizer (g/m^2)'] / (self.photosensitizer['mw']['value']/N_A) * self.parameters['surface_area (m^2)']
            self.parameters['solution_volume (m^3)'] = self.parameters['surface_area (m^2)']*photosensitizer_length
            self.parameters['photosensitizer_molar'] = self.variables['photosensitizers']/N_A / (self.parameters['solution_volume (m^3)']/liter)            

        # calculate the volume proportion that is constituted by the photosensitizer, in the volume where the photosensitizer is present        
        if photosensitizer_shape == 'cylinder':
            self.variables['molecular_volume (m^3)'] = cylinder_volume(photosensitizer_length/2, photosensitizer_depth)
        molecules_volume = self.variables['photosensitizers'] * self.variables['molecular_volume (m^3)']
        self.variables['volume_proportion'] = molecules_volume / self.parameters['solution_volume (m^3)']
            
        # print calculated content
        self.defined_model = True
        if self.verbose:
            message5 = ''
            message1 = 'The photosensitizer dimensions are {} x {} x {} m, in a {} shape.'.format(photosensitizer_length, photosensitizer_width, photosensitizer_depth, photosensitizer_shape)
            if not self.surface_system:
                message5 = 'The {} m deep solution was divided into {} layers'.format(self.parameters['solution_depth (m)'], ceil(layers))
            message3 = 'The molecular volume is {} cubic meters'.format(sigfigs_conversion(self.variables['molecular_volume (m^3)']))
            message4 = 'The volume proportion of {} photosensitizers is ({} m^2)/({} m^2) = {}'.format(sigfigs_conversion(self.variables['photosensitizers']), sigfigs_conversion(molecules_volume), sigfigs_conversion(self.parameters['solution_volume (m^3)']), sigfigs_conversion(self.variables['volume_proportion']))
            
            messages = [message1, message3, message4, message5]
            self.messages.extend(messages)
            for message in messages:            
                print(message)
                

    def define_light(self, light_source, irradiance = None, exposure = None, simulation_time = None, lux = None, lumens = None,):
        # define properties of the light source
        self.parameters['visible (nm)'] = {'upper': 780 * nano, 'lower': 390 * nano}
        if light_source in self.light_parameters.keys(): 
            lumens_per_watt = self.light_parameters[light_source]['lumens_per_watt']['value']
            self.parameters['visible_proportion'] = self.light_parameters[light_source]['visible_proportion']['value']
        elif type(light_source) is dict:
            lumens_per_watt = light_source['lumens_per_watt']
            self.parameters['visible_proportion'] = light_source['visible_proportion']
        else:
            error = f'--> ERROR: The light source {light_source} is neither a predefined option nor is a customized dictionary.'
            self.messages.append(error)
            print(error)
        
        # define the available light
        if irradiance is not None: # mW / cm^2
            self.parameters['watts'] = (irradiance*milli/centi**2) * self.parameters['surface_area (m^2)']
        elif exposure is not None:  # J / cm^2
            simulation_time = simulation_time * minute
            self.parameters['watts'] = (exposure/centi**2) / simulation_time * self.parameters['surface_area (m^2)']
        elif lux is not None: # lumen / m^2
            self.parameters['watts'] = lux / lumens_per_watt * self.parameters['surface_area (m^2)']
        elif lumens is not None: # lumen
            self.parameters['watts'] = lumens / lumens_per_watt
        else:
            error = '--> ERROR: The light source has an unrecognized dimension.'
            self.messages.append(error)
            print(error)

            
    def singlet_oxygen_calculations(self, excitation_calculation = False):
        if not self.defined_model:
            import sys
            error = 'ERROR: The model must first be defined.'
            self.messages.append(error)
            sys.exit(error)
            
        if 'watts' in self.parameters:                   
            # define the light watts
            effective_visible_light_watts = self.parameters['watts'] * self.parameters['visible_proportion']
            visible_region = self.parameters['visible (nm)']['upper'] - self.parameters['visible (nm)']['lower']
            excitation_visible_proportion = (self.parameters['q (m)']['upper'] - self.parameters['soret (m)']['lower']) / visible_region
            effective_excitation_watts = excitation_visible_proportion * effective_visible_light_watts  # homogeneous light intesity throughout the visible spectrum is assumed
            
            # photonic calculations
            average_excitation_wavelength = (self.parameters['q (m)']['upper'] + self.parameters['soret (m)']['lower']) / 2
            joules_per_photon = (h * c) / average_excitation_wavelength
            non_reflected_photons = 0.96 # “SINGLET OXYGEN GENERATION BY PORPHYRINS AND THE PHOTOSENSITIZATION IN LIPOSOMES KINETICS OF 9,lO-DIMETHYLANTHRACENE” by Gross et al., 1993
            self.variables['photon_moles_per_timestep'] = effective_excitation_watts*non_reflected_photons / joules_per_photon / N_A * self.parameters['timestep (s)']
            
            # singlet oxygen calculations
            '''so_from_light = photons_per_second * molecules_dissolved_oxygen * excitation_constant'''
            mw_molecular_oxygen = periodic_table.O.MW*2*kilo          #mg / mole
            self.variables['dissolved_mo_molar'] = 9*milli/mw_molecular_oxygen    # mg / L, ambient water quality criteria for Dissolved Oxygen, EPA   

            # calculate the proportion of excited photosensitizers
            self.excitation_calculation = excitation_calculation
            if self.excitation_calculation:
                self.variables['e_ps_calculated'] = self.variables['photon_moles_per_timestep']*self.variables['volume_proportion']*self.photosensitizer['quantum_yield']['value']
                print('e_ps times greater than the photosensitizer molar', self.variables['e_ps_calculated']/self.parameters['photosensitizer_molar'])
                if self.variables['e_ps_calculated'] > self.parameters['photosensitizer_molar']:
                    self.variables['e_ps_calculated'] = self.parameters['photosensitizer_molar']
                
            if self.verbose:
                message1 = 'photons per timestep: ', self.variables['photon_moles_per_timestep']
                message2 = 'molecular oxygen molecules: ', sigfigs_conversion(self.variables['dissolved_mo_molar'])
                message4 = 'effective excitation watts: ', sigfigs_conversion(effective_excitation_watts)
                messages = [message1, message2, message4]
                self.messages.extend(messages)
                for message in messages:            
                    print(message)
        else:
            error = '--> ERROR: The singlet oxygen generation cannot be calculated from the light source'
            self.messages.append(error)
            print(error)
        
        # eps oxidation
        if self.biofilm:
            self.variables['eps_oxidation'] = 1e10 # empirical rate constant that represents a biofilm system from the single cellular kinetic model via a competiting oxidation reaction of EPS versus fatty acids   
        
    def kinetic_calculation(self): 
        # ============= define photosensitizer excitation ==============
        ps = self.parameters['photosensitizer_molar']
        if not self.excitation_calculation:
            e_fraction = (self.variables['volume_proportion']*self.variables['photon_moles_per_timestep'])/(self.variables['photosensitizers']/N_A)
            if e_fraction > 1:
                e_fraction = 1

            self.variables['e_ps_decay_time (s)'] = 1.5 * nano # “Ultrafast excitation transfer and relaxation inlinear and crossed-linear arrays of porphyrins” by Akimoto et al., 1999
            self.variables['ps_excitation (s)'] = 50*femto # an estimated time that is below the detection limit of femto-second laser spectrophotometers ; literature has not been discovered that illuminates this time. # self.parameters['photosensitizer_molar']/(self.variables['photon_moles_per_timestep']*self.variables['volume_proportion']) 
            k_e_ps = 1/self.variables['ps_excitation (s)']
            k_ps_rlx = 1/self.variables['e_ps_decay_time (s)']        
            qy_e = self.photosensitizer['quantum_yield']['value']
            
            photosensitizer = f'{k_e_ps}*{e_fraction}*{qy_e}*ps - {k_ps_rlx}*e_ps'
        else:
            photosensitizer = self.variables['e_ps_calculated']
        
        # ============== define photobleaching ==============
        k_b_ps = self.variables['hv_photobleaching'] = self.photosensitizer['photobleaching_constant (cm^2/J)']['value'] * (self.parameters['watts']/(self.parameters['surface_area (m^2)']/centi**2))
        
        # ============== define singlet oxygen generation ==============
        mo = self.variables['dissolved_mo_molar']
        qy = self.variables['so_qy'] = self.photosensitizer['quantum_yield']['value'] * self.photosensitizer['so_specificity']['value'] 
        self.variables['e_ps_charge_transfer (s)'] = 500 * nano   # “Kinetics and efficiency of excitation energy transfer from chlorophylls, their heavy metal-substituted derivatives, and pheophytins to singlet oxygen” by Küpper et al., 2002  & “The role of singlet oxygen and oxygen concentration in photodynamic inactivation of bacteria” by Maisch et al., 2007 
        k_so = 1/self.variables['e_ps_charge_transfer (s)']
        
        cfu_lifetime_slope = (40-10) / (1E8-1E4)               # “The role of singlet oxygen and oxygen concentration in photodynamic inactivation of bacteria” by Maisch et al., 2007
        initial_lifetime = 10
        self.variables['so_decay_time (s)'] = (cfu_lifetime_slope*self.bacterial_cfu_ml+initial_lifetime)*micro
        if self.variables['so_decay_time (s)'] < 4*micro:
            self.variables['so_decay_time (s)'] = 4*micro  # a minimum lifetime of 4 microseconds is parameterized for the default in a complete aqueous solution 
        k_rlx_so = 1/self.variables['so_decay_time (s)']
        
        self.variables['so_rise_time (s)'] = 2 * micro        # “Time-Resolved Investigations of Singlet Oxygen Luminescence in Water, in Phosphatidylcholine, and in Aqueous Suspensions of Phosphatidylcholine or HT29 Cells” by Baier et al., 2005
        
        # ============== define fatty acid and EPS oxidation ==============
        k_fa = self.variables['k_fa'] = 7.69E2 * self.variables['fa_g/L_conc'] # https://www.jstage.jst.go.jp/article/jos/68/1/68_ess18179/_pdf/-char/ja
        fa = self.variables['fa_molar']
        biofilm = ''
        if not self.individual_cell:
            k_fa_reduction = (self.bacterial_cfu_ml/1E6)**0.1
            k_fa /= k_fa_reduction
            if self.biofilm:
                biofilm = 'so => o_eps + mo; {}*so'.format(self.variables['eps_oxidation'])            

        # ============== SBML kinetic model ==============
        self.model = (f'''
          model pdipy_oxidation
            # kinetic expressions
            ps -> e_ps; {photosensitizer}
            ps => b_ps ; {k_b_ps}*ps
            e_ps + mo => so + ps;  {qy}*{k_so}*e_ps*mo
            so => mo; {k_rlx_so}*so
            so + fa => o_fa + mo; {k_fa}*so*fa
            {biofilm}

            # define concentrations
            ps = {ps}
            e_ps = 0
            b_ps = 0
            mo = {mo}
            so = 0
            o_fa = 0
            fa = {fa}
            o_eps = 0
            
            # calculate the oxidation proportion
            oxidation := o_fa / (o_fa + fa);
            
          end
        ''')       
        # ============== model SED-ML plot ==============
        total_points = self.parameters['total_time (s)'] / self.parameters['timestep (s)'] 
        self.phrasedml_str = '''
          model1 = model "pdipy_oxidation"
          sim1 = simulate uniform({}, {}, {})
          task1 = run sim1 on model1
          plot "Oxidation proportion of prokaryotic membrane fatty acids" time vs oxidation
        '''.format(self.parameters['initial_time'], self.parameters['total_time (s)'], total_points)
        
        # ============== execute the model ==============
        tellurium_model = tellurium.loada(self.model)    
        result = tellurium_model.simulate(self.parameters['initial_time'], int(self.parameters['total_time (s)']), int(total_points))
        
        # ============== process the data ==============
        self.result_df = pandas.DataFrame(result)
        self.result_df.index = self.result_df[0]
        del self.result_df[0]
        self.result_df.index.name = 'Time (s)'
        if self.biofilm:
            self.result_df.columns = ['[ps]', '[e_ps]', '[b_ps]', '[mo]', '[so]', '[fa]', '[ofa]', '[oeps]']
        else:
            self.result_df.columns = ['[ps]', '[e_ps]', '[b_ps]', '[mo]', '[so]', '[fa]', '[ofa]']
        
        if self.verbose:
            message1 = tellurium_model.getCurrentAntimony()
            message2 = '\nCurrent integrator:', '\n', tellurium_model.integrator
            message3 = 'k_fa reduction factor:', k_fa_reduction

            messages = [message1, message2, message3]
            self.messages.extend(messages)
            for message in messages:            
                print(message)
            
            if self.jupyter:
                display(self.result_df)
            else:
                print(self.result_df)
            
        return self.result_df
    
        
    def export(self, simulation_path = None, simulation_name = None, figure_title = 'Oxidation proportion of prokaryotic membrane fatty acids', x_label = 'Time (min)', y_label = 'proportion of total', display_excitation_proportion = False):
        # parse the simulation results
        x_values = []
        oxidation_y_values = []
        excitation_y_values = []
        processed_data_dictionary = {}
        for index, point in self.result_df.iterrows():
            x_values.append(index/minute) # sigfigs_conversion(index))
            oxidation_proportion = point['[ofa]'] / (point['[ofa]']+point['[fa]']) 
            oxidation_y_values.append(oxidation_proportion)
            processed_data_dictionary[index/hour] = oxidation_proportion
            
            excitation_proportion = point['[e_ps]'] / (point['[ps]']+point['[e_ps]']+point['[b_ps]']) 
            excitation_y_values.append(excitation_proportion)
        
        # define the simulation_path
        self.simulation_path = simulation_path
        if self.simulation_path is None:
            if simulation_name is None:
                simulation_name = '-'.join([re.sub(' ', '_', str(x)) for x in [date.today(), 'PDIpy', self.parameters['photosensitizer_selection'], self.parameters['bacterial_specie']]])
            count = 0
            while os.path.exists(simulation_name):
                count += 1
                simulation_name = re.sub('([0-9]+)$', str(count), simulation_name)
                if not re.search('(-[0-9]+$)', simulation_name):
                    simulation_name += f'-{count}'
            self.simulation_path = os.path.join(os.getcwd(), simulation_name)        
        os.mkdir(self.simulation_path)
        
        # create and export the OMEX file
        inline_omex = '\n'.join([self.model, self.phrasedml_str])  
        omex_file_name = 'input.omex'
        tellurium.exportInlineOmex(inline_omex, os.path.join(self.simulation_path, omex_file_name))
        self.result_df.to_csv(os.path.join(self.simulation_path, 'raw_data.csv'))
        
        # export the processed data
        self.processed_data = pandas.DataFrame(list(processed_data_dictionary.items()), columns = ['time (hr)','oxidation_proportion'])
        self.processed_data.index = self.processed_data['time (hr)']
        del self.processed_data['time (hr)']
        self.processed_path = os.path.join(self.simulation_path, 'processed_data.csv')
        self.processed_data.to_csv(self.processed_path)
        
        # export the figure
        figure_path = os.path.join(self.simulation_path, 'output.svg')
        pyplot.rcParams['figure.figsize'] = (11, 7)
        pyplot.rcParams['figure.dpi'] = 150
        figure, ax = pyplot.subplots()
        ax.plot(x_values, oxidation_y_values, label = 'oxidation_proportion')
        if display_excitation_proportion:
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
                
        return self.processed_data
    

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
        
        
    def data_processing(self, target_reduction = None, target_time = None):
        if isnumber(target_reduction):
            if target_reduction > self.processed_data.iat[-1,0]:
                message = '--> ERROR: The inquired reduction is never reached.'
            else:
                for index, point in self.processed_data.iterrows():
                    if point['oxidation_proportion'] >= target_reduction:
                        quantity = index
                        unit = 'hour'
                        message = f'time to target ({unit}): {quantity}'
                        break
                    
        elif isnumber(target_time):
            if target_time > self.processed_data.index[-1]:
                message = '--> ERROR: The inquired time is never reached.'
            else:
                for index, point in self.processed_data.iterrows():
                    if index >= target_time:
                        quantity = point['oxidation_proportion']*100
                        unit = '%'
                        message = f'reduction at target ({unit}): {quantity}'
                        break
        else:
            message = '--> ERROR: Neither the target_time nor the target_reduction are parameterized as numbers.'
            quantity = None
            
        print(message)
        self.messages.append(message)
                
        return quantity, unit