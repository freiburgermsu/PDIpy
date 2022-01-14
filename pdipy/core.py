# import libraries
from scipy.constants import femto, pico, angstrom, nano, micro, milli, centi, kilo, liter, N_A, h, c, minute, hour
from math import pow, pi, sin, cos, radians, ceil, log
from chemicals.elements import periodic_table
from matplotlib import pyplot
from datetime import date
from pprint import pprint
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
        self.defined_model = {}
        self.verbose = verbose
        self.jupyter = jupyter
        self.messages = []
        
        # initial parameters
        self.parameters['singlet_oxygen_diffusion_distance (m)'] = 80 * nano       # Moan1990 
#         self.parameters['oxidation_angle'] = 5           # degrees
        self.parameters['root_path'] = os.path.dirname(__file__)
        
        # identify options
        self.bacteria = [re.search('(?<=bacteria\\\\)(.+)(?=\.json)', x).group() for x in glob(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'*.json'))]
        self.light_parameters = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'light_source.json')))
        self.photosensitizers = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'photosensitizers.json'), encoding = 'utf-8'))
        
    def define_system(self, timestep, total_time, initial_time = 0, bacterial_cfu_ml = 1E6, solution_dimensions = {}, photosensitizer_mg_per_disc = 0.09, cm2_per_disc = 1.91134, surface_system = False, biofilm = False, well_count = '24', individual_cell = False, medium = 'water'):
        # parameterize the photosensitizer system
        self.surface_system = surface_system
        self.biofilm = biofilm
        self.individual_cell = individual_cell
        self.parameters['medium'] = medium
        
        # timestep conditions
        self.bacterial_cfu_ml = bacterial_cfu_ml
        self.parameters['initial_time'] = initial_time
        self.parameters['total_time (s)'] = int(total_time * minute)
        self.parameters['timestep (s)'] = timestep * minute
        
        # define the solution 
        self.solution = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'wells.json')))[well_count]
        if solution_dimensions != {}:
            for key, value in solution_dimensions.items():
                self.solution[key] = value
            
        operational_depth = 0.7
        self.parameters['solution_depth (m)'] = self.solution['depth (cm)']*operational_depth*centi
        self.parameters['solution_area (m^2)'] = self.solution['area (cm^2)']*centi**2
        self.parameters['solution_volume (m^3)'] = self.parameters['solution_depth (m)']*self.parameters['solution_area (m^2)']
        
        # define the photosensitizing surface
        self.area = self.parameters['solution_area (m^2)']
        if self.surface_system:            
            self.area = self.parameters['surface_area (m^2)'] = cm2_per_disc*centi**2
            self.parameters['photosensitizer (g/m^2)'] = photosensitizer_mg_per_disc*milli/(cm2_per_disc*centi**2)
            
        # confirm the conpletion of the function
        self.defined_model.update({'define_system':True})
        
                
    def define_bacterium(self, bacterial_specie):    
        self.results['cellular_vitality'] = True
        self.results['biofilm_growth'] = True
        
        # load the bacterial parameters
        self.bacterium = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', 'S_aureus.json')))
        if bacterial_specie in self.bacteria:
            self.parameters['bacterial_specie'] = bacterial_specie
            self.bacterium = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'{bacterial_specie}.json')))
        elif type(bacterial_specie) is dict:
            for key, value in bacterial_specie.items():
                if key == 'name':
                    self.parameters['bacterial_specie'] = value
                else:
                    for key2, value2 in value.items():
                        self.bacterial_specie[key][key2] = value2 
        else:
            error = f'--> ERROR: The {bacterial_specie} bacterial specie is neither a predefined option nor is a customized dictionary.'
            self.messages.append(error)
            print(error)
        
        # define qualities of the bacterium
        self.membrane_chemicals = self.bacterium['membrane_chemicals']
                                          
        # fatty acid concentrations in the cytoplasmic membrane
        self.variables['fa_g/L_conc'] = total_proportion = 0
        fa_mws = []
        for chemical in self.membrane_chemicals:
            if re.search('FA', chemical):
                total_proportion += self.membrane_chemicals[chemical]['proportion']['value']
        for chemical in self.membrane_chemicals:
            if re.search('FA', chemical):
                fa_density_proportion = self.membrane_chemicals[chemical]['density (g/L)']['value'] * self.membrane_chemicals[chemical]['proportion']['value']
                self.variables['fa_g/L_conc'] += fa_density_proportion
                mw = average(self.membrane_chemicals[chemical]['mw'])
                fa_mws.append(mw*fa_density_proportion/total_proportion)
        self.variables['fa_molar'] = self.variables['fa_g/L_conc'] / average(fa_mws)
        
        # confirm the conpletion of the function
        self.defined_model.update({'define_bacterium':True})
                

    def define_photosensitizer(self, photosensitizer = 'A3B_4Zn', photosensitizer_molar = None, molecular_proportion = None):
        def cylinder_volume(radius, height):
            return (radius**2)*pi*height
        
        # load the photosensitizer parameters
        self.photosensitizer = self.photosensitizers['A3B_4Zn']
        if type(photosensitizer) is dict:
            for key, value in photosensitizer.items():
                if key == 'name':
                    self.parameters['photosensitizer_selection'] = value
                else:
                    for key2, value2 in value.items():
                        self.photosensitizer[key][key2] = value2
        elif photosensitizer in self.photosensitizers:
            self.parameters['photosensitizer_selection'] = photosensitizer
            self.photosensitizer = self.photosensitizers[photosensitizer]                    
        else: #if photosensitizer not in self.photosensitizers:
            error = f'--> ERROR: The photosensitizer {photosensitizer} is neither a predefined option nor a customized dictionary.'
            self.messages.append(error)
            print(error)      
                                                       
        # absorption characteristics
        upper_soret = self.photosensitizer['soret (nm)']['value'][1]*nano
        lower_soret = self.photosensitizer['soret (nm)']['value'][0]*nano
        upper_q = self.photosensitizer['q (nm)']['value'][1]*nano
        lower_q = self.photosensitizer['q (nm)']['value'][0]*nano
        self.parameters['soret (m)'] = {'upper': upper_soret, 'lower': lower_soret}
        self.parameters['q (m)'] = {'upper': upper_q, 'lower': lower_q}

        # physical dimensions
        photosensitizer_length = self.photosensitizer['dimensions']['length (A)']*angstrom
        photosensitizer_width = self.photosensitizer['dimensions']['width (A)']*angstrom
        photosensitizer_depth = self.photosensitizer['dimensions']['depth (A)']*angstrom
        photosensitizer_shape = self.photosensitizer['dimensions']['shape']
        photosensitizer_mw = self.photosensitizer['mw']['value']
                                          
        if not self.surface_system:
            self.variables['photosensitizer_molar'] = photosensitizer_molar
            if self.variables['photosensitizer_molar'] is None:
                error = '--> ERROR: A molar photosensitizer concentration must be defined for systems of dissolved photosensitizers.'
                self.messages.append(error)
                print(error)
            self.variables['photosensitizers'] = (self.variables['photosensitizer_molar']*N_A) * (self.parameters['solution_volume (m^3)']/liter)
        else:
            self.variables['photosensitizers'] = self.parameters['photosensitizer (g/m^2)'] / (photosensitizer_mw/N_A) * self.parameters['surface_area (m^2)']
            self.parameters['solution_volume (m^3)'] = self.parameters['surface_area (m^2)']*photosensitizer_length
            self.variables['photosensitizer_molar'] = self.variables['photosensitizers']/N_A / (self.parameters['solution_volume (m^3)']/liter)  

        # calculate the volume proportion that is constituted by the photosensitizer, in the volume where the photosensitizer is present        
        if photosensitizer_shape == 'disc':
            self.variables['molecular_volume (m^3)'] = cylinder_volume(photosensitizer_length/2, photosensitizer_depth)
        else:
            error = f'--> ERROR: The volume formula for the {photosensitizer_shape} shape is not available.'
            self.messages.append(error)
            print(error)
        molecules_volume = self.variables['photosensitizers'] * self.variables['molecular_volume (m^3)']
        self.variables['volume_proportion'] = molecules_volume / self.parameters['solution_volume (m^3)']
            
        # print calculated content
        if self.verbose:
            message2 = ''
            message1 = 'The photosensitizer dimensions as a {} = {} m x {} m x {} m.'.format(photosensitizer_shape, photosensitizer_length, photosensitizer_width, photosensitizer_depth)
#             if not self.surface_system:
#                 message2 = 'The {} m deep solution was divided into {} layers'.format(self.parameters['solution_depth (m)'], ceil(layers))
            message3 = 'Photosensitizer volume = {} m^3'.format(sigfigs_conversion(self.variables['molecular_volume (m^3)']))
            message4 = 'The volume proportion of {} photosensitizers = ({} m^3 of photosensitizer)/({} m^3 of solution) = {}'.format(sigfigs_conversion(self.variables['photosensitizers']), sigfigs_conversion(molecules_volume), sigfigs_conversion(self.parameters['solution_volume (m^3)']), sigfigs_conversion(self.variables['volume_proportion']))
            
            messages = [message1, message3, message4, message2]
            self.messages.extend(messages)
            for message in messages:            
                print(message)
                
        # confirm the conpletion of the function
        self.defined_model.update({'define_photosensitizer':True})
                

    def define_light(self, light_source, irradiance = None, exposure = None, simulation_time = None, lux = None, lumens = None,):
        # define properties of the light source
        self.parameters['visible (nm)'] = {'upper': 780 * nano, 'lower': 390 * nano}
        self.light = self.light_parameters['LED']
        if type(light_source) is dict:
            for key, value in light_source.items():
                if key == 'name':
                    self.parameters['light_source'] = value                        
                else:
                    for key2, value2 in value.items():
                        self.light[key][key2] = value2
        elif light_source in self.light_parameters: 
            self.parameters['light_source'] = light_source
            self.light = self.light_parameters[light_source]
        else:
            error = f'--> ERROR: The light source {light_source} is neither a predefined option nor a customized dictionary.'
            self.messages.append(error)
            print(error)
        
        # define the available light
        lumens_per_watt = self.light['lumens_per_watt']['value']
        if irradiance is not None: # mW / cm^2
            self.parameters['watts'] = (irradiance*milli/centi**2) * self.area
        elif exposure is not None:  # J / cm^2
            simulation_time = simulation_time * minute
            self.parameters['watts'] = (exposure/centi**2) / simulation_time * self.area
        elif lux is not None: # lumen / m^2
            self.parameters['watts'] = lux / lumens_per_watt * self.area
        elif lumens is not None: # lumen
            self.parameters['watts'] = lumens / lumens_per_watt
        else:
            error = '--> ERROR: The light source has an unrecognized dimension.'
            self.messages.append(error)
            print(error)
            
        # confirm the conpletion of the function
        self.defined_model.update({'define_light':True})

            
    def singlet_oxygen_calculations(self,):
        for function in self.defined_model:
            if not self.defined_model[function]:
                import sys
                error = f'--> ERROR: The {function} function must be defined before the simulation can execute.'
                self.messages.append(error)
                sys.exit(error)
                
        # define the light watts
        visible_region = self.parameters['visible (nm)']['upper'] - self.parameters['visible (nm)']['lower']
        excitation_visible_proportion = (self.parameters['q (m)']['upper'] - self.parameters['soret (m)']['lower']) / visible_region
        visible_light_watts = self.parameters['watts'] * self.light['visible_proportion']['value']
        effective_excitation_watts = excitation_visible_proportion * visible_light_watts  # homogeneous light intesity throughout the visible spectrum is assumed

        # photonic calculations
        relative_soret_excitation = 10
        weighted_average_excitation_wavelength = (self.parameters['q (m)']['upper'] + self.parameters['soret (m)']['lower']*relative_soret_excitation) / (relative_soret_excitation+1)
        joules_per_photon = (h*c) / weighted_average_excitation_wavelength
        non_reflected_photons = 0.96 # “SINGLET OXYGEN GENERATION BY PORPHYRINS AND THE PHOTOSENSITIZATION IN LIPOSOMES KINETICS OF 9,lO-DIMETHYLANTHRACENE” by Gross et al., 1993
        self.variables['photon_moles_per_timestep'] = effective_excitation_watts/joules_per_photon*non_reflected_photons/N_A * self.parameters['timestep (s)']

        # singlet oxygen calculations
        '''so_from_light = photons_per_second * molecules_dissolved_oxygen * excitation_constant'''
        mw_molecular_oxygen = periodic_table.O.MW*2*kilo          #mg/mole
        self.variables['dissolved_mo_molar'] = 9/mw_molecular_oxygen    # mg/L -> M, ambient water quality criteria for Dissolved Oxygen, EPA   

        if self.verbose:
            message1 = 'photons per timestep: ', self.variables['photon_moles_per_timestep']
            message2 = 'molecular oxygen molecules: ', sigfigs_conversion(self.variables['dissolved_mo_molar'])
            message4 = 'effective excitation watts: ', sigfigs_conversion(effective_excitation_watts)
            messages = [message1, message2, message4]
            self.messages.extend(messages)
            for message in messages:            
                print(message)
        
        # eps oxidation
        if self.biofilm:
            self.variables['eps_oxidation'] = 1e10 # empirical rate constant that represents a biofilm system from the single cellular kinetic model via a competiting oxidation reaction of EPS versus fatty acids   
        
    def kinetic_calculation(self, calculate_excited_photosensitizer_conc = False, percent_oxidation_threshold = 8):
        self.parameters['percent_oxidation_threshold'] = percent_oxidation_threshold
        # ============= define photosensitizer excitation ==============
        ps = self.variables['photosensitizer_molar']
        if not calculate_excited_photosensitizer_conc:
            e_fraction = (self.variables['photon_moles_per_timestep']*self.variables['volume_proportion']) / (self.variables['photosensitizers']/N_A)
            if e_fraction > 1:
                e_fraction = 1

            self.variables['e_ps_decay_time (s)'] = average(1.5,15)*nano # “Ultrafast excitation transfer and relaxation inlinear and crossed-linear arrays of porphyrins” by Akimoto et al., 1999  ; "Kinetics and efficiency of excitation energy transfer from chlorophylls, their heavy metal-substituted derivatives, and pheophytins to singlet oxygen" by Küpper et al., 2002
            k_ps_rlx = 1/self.variables['e_ps_decay_time (s)']     
            self.variables['ps_excitation (s)'] = 50*femto # an estimated time that is below the detection limit of femto-second laser spectrophotometers ; literature has not been discovered that illuminates this time. # self.parameters['photosensitizer_molar']/(self.variables['photon_moles_per_timestep']*self.variables['volume_proportion']) 
            k_e_ps = 1/self.variables['ps_excitation (s)']
            qy_e = self.photosensitizer['quantum_yield']['value']
            
            photosensitizer = f'{k_e_ps}*{e_fraction}*{qy_e}*ps - {k_ps_rlx}*e_ps'
        else:
            self.variables['e_ps_calculated'] = self.variables['photon_moles_per_timestep']*self.variables['volume_proportion']*self.photosensitizer['quantum_yield']['value']
            if self.variables['e_ps_calculated'] > self.parameters['photosensitizer_molar']:
                self.variables['e_ps_calculated'] = self.parameters['photosensitizer_molar']
                if self.verbose:
                    message = 'excited photosensitizer {} exceeds the photosensitizer concentration {}'.format(self.variables['e_ps_calculated'], self.parameters['photosensitizer_molar'])
                    print(message)
                    self.messages.append(message)
            photosensitizer = self.variables['e_ps_calculated']
            
        # ============== define photobleaching ==============
        k_b_ps = self.photosensitizer['photobleaching_constant (cm^2/J)']['value'] * (self.parameters['watts']/(self.area/centi**2))
        self.variables['hv_photobleaching (s)'] = 1/k_b_ps
        
        # ============== define singlet oxygen generation ==============
        mo = self.variables['dissolved_mo_molar']
        qy = self.variables['so_qy'] = self.photosensitizer['quantum_yield']['value'] * self.photosensitizer['so_specificity']['value'] 
        
        self.variables['e_ps_charge_transfer (s)'] = average(20,)*nano   # “Kinetics and efficiency of excitation energy transfer from chlorophylls, their heavy metal-substituted derivatives, and pheophytins to singlet oxygen” by Küpper et al., 2002 
        k_so = 1/self.variables['e_ps_charge_transfer (s)']
        
        cfu_log_lifetime_slope = (40-10)/(8-4)               # “The role of singlet oxygen and oxygen concentration in photodynamic inactivation of bacteria” by Maisch et al., 2007
        lifetime = cfu_log_lifetime_slope*log(self.bacterial_cfu_ml, 10)
        self.variables['so_decay_time (s)'] = lifetime*micro
        if self.variables['so_decay_time (s)'] < 4*micro:  # a minimum lifetime of 4 microseconds is parameterized for dilute aqueous solution, according to the concensus from literature: e.g. “The role of singlet oxygen and oxygen concentration in photodynamic inactivation of bacteria” by Maisch et al., 2007 
            self.variables['so_decay_time (s)'] = 4*micro  
        k_rlx_so = 1/self.variables['so_decay_time (s)']
        
#         self.variables['so_rise_time (s)'] = 2*micro        # “Time-Resolved Investigations of Singlet Oxygen Luminescence in Water, in Phosphatidylcholine, and in Aqueous Suspensions of Phosphatidylcholine or HT29 Cells” by Baier et al., 2005
        
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
        total_points = int(self.parameters['total_time (s)'] / self.parameters['timestep (s)'])
        self.phrasedml_str = '''
          model1 = model "pdipy_oxidation"
          sim1 = simulate uniform({}, {}, {})
          task1 = run sim1 on model1
          plot "Cytoplasmic oxidation of {}" time vs oxidation
        '''.format(self.parameters['initial_time'], self.parameters['total_time (s)'], total_points, self.parameters['bacterial_specie'])
        
        # ============== execute the model ==============
        tellurium_model = tellurium.loada(self.model)    
        result = tellurium_model.simulate(self.parameters['initial_time'], self.parameters['total_time (s)'], total_points)
        
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
    
        
    def export(self, simulation_name = None, simulation_path = None, exposure_axis = False, figure_title = None, x_label = 'Time (min)', y_label = 'proportion of total', fa_oxidation_proportion = True, display_excitation_proportion = False):
        # parse the simulation results
        x_values = []
        oxidation_y_values = []
        reduction_y_values = []
        excitation_y_values = []
        oxidations = []
        reductions = []
        indices = []
        for index, point in self.result_df.iterrows():      # cytoplasmic oxidation may not be 1:1 correlated with reduction
            x_value = index/minute
            if exposure_axis:
                x_value *= (self.parameters['watts']*minute/self.parameters['surface_area (m^2)'])
            x_values.append(x_value) # sigfigs_conversion(index))
            indices.append(index/hour)
            
            # calculate the oxidation_proportion of the cytoplasmic fatty acids
            oxidation_proportion = point['[ofa]'] / (point['[ofa]']+point['[fa]']) 
            oxidation_y_values.append(oxidation_proportion)
            oxidations.append(oxidation_proportion)
            
            # calculate the percent_reduction of the bacteria
            percent_reduction = oxidation_proportion * (100/self.parameters['percent_oxidation_threshold'])
            reduction_y_values.append(percent_reduction)
            reductions.append(percent_reduction)
            
            # calculate the excitation_proportion of the photosensitizers
            excitation_proportion = point['[e_ps]'] / (point['[ps]']+point['[e_ps]']+point['[b_ps]']) 
            excitation_y_values.append(excitation_proportion)            
        
        # define the simulation_path
        if simulation_path is None:
            if simulation_name is None:
                simulation_name = '-'.join([re.sub(' ', '_', str(x)) for x in [date.today(), 'PDIpy', self.parameters['photosensitizer_selection'], self.parameters['bacterial_specie']]])
            count = 0
            while os.path.exists(simulation_name):
                count += 1
                simulation_name = re.sub('([0-9]+)$', str(count), simulation_name)
                if not re.search('(-[0-9]+$)', simulation_name):
                    simulation_name += f'-{count}'
            simulation_path = os.path.join(os.getcwd(), simulation_name)        
            os.mkdir(simulation_path)
        else:
            count = 0
            while os.path.exists(simulation_path):
                if not re.search('-[0-9]+\..+', simulation_path):
                    simulation_path = re.sub('(\..+)', f'-{count}\..+', simulation_path)
                else:
                    simulation_path = re.sub('-[0-9]+', f'-{count}\..+', simulation_path)
                count += 1
            os.mkdir(simulation_path)
        
        # create and export the OMEX file
        inline_omex = '\n'.join([self.model, self.phrasedml_str])  
        omex_file_name = 'input.omex'
        tellurium.exportInlineOmex(inline_omex, os.path.join(simulation_path, omex_file_name))
        self.result_df.to_csv(os.path.join(simulation_path, 'raw_data.csv'))
        
        # export the processed data
        self.processed_data = pandas.DataFrame(index = indices)
        self.processed_data.index.name = 'time (hr)'
        self.processed_data['oxidation_proportion'] = oxidations
        self.processed_data['percent_reduction'] = reductions
        display(self.processed_data)
        self.processed_data.to_csv(os.path.join(simulation_path, 'processed_data.csv'))
        
        # export the figure
        figure_path = os.path.join(simulation_path, 'output.svg')
        pyplot.rcParams['figure.figsize'] = (11, 7)
        pyplot.rcParams['figure.dpi'] = 150
        figure, ax = pyplot.subplots()
        ax.plot(x_values, reduction_y_values, label = 'Reduction proportion')
        if fa_oxidation_proportion:
            ax.plot(x_values, oxidation_y_values, label = 'Fatty acid oxidation')
        if display_excitation_proportion:
            ax.plot(x_values, excitation_y_values, label = 'Photosensitizer excitation')
        if exposure_axis:
            x_label = 'exposure (J/m^2)'
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        if figure_title is None:
            figure_title = 'PDI oxidation of the {} cytoplasm'.format(self.parameters['bacterial_specie'])
        ax.set_title(figure_title)
        ax.legend()
        figure.savefig(figure_path)
        if self.verbose:
            if not self.jupyter:
                figure.show()
        
        # define a table of parameters
        parameters = {'parameter':[], 'value':[]}
        parameters['parameter'].append('simulation_path')
        parameters['value'].append(simulation_path)
        for parameter, value in self.parameters.items():
            parameters['parameter'].append(parameter)
            if isnumber(value):
                parameters['value'].append(sigfigs_conversion(value, 5))
            else:
                parameters['value'].append(value)            
        parameters_table = pandas.DataFrame(parameters)
        parameters_path = os.path.join(simulation_path, 'parameters.csv')
        parameters_table.to_csv(parameters_path)
        if self.verbose:
            if self.jupyter:
                display(parameters_table)
            else:
                print(parameters_table)
        
        # define a table of variables
        variables = {'variable':[], 'value':[]}
        variables['variable'].append('simulation_path')
        variables['value'].append(simulation_path)
        for variable, value in self.variables.items():
            variables['variable'].append(variable)
            if isnumber(value):
                variables['value'].append(sigfigs_conversion(value, 5))
            else:
                variables['value'].append(value)
        variables_table = pandas.DataFrame(variables)        
        variables_path = os.path.join(simulation_path, 'variables.csv')
        variables_table.to_csv(variables_path)
        if self.verbose:
            if self.jupyter:
                display(variables_table)
            else:
                print(variables_table)     
                
        return self.processed_data
        
        
    def data_processing(self, target_reduction = None, target_time = None):
        quantity = unit = None
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
            
        print(message)
        self.messages.append(message)
                
        return quantity, unit