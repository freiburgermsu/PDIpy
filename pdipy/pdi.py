from scipy.constants import femto, pico, angstrom, nano, micro, milli, centi, kilo, liter, N_A, h, c, minute, hour
from math import pi, cos, ceil, log, pow
from chemicals.elements import periodic_table
from numpy import array, log10
from matplotlib import pyplot
from hillfit import HillFit
from datetime import date
from sigfig import round
from glob import glob
from chemw import ChemMW
import tellurium
import pandas
import json, sys, re, os, sys

pandas.set_option('max_colwidth', None)

# define useful functions
def sigfigs_conversion(num, sigfigs_in = 2):
    return round(num, sigfigs=sigfigs_in, notation = 'sci')

def average(num_1, num_2 = None):
    if isnumber(num_1): 
        if isnumber(num_2):
            numbers = [num_1, num_2]
            return sum(numbers) / len(numbers)
        else:
            return num_1
    elif type(num_1) is list:
        summation = total = 0
        for num in num_1:
            if num is not None:
                summation += num
                total += 1
        if total > 0:
            return summation/total
        raise None # ValueError('The arguments must be numbers or a list of numbers')
    elif isnumber(num_2):
        return num_2
    else:
        raise None # ValueError('The arguments must be numbers or a list of numbers')
    
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

class PDI():
    def __init__(self, 
                 total_time: float,              # The total simulation time
                 solution_dimensions: dict = {}, # defines dimensions of the simulated solution 
                 surface_system: bool = False,   # specifies whether cross-linked photosensitizer will be simulated
                 well_count: float = 24,        # The petri dish well count that will be simulated, which begets default dimensions of the simulated solution
                 timestep: float = 3,            # The simulation timestep value, which subtly affects the simulation predictions
                 verbose: bool = False,     
                 jupyter: bool = False
                 ):
        # define base organizations
        self.parameters = {}
        self.variables = {}
        self.results = {}
        self.defined_model = {}
        self.messages = []
        self.surface_system = surface_system
        self.verbose = verbose
        self.jupyter = jupyter
        self.chem_mw = ChemMW()
        
        # initial parameters
        self.parameters['so_diffusion (m)'] = 80 * nano       # Moan1990 
#         self.parameters['oxidation_angle'] = 5           # degrees
        self.parameters['root_path'] = os.path.dirname(__file__)
        
        # identify options
        self.bacteria = [re.search('(?<=bacteria\\\\)(.+)(?=\.json)', x).group() for x in glob(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'*.json'))]
        self.light_parameters = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'light_source.json'), encoding = 'utf-8'))
        self.photosensitizers = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'photosensitizers.json'), encoding = 'utf-8'))        
        
        # time conditions
        self.parameters['total_time (s)'] = int(total_time * minute)
        self.parameters['timestep (s)'] = timestep * minute
        
        # define the solution 
        self.solution = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'wells.json')))[well_count]
        if solution_dimensions != {}:
            for key, value in solution_dimensions.items():
                self.solution[key] = value
            
        operational_depth = 0.7  # proportion of well volume that is filled with solution
        self.parameters['solution_depth (m)'] = self.solution['depth (cm)']*operational_depth*centi
        self.parameters['solution_area (m^2)'] = self.solution['area (cm^2)']*centi**2
        self.parameters['solution_volume (m^3)'] = self.parameters['solution_depth (m)']*self.parameters['solution_area (m^2)']
            
        # completion of the function
        self.defined_model.update({'define_system':True})
        
                
    def define_bacterium(self, 
                         bacterial_specie: str = None,    # specifies one of the bacteria in the parameters directory of PDI to simulate
                         bacterial_characteristics: dict = {},      # passes a custom dictionary of characteristics of the simulated bacterium, which can refine characteristics from the bacterial_specie argument
                         bacterial_cfu_ml: float = 1E6,   # specifies the solution concentration of the simulated bacterium for solution simulations
                         biofilm: bool = False,          # specifies whether a biofilm simulation will be conducted
                         ):    
        def shell_volume(outer_radi, inner_radi_2):
            return (4*pi/3)*(outer_radi**3 - inner_radi_2**3)  
        def shell_radius(volume):
            return pow((volume*3)/(4*pi), 1/3)
    
        self.bacterial_cfu_ml = bacterial_cfu_ml
        self.biofilm = biofilm
        
        # load the bacterial parameters
        self.bacterium = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', 'S_aureus.json')))
        if bacterial_specie in self.bacteria:
            self.parameters['bacterial_specie'] = bacterial_specie
            self.bacterium = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'{bacterial_specie}.json')))
        if bacterial_characteristics != {}:
            for key, value in bacterial_specie.items():
                if key == 'name':
                    self.parameters['bacterial_specie'] = value
                else:
                    for key2, value2 in value.items():
                        self.bacterial_specie[key][key2] = value2 
       
        self.membrane_chemicals = self.bacterium['membrane_chemicals']
                                          
        # the fatty acid concentrations of the cytoplasmic membrane are determined
        self.variables['fa_g/L_conc'] = total_proportion = 0
        weighed_mws = []
        for chemical in self.membrane_chemicals:
            if re.search('FA', chemical):
                total_proportion += self.membrane_chemicals[chemical]['proportion']['value']
        for chemical in self.membrane_chemicals:
            if re.search('FA', chemical):
                fa_density_proportion = self.membrane_chemicals[chemical]['density (g/L)']['value'] * self.membrane_chemicals[chemical]['proportion']['value']
                self.variables['fa_g/L_conc'] += fa_density_proportion
                mw = average(self.membrane_chemicals[chemical]['mw'])
                weighed_mws.append(mw*fa_density_proportion/total_proportion)
        self.variables['fa_molar'] = self.variables['fa_g/L_conc'] / average(weighed_mws)
        
        # calculate the mass proportion of fatty acids in the membrane
        self.variables['cell_radius (m)'] = shell_radius(self.bacterium['cell_volume (fL)']['value']*femto/liter)
        membrane_inner_radius = self.variables['cell_radius (m)'] - self.bacterium['membrane_thickness (nm)']['value']*nano
        self.variables['membrane_volume'] = shell_volume(self.variables['cell_radius (m)'], membrane_inner_radius)
        
        fa_g = self.variables['fa_g/L_conc']*self.variables['membrane_volume']
        self.variables['fa_mass_proportion'] = fa_g/self.bacterium['cell_mass (pg)']['value']*pico    # What is the purpose of this proportion?
        
        # completion of the function
        self.defined_model.update({'define_bacterium':True})
                

    def define_photosensitizer(self, 
                               photosensitizer: str = 'A3B_4Zn',     # specify which photosensitizer from the predefined options will be simulated
                               photosensitizer_characteristics: dict = {},   # Custom specifications of the simulation photosensitizer, which can be used to refine the parameterized photosensitizer
                               photosensitizer_molar: float = None,    # specify the photosensitizer molar concentration for solution simulations
                               photosensitizer_g: float = 90e-9,        # specify the mass of photosensitizer for cross-linked simulations
                               cross_linked_sqr_m: float = 0.0191134,  # define the cross-linked square-meters for surface simulations
                               molecular_proportion: float = None     #!!! What is the purpose or intention of this argument?
                               ):
        def cylinder_volume(radius, height):
            return (radius**2)*pi*height
        
        # define the photosensitizing surface
        self.area = self.parameters['solution_area (m^2)']
        if self.surface_system:            
            self.area = self.parameters['cross_linked_surface_area (m^2)'] = cross_linked_sqr_m
            self.parameters['photosensitizer (g/m^2)'] = photosensitizer_g/cross_linked_sqr_m
            
        # load and refine photosensitizer parameters
        self.photosensitizer = self.photosensitizers['A3B_4Zn']
        if photosensitizer in self.photosensitizers:
            self.parameters['photosensitizer_selection'] = photosensitizer
            self.photosensitizer = self.photosensitizers[photosensitizer] 
        else: #if photosensitizer not in self.photosensitizers:
            error = f'--> ERROR: The photosensitizer {photosensitizer} is neither a predefined option nor a customized dictionary.'
            self.messages.append(error)
            raise ImportError(error)    
            
        if photosensitizer_characteristics != {}:
            for key, value in photosensitizer_characteristics.items():
                if key == 'name':
                    self.parameters['photosensitizer_selection'] = value
                else:
                    for key2, value2 in value.items():
                        self.photosensitizer[key][key2] = value2
                                                       
        # classify the photonic absorption characteristics of the photosensitizer
        upper_soret = self.photosensitizer['soret (nm)']['value'][1]*nano
        lower_soret = self.photosensitizer['soret (nm)']['value'][0]*nano
        upper_q = self.photosensitizer['q (nm)']['value'][1]*nano
        lower_q = self.photosensitizer['q (nm)']['value'][0]*nano
        self.parameters['soret (m)'] = {'upper': upper_soret, 'lower': lower_soret}
        self.parameters['q (m)'] = {'upper': upper_q, 'lower': lower_q}

        # calculate physical dimensions of the photosensitizer
        photosensitizer_length = self.photosensitizer['dimensions']['length (A)']*angstrom
        photosensitizer_width = self.photosensitizer['dimensions']['width (A)']*angstrom
        photosensitizer_depth = self.photosensitizer['dimensions']['depth (A)']*angstrom
        photosensitizer_shape = self.photosensitizer['dimensions']['shape']
        photosensitizer_mw = self.chem_mw.mass(self.photosensitizer['formula']['value'])
                                          
        if not self.surface_system:
            if self.variables['photosensitizer_molar'] is None:
                error = '--> ERROR: A molar photosensitizer concentration must be defined for systems of dissolved photosensitizers.'
                self.messages.append(error)
                raise ImportError(error)
            self.variables['photosensitizer_molar'] = photosensitizer_molar
            self.variables['photosensitizers'] = (self.variables['photosensitizer_molar']*N_A) * (self.parameters['solution_volume (m^3)']/liter)
        else:
            self.variables['photosensitizers'] = self.parameters['photosensitizer (g/m^2)'] / (photosensitizer_mw/N_A) * self.parameters['cross_linked_surface_area (m^2)']
            self.parameters['solution_volume (m^3)'] = self.parameters['cross_linked_surface_area (m^2)']*photosensitizer_length
            self.variables['photosensitizer_molar'] = self.variables['photosensitizers']/N_A / (self.parameters['solution_volume (m^3)']/liter)  

        # calculate the volume proportion that is constituted by the photosensitizer, in the volume where the photosensitizer is present        
        if photosensitizer_shape == 'disc':
            self.variables['molecular_volume (m^3)'] = cylinder_volume(photosensitizer_length/2, photosensitizer_depth)
        else:
            error = f'--> ERROR: The volume formula for the {photosensitizer_shape} shape is not available.'
            self.messages.append(error)
            raise ValueError(error)
            
        molecules_volume = self.variables['photosensitizers'] * self.variables['molecular_volume (m^3)']
        self.variables['volume_proportion'] = molecules_volume / self.parameters['solution_volume (m^3)']
            
        # print calculated content
        if self.verbose:
            messages = [
                    'The photosensitizer dimensions as a {} = {} m x {} m x {} m.'.format(photosensitizer_shape, photosensitizer_length, photosensitizer_width, photosensitizer_depth), 
                    'Photosensitizer volume = {} m\N{superscript three}'.format(sigfigs_conversion(self.variables['molecular_volume (m^3)'])), 
                    'The volume proportion of {} photosensitizers = ({} m\N{superscript three} of photosensitizer)/({} m\N{superscript three} of solution) = {}'.format(
                            sigfigs_conversion(self.variables['photosensitizers']), 
                            sigfigs_conversion(molecules_volume), 
                            sigfigs_conversion(self.parameters['solution_volume (m^3)']), 
                            sigfigs_conversion(self.variables['volume_proportion'])
                            )
                    ]
            
            self.messages.extend(messages)
            for message in messages:            
                print(message)
                
        # completion of the function
        self.defined_model.update({'define_photosensitizer':True})
                

    def define_light(self, 
                     light_source: str = None,             # specifies a light source from the predefined options
                     light_characteristics: dict = {},      # specifies custom characteristics of the light source, in addition to or substitute of a predefined option
                     measurement: dict = None,              # provides the measurement, in the proper respective units, for the photonic intensity of the light source
                     ):
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
        if 'irradiance' in measurement: # mW / cm^2
            self.parameters['watts'] = (measurement['irradiance']*milli/centi**2) * self.area
        elif 'exposure' is not measurement:  # J / cm^2
            simulation_time = simulation_time * minute
            self.parameters['watts'] = (measurement['exposure']/centi**2) / simulation_time * self.area
        elif 'lux' is not measurement: # lumen / m^2
            self.parameters['watts'] = measurement['lux'] / lumens_per_watt * self.area
        elif 'lumens' is not measurement: # lumen
            self.parameters['watts'] = measurement['lumens'] / lumens_per_watt
        else:
            error = '--> ERROR: The light source has an unrecognized dimension.'
            self.messages.append(error)
            raise ValueError(error)
            
        # completion of the function
        self.defined_model.update({'define_light':True})
        
    def _completed_functions(self,):
        missing_functions = []
        for function in self.defined_model:
            if not self.defined_model[function]:
                missing_functions.append(function)
            if missing_functions != []:
                error = f'--> ERROR: The {function} function must be defined before the simulation can execute.'
                self.messages.append(error)
                sys.exit(error)
            
    def _singlet_oxygen_calculations(self,):
        # ensure that all prior functions have completed within this instance
        self._completed_functions()
                
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
            messages = [
                    'photons per timestep: ', self.variables['photon_moles_per_timestep'], 
                    'molecular oxygen molecules: ', sigfigs_conversion(self.variables['dissolved_mo_molar']), 
                    'effective excitation watts: ', sigfigs_conversion(effective_excitation_watts)
                    ]
            
            self.messages.extend(messages)
            for message in messages:            
                print(message)
        
        # eps oxidation
        if self.biofilm:
            self.variables['eps_oxidation'] = self.bacterium['eps_oxidation_rate_constant']['value']
            self.variables['bacterial_biofilm_mass_proportion'] = 
        
    def _kinetic_calculation(self, calculate_excited_photosensitizer_conc):        
        # ============= define photosensitizer excitation ==============
        ps = self.variables['photosensitizer_molar']
        if not calculate_excited_photosensitizer_conc:
            e_fraction = (self.variables['photon_moles_per_timestep']*self.variables['volume_proportion']) / (self.variables['photosensitizers']/N_A)
            if e_fraction > 1:
                e_fraction = 1

            self.variables['e_ps_decay_time (s)'] = average(1.5,15)*nano # “Ultrafast excitation transfer and relaxation inlinear and crossed-linear arrays of porphyrins” by Akimoto et al., 1999  ; "Kinetics and efficiency of excitation energy transfer from chlorophylls, their heavy metal-substituted derivatives, and pheophytins to singlet oxygen" by Küpper et al., 2002
            k_ps_rlx = 1/self.variables['e_ps_decay_time (s)']     
            self.variables['ps_excitation (s)'] = 50*femto # an estimated time that is below the detection limit of femto-second laser spectrophotometers; specific measurements have not been discovered.
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
        
        # ============== define fatty acid and EPS oxidation ==============
        k_fa = self.variables['k_fa'] = 7.69E2 * self.variables['fa_g/L_conc'] # https://www.jstage.jst.go.jp/article/jos/68/1/68_ess18179/_pdf/-char/ja
        fa = self.variables['fa_molar']
        k_fa_reduction = (self.bacterial_cfu_ml/1E6)**0.1
        k_fa /= k_fa_reduction
        
        biofilm = eps = ''
        if self.biofilm:
            biofilm = 'so + eps => o_eps + mo; {}*so*eps'.format(self.variables['eps_oxidation'])            
            eps = 'eps = {}'.format((fa/self.variables['fa_mass_proportion'])/self.bacterium['cellular_dry_mass_proportion_biofilm']['value'])

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
            {eps}
            
            # calculate the oxidation proportion
            oxidation := o_fa / (o_fa + fa);
            
          end
        ''')       
        # ============== SED-ML plot ==============
        total_points = int(self.parameters['total_time (s)'] / self.parameters['timestep (s)'])
        self.phrasedml_str = '''
          model1 = model "pdipy_oxidation"
          sim1 = simulate uniform(0, {}, {})
          task1 = run sim1 on model1
          plot "Cytoplasmic oxidation of {}" time vs oxidation
        '''.format(self.parameters['total_time (s)'], total_points, self.parameters['bacterial_specie'])
        
        # ============== execute the model ==============
        tellurium_model = tellurium.loada(self.model)    
        result = tellurium_model.simulate(0, self.parameters['total_time (s)'], total_points)
        
        # ============== process the data ==============
        self.raw_data = pandas.DataFrame(result)
        self.raw_data.index = self.raw_data[0]
        del self.raw_data[0]
        self.raw_data.index.name = 'Time (s)'
        
        if self.biofilm:
            self.raw_data.columns = ['[ps]', '[e_ps]', '[b_ps]', '[mo]', '[so]', '[fa]', '[ofa]', '[oeps]']
        else:
            self.raw_data.columns = ['[ps]', '[e_ps]', '[b_ps]', '[mo]', '[so]', '[fa]', '[ofa]']
        
        if self.verbose:
            message1 = tellurium_model.getCurrentAntimony()
            message2 = '\nCurrent integrator:', '\n', tellurium_model.integrator
            message3 = 'k_fa reduction factor:', k_fa_reduction

            messages = [message1, message2, message3]
            self.messages.extend(messages)
            for message in messages:            
                print(message)
            
            if self.jupyter:
                display(self.raw_data)
            else:
                print(self.raw_data)
    
    def simulate(self, figure_title = None, y_label = 'log10', exposure_axis = False, display_fa_oxidation = False, calculate_excited_photosensitizer_conc = False):
        def approach_1(xs, limit, inactivation_y_values, top_increment, top_increment_change, count, relative_to_1 = 'greater'):
            if relative_to_1 == 'greater':
                while inactivation_y_values[-1] > limit:
                    inactivation_y_values = list(eval(f'{self.hf.bottom} + ({self.hf.top}-{top_increment}-{self.hf.bottom})*xs**({self.hf.nH}+{nH_change}) / (({self.hf.ec50}+{ec50_change})**({self.hf.nH}+{nH_change}) + xs**({self.hf.nH}+{nH_change}))'))
                    top_increment += top_increment_change
                    count += 1
            elif relative_to_1 == 'lesser':
                while inactivation_y_values[-1] < limit:
                    inactivation_y_values = list(eval(f'{self.hf.bottom} + ({self.hf.top}-{top_increment}-{self.hf.bottom})*xs**({self.hf.nH}+{nH_change}) / (({self.hf.ec50}+{ec50_change})**({self.hf.nH}+{nH_change}) + xs**({self.hf.nH}+{nH_change}))'))
                    top_increment += top_increment_change
                    count += 1
            return count, inactivation_y_values, top_increment
        
        # calculate the kinetics of the simulation 
        self._singlet_oxygen_calculations()
        self._kinetic_calculation(calculate_excited_photosensitizer_conc)
                    
        # parse the simulation results
        x_values = []
        oxidation_y_values = []
        excitation_y_values = []
        first = True
        for index, point in self.raw_data.iterrows():      # cytoplasmic oxidation may not be 1:1 correlated with reduction
            if first:   # The inital point of 0,0 crashes the regression of HillFit, and thus it is skipped
                first = False
                continue
            x_value = index/hour
            if exposure_axis:
                x_value *= (self.parameters['watts']*hour/self.parameters['cross_linked_surface_area (m^2)'])
            x_values.append(x_value) # sigfigs_conversion(index))
            
            # calculate the oxidation_proportion of the cytoplasmic fatty acids
            oxidation_proportion = point['[ofa]'] / (point['[ofa]']+point['[fa]']) 
            oxidation_y_values.append(oxidation_proportion)

        # determine the regression equation for the fitted curve via the Hill equation
        xs = array(x_values)
        ys = array(oxidation_y_values)
        self.hf = HillFit(xs, ys)
        self.hf.fitting(x_label = 'time (hr)', y_label = 'oxidation proportion', view_figure = False)
        
        # define and refine the fitted Hill equation parameters
        increments = [0.01*self.hf.top, -0.00001*self.hf.top, 0.000001*self.hf.top, -0.0000001*self.hf.top, 0.00000001*self.hf.top]
        top_increment = 0.01*self.hf.top
        ec50_change = -.76*self.hf.ec50
        nH_change = 1.24*self.hf.nH
        limit = 1
        if self.biofilm:
            increments = [0.01*self.hf.top, -0.00001*self.hf.top, 0.000001*self.hf.top]
            top_increment = 0.01*self.hf.top
            ec50_change = -.8*self.hf.ec50
            nH_change = -0.1*self.hf.nH
            limit = 1-10**-7

        inactivation_y_values = list(eval(f'{self.hf.bottom} + ({self.hf.top}-{top_increment}-{self.hf.bottom})*xs**({self.hf.nH}+{nH_change}) / (({self.hf.ec50}+{ec50_change})**({self.hf.nH}+{nH_change}) + xs**({self.hf.nH}+{nH_change}))'))
        count = 0
        relative_to_1 = 'greater'
        for top_increment_change in increments:
            count, inactivation_y_values, top_increment = approach_1(xs, limit, inactivation_y_values, top_increment, top_increment_change, count, relative_to_1)
            print('refinement loop: ', count)
            if relative_to_1 == 'greater':
                relative_to_1 = 'lesser'
            else:
                relative_to_1 = 'greater'

        if self.verbose:
            message = f'The oxidation data was distilled into inactivation data in {count} loops'
            print(message)
            self.messages.append(message)
                
        # create the DataFrame of processed data        
        log_oxidation = -log10(1-array(oxidation_y_values))
        log_inactivation = -log10(1-array(inactivation_y_values))
        self.processed_data = pandas.DataFrame(index = list(xs))
        index_label = 'time (hr)'
        if exposure_axis:
            index_label = 'exposure (J/cm\N{superscript two})'
        self.processed_data.index.name = index_label
        self.processed_data['oxidation'] = oxidation_y_values
        self.processed_data['inactivation'] = inactivation_y_values
        self.processed_data['log-oxidation'] = log_oxidation
        self.processed_data['log-inactivation'] = log_inactivation
        
        # create the data figure
        pyplot.rcParams['figure.figsize'] = (11, 7)
        pyplot.rcParams['figure.dpi'] = 150
        self.figure, ax = pyplot.subplots()
        ax.plot(xs, log_inactivation, label = 'Inactivation')
        if display_fa_oxidation:
            ax.plot(xs, log_oxidation, label = 'Oxidation')
        ax.set_ylabel(y_label)
        ax.set_xlabel(index_label)
        if figure_title is None:
            figure_title = 'Cytoplasmic oxidation and inactivation of {} via PDI'.format(self.parameters['bacterial_specie'])
        ax.set_title(figure_title)
        ax.legend()    
        
    def export(self, export_name = None, export_directory = None):       
        # define the simulation_path
        if export_directory is None:
            export_directory = os.getcwd()
        elif not os.path.exists(export_directory):
            error = '--> ERROR: The provided directory does not exist'
            print(error)
            self.messages.append(error)

        if export_name is None:
            export_name = '-'.join([re.sub(' ', '_', str(x)) for x in [date.today(), 'PDIpy', self.parameters['photosensitizer_selection'], self.parameters['bacterial_specie']]])
            
        count = 0
        export_path = os.path.join(export_directory, export_name)
        while os.path.exists(export_path):
            if not re.search('(-[0-9]+$)', export_path):
                export_path += f'-{count}'   
            else:
                export_path = re.sub('([0-9]+)$', str(count), export_name)
            count += 1
            
        os.mkdir(export_path)
            
        # create and export the OMEX file
        inline_omex = '\n'.join([self.model, self.phrasedml_str])  
        omex_file_name = 'input.omex'
        tellurium.exportInlineOmex(inline_omex, os.path.join(export_path, omex_file_name))
        self.raw_data.to_csv(os.path.join(export_path, 'raw_data.csv'))
        
        # export the processed data and the regression content
        self.processed_data.to_csv(os.path.join(export_path, 'processed_data.csv'))
        if not self.biofilm:
            self.hf.export(export_path, 'hillfit-regression')
        
        # export the figure
        self.figure.savefig(os.path.join(export_path, 'inactivation.svg'))
        if self.verbose:
            if not self.jupyter:
                self.figure.show()
        
        # define and export a table of parameters
        parameters = {'parameter':[], 'value':[]}
        parameters['parameter'].append('simulation_path')
        parameters['value'].append(export_path)
        for parameter, value in self.parameters.items():
            parameters['parameter'].append(parameter)
            if isnumber(value):
                parameters['value'].append(sigfigs_conversion(value, 5))
            else:
                parameters['value'].append(value)      
                
        parameters_table = pandas.DataFrame(parameters)
        parameters_table.to_csv(os.path.join(export_path, 'parameters.csv'))
        if self.verbose:
            if self.jupyter:
                display(parameters_table)
            else:
                print(parameters_table)
        
        # define and export a table of variables
        variables = {'variable':[], 'value':[]}
        variables['variable'].append('simulation_path')
        variables['value'].append(export_path)
        for variable, value in self.variables.items():
            variables['variable'].append(variable)
            if isnumber(value):
                variables['value'].append(sigfigs_conversion(value, 5))
            else:
                variables['value'].append(value)
                
        variables_table = pandas.DataFrame(variables)   
        variables_table.to_csv(os.path.join(export_path, 'variables.csv'))
        if self.verbose:
            if self.jupyter:
                display(variables_table)
            else:
                print(variables_table)     
        
        
    def data_parsing(self, log_reduction = None, target_time = None):
        value = unit = None
        if isnumber(log_reduction):
            if log_reduction > self.processed_data['inactivation'].iloc[-1]:
                message = '--> ERROR: The inquired reduction is never reached.'
            else:
                for index, point in self.processed_data.iterrows():
                    if point['inactivation'] >= log_reduction:
                        value = index
                        unit = 'hours'
                        message = f'{unit} to target: {value}'
                        break
        elif isnumber(target_time):
            if target_time > self.processed_data.index[-1]:
                message = '--> ERROR: The inquired time is never reached.'
            else:
                for index, point in self.processed_data.iterrows():
                    if index >= target_time:
                        value = point['log-inactivation']
                        unit = 'log10-inactivation'
                        message = '{} at {}: {}'.format(unit, target_time, value)
                        break
        else:
            message = '--> ERROR: Neither the target_time nor the target_reduction are parameterized as numbers.'
            
        print(message)
        self.messages.append(message)    
            
        return value, unit