from scipy.constants import femto, pico, angstrom, nano, micro, milli, centi, liter, N_A, h, c, minute, hour
from math import pi, cos, ceil
from numpy import array, log10, diff
from matplotlib import pyplot
from hillfit import HillFit
from sigfig import round
from glob import glob
from chemw import ChemMW
import tellurium
import pandas
import warnings, json, re, os, sys

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
    remainder = re.sub('[0-9e.-]', '', str(num))
    if remainder != '':
        return False
    try:
        float(num)
    except:
        try:
            int(num)
        except:
            return False
    return True

class PDI():
    def __init__(self, 
                 total_time: float,              # The total simulation time
                 solution_dimensions: dict = {}, # defines dimensions of the simulated solution 
                 surface_system: bool = False,   # specifies whether cross-linked photosensitizer will be simulated
                 well_count: str = 24,        # The petri dish well count that will be simulated, which begets default dimensions of the simulated solution
                 timestep: float = 3,            # The simulation timestep value, which subtly affects the simulation predictions
                 verbose: bool = False,     
                 jupyter: bool = False
                 ):
        # define base organizations
        self.parameters = {}
        self.variables = {}
        self.results = {}
        self.paths = {}
        self.defined_model = {}
        self.messages = []
        self.verbose = verbose
        self.jupyter = jupyter
        self.chem_mw = ChemMW(printing = False)
        
        # initial parameters
        self.parameters['surface_system'] = surface_system
        self.parameters['so_diffusion_m'] = 80 * nano       # Moan1990 
#         self.parameters['oxidation_angle'] = 5           # degrees
        self.parameters['root_path'] = os.path.dirname(__file__)
        
        # identify options
        self.bacteria = [re.search('(?<=bacteria\\\\)(.+)(?=\.json)', x).group() for x in glob(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'*.json'))]
        self.light_parameters = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'light_source.json'), encoding = 'utf-8'))
        self.photosensitizers = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'photosensitizers.json'), encoding = 'utf-8'))        
        
        # time conditions
        self.parameters['total_time_s'] = int(total_time * minute)
        self.parameters['timestep_s'] = timestep * minute
        
        # define the solution 
        well_count = str(well_count)
        self.solution = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'wells.json')))[well_count]
        if solution_dimensions != {}:
            for key, value in solution_dimensions.items():
                self.solution[key] = value
            
        operational_depth = 0.7  # proportion of well volume that is filled with solution
        self.parameters['solution_depth_m'] = self.solution['depth_cm']*operational_depth*centi
        self.parameters['solution_sqr_m'] = self.solution['area_sqr_cm']*centi**2
        self.parameters['solution_cub_m'] = self.parameters['solution_depth_m']*self.parameters['solution_sqr_m']
            
        # completion of the function
        self.defined_model.update({'define_system':True})
        
                
    def _define_bacterium(self, 
                         bacterial_specie: str = None,    # specifies one of the bacteria in the parameters directory of PDI to simulate
                         bacterial_characteristics: dict = {},      # passes a custom dictionary of characteristics of the simulated bacterium, which can refine characteristics from the bacterial_specie argument
                         bacterial_cfu_ml: float = 1E6,   # specifies the solution concentration of the simulated bacterium for solution simulations
                         biofilm: bool = False,          # specifies whether a biofilm simulation will be conducted
                         ):    
        def shell_volume(outer_radi, inner_radi_2):
            return (4*pi/3)*(outer_radi**3 - inner_radi_2**3)  
        def shell_radius(volume):
            return pow((volume*3)/(4*pi), 1/3)
    
        self.parameters['bacterial_cfu_ml'] = bacterial_cfu_ml
        self.parameters['biofilm'] = biofilm
        
        # load the bacterial parameters
        self.bacterium = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', 'S_aureus.json')))
        if bacterial_specie in self.bacteria:
            self.parameters['bacterial_specie'] = bacterial_specie
            self.bacterium = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'{bacterial_specie}.json')))
        if bacterial_characteristics != {}:
            for key, value in bacterial_characteristics.items():
                if key == 'name':
                    self.parameters['bacterial_specie'] = value
                else:
                    for key2, value2 in value.items():
                        self.bacterium[key][key2] = value2 
       
        # the fatty acid concentrations of the cytoplasmic membrane are determined
        self.variables['fa_gL_conc'] = total_proportion = 0
        weighed_mws = []
        for chemical in self.bacterium['membrane_chemicals']:
            if re.search('FA', chemical):
                total_proportion += self.bacterium['membrane_chemicals'][chemical]['proportion']['value']
        for chemical in self.bacterium['membrane_chemicals']:
            if re.search('FA', chemical):
                fa_density_proportion = self.bacterium['membrane_chemicals'][chemical]['density_gL']['value'] * self.bacterium['membrane_chemicals'][chemical]['proportion']['value']
                self.variables['fa_gL_conc'] += fa_density_proportion
                mw = average([self.chem_mw.mass(formula) for formula in self.bacterium['membrane_chemicals'][chemical]['formula']])
                weighed_mws.append(mw*fa_density_proportion/total_proportion)
        self.variables['fa_molar'] = self.variables['fa_gL_conc'] / average(weighed_mws)
        
        # calculate the mass proportion of fatty acids in the membrane
        self.variables['cell_radius_m'] = shell_radius(self.bacterium['cell_volume_fL']['value']*femto/liter)
        membrane_inner_radius = self.variables['cell_radius_m'] - self.bacterium['membrane_thickness_nm']['value']*nano
        self.variables['membrane_cub_m'] = shell_volume(self.variables['cell_radius_m'], membrane_inner_radius)
        
        fa_g = self.variables['fa_gL_conc']*self.variables['membrane_cub_m']
        self.variables['fa_mass_proportion'] = fa_g/self.bacterium['cell_mass_pg']['value']*pico    # What is the purpose of this proportion?
        
        # completion of the function
        self.defined_model.update({'define_bacterium':True})
                

    def _define_photosensitizer(self, 
                               photosensitizer: str = 'A3B_4Zn',     # specify which photosensitizer from the predefined options will be simulated
                               photosensitizer_characteristics: dict = {},   # Custom specifications of the simulation photosensitizer, which can be used to refine the parameterized photosensitizer
                               photosensitizer_molar: float = None,    # specify the photosensitizer molar concentration for solution simulations
                               photosensitizer_g: float = 90e-9,        # specify the mass of photosensitizer for cross-linked simulations
                               cross_linked_sqr_m: float = 0.0191134,  # define the cross-linked square-meters for surface simulations
                               # molecular_proportion: float = None     #!!! What is the purpose or intention of this argument?
                               ):
        def cylinder_volume(radius, height):
            return (radius**2)*pi*height
        
        # define the photosensitizing surface
        self.area = self.parameters['solution_sqr_m']
        if self.parameters['surface_system']:            
            self.area = self.parameters['cross_linked_surface_sqr_m'] = cross_linked_sqr_m
            self.parameters['photosensitizer_g_sqr_m'] = photosensitizer_g/cross_linked_sqr_m
            
        # load and refine photosensitizer parameters
        self.photosensitizer = self.photosensitizers['A3B_4Zn']
        if photosensitizer in self.photosensitizers:
            self.parameters['photosensitizer_selection'] = photosensitizer
            self.photosensitizer = self.photosensitizers[photosensitizer] 
        else: #if photosensitizer not in self.photosensitizers:
            error = f'--> ERROR: The photosensitizer {photosensitizer} is neither a predefined option nor a customized dictionary.'
            self.messages.append(error)
            raise TypeError(error)    
            
        if photosensitizer_characteristics != {}:
            for key, value in photosensitizer_characteristics.items():
                if key == 'name':
                    self.parameters['photosensitizer_selection'] = value
                else:
                    for key2, value2 in value.items():
                        self.photosensitizer[key][key2] = value2
                                                       
        # classify the photonic absorption characteristics of the photosensitizer
        upper_soret = self.photosensitizer['soret_nm']['value'][1]*nano
        lower_soret = self.photosensitizer['soret_nm']['value'][0]*nano
        upper_q = self.photosensitizer['q_nm']['value'][1]*nano
        lower_q = self.photosensitizer['q_nm']['value'][0]*nano
        self.parameters['soret_m'] = {'upper': upper_soret, 'lower': lower_soret}
        self.parameters['q_m'] = {'upper': upper_q, 'lower': lower_q}
        self.parameters['excitation_range_m'] = self.parameters['q_m']['upper'] - self.parameters['soret_m']['lower']

        # calculate physical dimensions of the photosensitizer
        photosensitizer_length = self.photosensitizer['dimensions']['length_A']*angstrom
        photosensitizer_width = self.photosensitizer['dimensions']['width_A']*angstrom
        photosensitizer_depth = self.photosensitizer['dimensions']['depth_A']*angstrom
        photosensitizer_shape = self.photosensitizer['dimensions']['shape']
        photosensitizer_mw = self.chem_mw.mass(self.photosensitizer['formula']['value'])
                                          
        if not self.parameters['surface_system']:
            if photosensitizer_molar is None:
                error = '--> ERROR: A molar photosensitizer concentration must be defined for systems of dissolved photosensitizers.'
                self.messages.append(error)
                raise TypeError(error)
            self.variables['photosensitizer_molar'] = photosensitizer_molar
            self.variables['photosensitizers'] = (self.variables['photosensitizer_molar']*N_A) * (self.parameters['solution_cub_m']/liter)
        else:
            self.variables['photosensitizers'] = self.parameters['photosensitizer_g_sqr_m'] / (photosensitizer_mw/N_A) * self.parameters['cross_linked_surface_sqr_m']
            self.parameters['solution_cub_m'] = self.parameters['cross_linked_surface_sqr_m']*photosensitizer_length
            self.variables['photosensitizer_molar'] = self.variables['photosensitizers']/N_A / (self.parameters['solution_cub_m']/liter)  

        # calculate the volume proportion that is constituted by the photosensitizer, in the volume where the photosensitizer is present        
        if photosensitizer_shape == 'disc':
            self.variables['molecular_volume_cub_m'] = cylinder_volume(photosensitizer_length/2, photosensitizer_depth)
        else:
            error = f'--> ERROR: The volume formula for the {photosensitizer_shape} shape is not available.'
            self.messages.append(error)
            raise ValueError(error)
            
        molecules_volume = self.variables['photosensitizers'] * self.variables['molecular_volume_cub_m']
        self.variables['volume_proportion'] = molecules_volume / self.parameters['solution_cub_m']
            
        # print calculated content
        if self.verbose:
            messages = [
                    'The photosensitizer dimensions as a {} = {} m x {} m x {} m.'.format(photosensitizer_shape, photosensitizer_length, photosensitizer_width, photosensitizer_depth), 
                    'Photosensitizer volume = {} m\N{superscript three}'.format(sigfigs_conversion(self.variables['molecular_volume_cub_m'])), 
                    'The volume proportion of {} photosensitizers = ({} m\N{superscript three} of photosensitizer)/({} m\N{superscript three} of solution) = {}'.format(
                            sigfigs_conversion(self.variables['photosensitizers']), 
                            sigfigs_conversion(molecules_volume), 
                            sigfigs_conversion(self.parameters['solution_cub_m']), 
                            sigfigs_conversion(self.variables['volume_proportion'])
                            )
                    ]
            
            self.messages.extend(messages)
            for message in messages:            
                print(message)
                
        # completion of the function
        self.defined_model.update({'define_photosensitizer':True})
                

    def _define_light(self, 
                     light_source: str = None,             # specifies a light source from the predefined options
                     light_characteristics: dict = {},      # specifies custom characteristics of the light source, in addition to or substitute of a predefined option
                     measurement: dict = None,              # provides the measurement, in the proper respective units, for the photonic intensity of the light source
                     ):
        # define properties of the light source
        self.parameters['visible_m'] = [390*nano, 780*nano]
        self.parameters['light_source'] = 'LED'
        self.light = self.light_parameters[self.parameters['light_source']]
        if light_source is not None:
            self.parameters['light_source'] = light_source
            if light_source in list(self.light_parameters.keys()): 
                self.light = self.light_parameters[self.parameters['light_source']]
#        else:
#            error = f'--> ERROR: The light source {light_source} is neither a predefined option nor a customized dictionary.'
#            self.messages.append(error)
#            ValueError(error)
            
        if light_characteristics != {}:
            for key, value in light_characteristics.items():                 
                self.light[key] = value
    
        # define the available light
        lumens_per_watt = self.light['lumens_per_watt']['value']
        if 'irradiance' in measurement: # mW / cm^2
            self.parameters['watts'] = (measurement['irradiance']*milli/centi**2) * self.area
        elif 'exposure' not in measurement:  # J / cm^2
            simulation_time = self.parameters['total_time_s'] * minute
            self.parameters['watts'] = (measurement['exposure']/centi**2) / simulation_time * self.area
        elif 'lux' not in measurement: # lumen / m^2
            self.parameters['watts'] = measurement['lux'] / lumens_per_watt * self.area
        elif 'lumens' not in measurement: # lumen
            self.parameters['watts'] = measurement['lumens'] / lumens_per_watt
        else:
            error = '--> ERROR: The light source has an unrecognized dimension.'
            self.messages.append(error)
            raise ValueError(error)
            
        # completion of the function
        self.defined_model.update({'define_light':True})
        
    def define_conditions(self,
                         bacterial_specie: str = None,    # specifies one of the bacteria in the parameters directory of PDI to simulate
                         bacterial_characteristics: dict = {},      # passes a custom dictionary of characteristics of the simulated bacterium, which can refine characteristics from the bacterial_specie argument
                         bacterial_cfu_ml: float = 1E6,   # specifies the solution concentration of the simulated bacterium for solution simulations
                         biofilm: bool = False,          # specifies whether a biofilm simulation will be conducted
                         photosensitizer: str = 'A3B_4Zn',     # specify which photosensitizer from the predefined options will be simulated
                         photosensitizer_characteristics: dict = {},   # Custom specifications of the simulation photosensitizer, which can be used to refine the parameterized photosensitizer
                         photosensitizer_molar: float = None,    # specify the photosensitizer molar concentration for solution simulations
                         photosensitizer_g: float = 90e-9,        # specify the mass of photosensitizer for cross-linked simulations
                         cross_linked_sqr_m: float = 0.0191134,  # define the cross-linked square-meters for surface simulations
                         light_source: str = None,             # specifies a light source from the predefined options
                         light_characteristics: dict = {},      # specifies custom characteristics of the light source, in addition to or substitute of a predefined option
                         measurement: dict = None,              # provides the measurement, in the proper respective units, for the photonic intensity of the light source
                          ):
        self._define_bacterium(bacterial_specie, bacterial_characteristics, bacterial_cfu_ml, biofilm)
        self._define_photosensitizer(photosensitizer, photosensitizer_characteristics, photosensitizer_molar, photosensitizer_g)                         
        self._define_light(light_source, light_characteristics, measurement)
            
    def _singlet_oxygen_calculations(self,):
        # ensure that all prior functions have completed within this instance
        missing_functions = []
        for function in self.defined_model:
            if not self.defined_model[function]:
                missing_functions.append(function)
            if missing_functions != []:
                error = f'--> ERROR: The {function} function must be defined before the simulation can execute.'
                self.messages.append(error)
                raise SystemError(error)
                
        # define the light watts
        excitation_visible_proportion = (self.parameters['excitation_range_m'] / diff(self.parameters['visible_m']))[0]
        visible_light_watts = self.parameters['watts'] * self.light['visible_proportion']['value']
        effective_excitation_watts = excitation_visible_proportion * visible_light_watts  # homogeneous light intesity throughout the visible spectrum is assumed

        # calculate the quantity of photons from the defined system
        relative_soret_excitation = 10
        weighted_average_excitation_wavelength = (self.parameters['q_m']['upper'] + self.parameters['soret_m']['lower']*relative_soret_excitation) / (1+relative_soret_excitation)
        joules_per_photon = (h*c) / weighted_average_excitation_wavelength
        non_reflected_photons = 0.96                   # “SINGLET OXYGEN GENERATION BY PORPHYRINS AND THE PHOTOSENSITIZATION IN LIPOSOMES KINETICS OF 9,lO-DIMETHYLANTHRACENE” by Gross et al., 1993
        self.variables['photon_moles_per_timestep'] = (effective_excitation_watts/joules_per_photon)*non_reflected_photons/N_A * self.parameters['timestep_s']

        # singlet oxygen calculations
        '''ambient_so = photons_per_second * molecules_dissolved_oxygen * excitation_constant'''
        self.variables['dissolved_mo_molar'] = 9*milli*self.chem_mw.mass('O2')    # mg/L -> M, ambient water quality criteria for Dissolved Oxygen, EPA   

        if self.verbose:
            messages = [
                    'photons per timestep: ', self.variables['photon_moles_per_timestep'], 
                    'molecular oxygen molecules: ', sigfigs_conversion(self.variables['dissolved_mo_molar']), 
                    'effective excitation watts: ', sigfigs_conversion(effective_excitation_watts)
                    ]
            
            self.messages.extend(messages)
            for message in messages:            
                print(message)
        
        # assign eps oxidation and attenuated bacteriual oxidation for biofilm simulations.
        if self.parameters['biofilm']:
            self.variables['eps_oxidation'] = self.bacterium['eps_oxidation_rate_constant']['value']         # for reference, 81% eps reduction versus 3-6-log bacterial reduction “Photodynamic Inactivation of Bacterial and Yeast Biofilms With a Cationic Porphyrin” by Beirao et al., 2014
            self.variables['bacterial_biofilm_mass_proportion'] = self.bacterium['cellular_dry_mass_proportion_biofilm']['value']
        
    def _kinetic_calculation(self,):      
        # ========== define photosensitizer excitation ===========
        ps = self.variables['photosensitizer_molar']
        e_fraction = min(
                1, ((self.variables['photon_moles_per_timestep']*self.variables['volume_proportion']) / (self.variables['photosensitizers']/N_A))
                )
        self.variables['e_ps_decay_time_s'] = average(1.5,15)*nano # “Ultrafast excitation transfer and relaxation inlinear and crossed-linear arrays of porphyrins” by Akimoto et al., 1999  ; "Kinetics and efficiency of excitation energy transfer from chlorophylls, their heavy metal-substituted derivatives, and pheophytins to singlet oxygen" by Küpper et al., 2002
        k_ps_rlx = 1/self.variables['e_ps_decay_time_s']     
        self.variables['ps_excitation_s'] = 50*femto # an estimated time that is below the detection limit of femto-second laser spectrophotometers; specific measurements have not been discovered.
        k_e_ps = 1/self.variables['ps_excitation_s']
        qy_e = self.photosensitizer['e_quantum_yield']['value']
        
        photosensitizer = f'{k_e_ps}*{e_fraction}*{qy_e}*ps - {k_ps_rlx}*e_ps'
            
        # ============== define photobleaching ==============
        k_b_ps = self.photosensitizer['photobleaching_constant (cm^2/J)']['value'] * (self.parameters['watts']/(self.area/centi**2))
        self.variables['hv_photobleaching_s'] = 1/k_b_ps
        
        # ============== define singlet oxygen generation ==============
        mo = self.variables['dissolved_mo_molar']
        qy = self.variables['so_qy'] = self.photosensitizer['e_quantum_yield']['value'] * self.photosensitizer['so_specificity']['value'] 
        
        self.variables['e_ps_charge_transfer_s'] = average(20,)*nano   # “Kinetics and efficiency of excitation energy transfer from chlorophylls, their heavy metal-substituted derivatives, and pheophytins to singlet oxygen” by Küpper et al., 2002 
        k_so = 1/self.variables['e_ps_charge_transfer_s']
        
        lifetime_logcfu_slope = (40-10)/(8-4)               # “The role of singlet oxygen and oxygen concentration in photodynamic inactivation of bacteria” by Maisch et al., 2007
        self.variables['so_decay_time_s'] = max(
                4, lifetime_logcfu_slope*log10(self.parameters['bacterial_cfu_ml'])
                )*micro  # a minimum lifetime of 4 microseconds is parameterized for dilute aqueous solution, according to the concensus from literature: e.g. “The role of singlet oxygen and oxygen concentration in photodynamic inactivation of bacteria” by Maisch et al., 2007 
        k_rlx_so = 1/self.variables['so_decay_time_s']
        
        # ============== define fatty acid and EPS oxidation ==============
        k_fa = self.variables['k_fa'] = 7.69E2 * self.variables['fa_gL_conc'] # https://www.jstage.jst.go.jp/article/jos/68/1/68_ess18179/_pdf/-char/ja
        fa = self.variables['fa_molar']
        k_fa_reduction = (self.parameters['bacterial_cfu_ml']/1E6)**0.1      # arbitrary reduction, dependent upon the CFU concentration
        k_fa /= k_fa_reduction
        
        biofilm = eps = ''
        if self.parameters['biofilm']:
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
        total_points = int(self.parameters['total_time_s'] / self.parameters['timestep_s'])
        self.phrasedml_str = '''
          model1 = model "pdipy_oxidation"
          sim1 = simulate uniform(0, {}, {})
          task1 = run sim1 on model1
          plot "Cytoplasmic oxidation of {}" time vs oxidation
        '''.format(self.parameters['total_time_s'], total_points, self.parameters['bacterial_specie'])
        
        # ============== execute the model in Tellurium ==============
        tellurium_model = tellurium.loada(self.model)    
        result = tellurium_model.simulate(0, self.parameters['total_time_s'], total_points)
        
        # ============== process the data ==============
        self.raw_data = pandas.DataFrame(result)
        self.raw_data.index = self.raw_data[0]
        del self.raw_data[0]
        self.raw_data.index.name = 'Time (s)'
        
        if self.parameters['biofilm']:
            self.raw_data.columns = ['[ps]', '[e_ps]', '[b_ps]', '[mo]', '[so]', '[fa]', '[ofa]', '[oeps]']
        else:
            self.raw_data.columns = ['[ps]', '[e_ps]', '[b_ps]', '[mo]', '[so]', '[fa]', '[ofa]']
        
        if self.verbose:
            messages = [
                    tellurium_model.getCurrentAntimony(), 
                    f'\nCurrent integrator:\n{tellurium_model.integrator}', 
                    f'k_fa reduction factor:{k_fa_reduction}', ]
            self.messages.extend(messages)
            for message in messages:            
                print(message)
            
            if self.jupyter:
                display(self.raw_data)
            else:
                print(self.raw_data)
                
    def _export(self, export_name, export_directory):       
        # define the simulation_path
        if export_directory is None:
            export_directory = os.getcwd()
        elif not os.path.exists(export_directory):
            error = '--> ERROR: The provided directory does not exist'
            self.messages.append(error)
            raise ImportError(error)
            
        if export_name is None:
            export_name = '-'.join([re.sub(' ', '_', str(x)) for x in ['PDIpy', self.parameters['photosensitizer_selection'], self.parameters['bacterial_specie']]])
            
        # make the simulation directory
        count = -1
        self.paths['export_path'] = os.path.join(export_directory, export_name)
        while os.path.exists(self.paths['export_path']):
            count += 1
            if not re.search('(-[0-9]+$)', self.paths['export_path']):
                self.paths['export_path'] += f'-{count}'   
            else:
                self.paths['export_path'] = re.sub('([0-9]+)$', str(count), self.paths['export_path'])
            
        os.mkdir(self.paths['export_path'])
            
        # create and export the OMEX file
        inline_omex = '\n'.join([self.model, self.phrasedml_str])  
        omex_file_name = 'pdipy_input.omex'
        tellurium.exportInlineOmex(inline_omex, os.path.join(self.paths['export_path'], omex_file_name))
        
        # export the CSVs and regression content
        self.raw_data.to_csv(os.path.join(self.paths['export_path'], 'raw_data.csv'))
        self.processed_data.to_csv(os.path.join(self.paths['export_path'], 'processed_data.csv'))
        if not self.parameters['biofilm']:                                   #!!! Why does the specification of a biofilm simulation significantly influence the ability to export the regression plot?
            self.hf.export(self.paths['export_path'], 'hillfit-regression')
        
        # export the figure
        self.figure.savefig(os.path.join(self.paths['export_path'], 'inactivation.svg'))
        if self.verbose:
            if not self.jupyter:
                self.figure.show()
        
        # define and export tables of parameters and variables
        self.parameters['simulation_path'] = self.variables['simulation_path'] = self.paths['export_path']
        parameters = {'parameter':[], 'value':[]}
        for parameter, value in self.parameters.items():
            parameters['parameter'].append(parameter)
            if isnumber(value):
                parameters['value'].append(sigfigs_conversion(value, 5))
            else:
                parameters['value'].append(value)      
                
        variables = {'variable':[], 'value':[]}
        for variable, value in self.variables.items():
            variables['variable'].append(variable)
            if isnumber(value):
                variables['value'].append(sigfigs_conversion(value, 5))
            else:
                variables['value'].append(value)
                
        parameters_table = pandas.DataFrame(parameters)
        parameters_table.to_csv(os.path.join(self.paths['export_path'], 'parameters.csv'))
        
        variables_table = pandas.DataFrame(variables)   
        variables_table.to_csv(os.path.join(self.paths['export_path'], 'variables.csv'))
        
        if self.verbose:
            if self.jupyter:
                display(parameters_table)
                display(variables_table)
            else:
                print(variables_table)     
                print(parameters_table)
    
    def simulate(self,
                 export_name: str = None,
                 export_directory: str = None,
                 figure_title: str = None,             # the figure title
                 y_label: str = 'log10',               # the label of the y-axis for the figure 
                 exposure_axis: bool = False,          # signifying exposure on the x-axis instead of time
                 display_fa_oxidation: bool = False,   # optionally overlaying the fatty acid oxidation proportion in the figure  
                 display_ps_excitation: bool = False,
                 export_contents: bool = True
                 ):
        def asymptote(xs, limit, inactivation_y_values, top_increment, top_increment_change, count, relative_to_1 = 'greater'):
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
        self._kinetic_calculation()
                    
        # parse the simulation results
        x_values = []
        oxidation_y_values = []
        excitation_y_values = []
        first = True
        for index, point in self.raw_data.iterrows():      # cytoplasmic oxidation may not be 1:1 correlated with reduction
            if first:                                      # The inital point of 0,0 crashes the regression of HillFit, and thus it is skipped
                first = False
                continue
            x_value = index/hour
            if exposure_axis:  # J/cm^2
                x_value *= self.parameters['watts']*hour/(self.area/centi**2)
            x_values.append(x_value)
            # calculate the oxidation_proportion of the cytoplasmic fatty acids
            oxidation_proportion = point['[ofa]'] / (point['[ofa]']+point['[fa]']) 
            oxidation_y_values.append(oxidation_proportion)
            
            # calculate the PS excitation proportion 
            excitation_proportion = point['[e_ps]'] / (point['[ps]'] + point['[e_ps]'] + point['[b_ps]'])
            excitation_y_values.append(excitation_proportion)

        # create the DataFrame of processed data  
        xs = array(x_values)   
        index_label = 'time (hr)'
        if exposure_axis:
            index_label = 'exposure (J/cm\N{superscript two})'
            
        self.processed_data = pandas.DataFrame(index = list(xs))
        self.processed_data.index.name = index_label
        self.processed_data['oxidation'] = oxidation_y_values
        self.processed_data['excitation'] = excitation_y_values
        self.processed_data['log10-oxidation'] = -log10(1-array(oxidation_y_values))
        self.processed_data['log10-excitation'] = -log10(1-array(excitation_y_values))
                        
        # determine the regression equation for the fitted curve via the Hill equation
        ys = array(oxidation_y_values)
        if ys[0] >= ys[-1]:
            raise ValueError(f'The last oxidation proportion {ys[-1]} is less than or equal to the first oxidation proportion {ys[0]}, which is non-physical. Change the simulation conditions and attempt another simulation.')  
        self.hf = HillFit(xs, ys)
        self.hf.fitting(x_label = index_label, y_label = 'oxidation proportion', view_figure = False)
        
        # define and refine the fitted Hill equation parameters
        increments = [0.01*self.hf.top, -0.00001*self.hf.top, 0.000001*self.hf.top, -0.0000001*self.hf.top, 0.00000001*self.hf.top]
        top_increment = 0.01*self.hf.top
        ec50_change = -.76*self.hf.ec50
        nH_change = 1.24*self.hf.nH
        limit = 1
        if self.parameters['biofilm']:
            increments = [0.01*self.hf.top, -0.00001*self.hf.top, 0.000001*self.hf.top]
            top_increment = 0.01*self.hf.top
            ec50_change = -.8*self.hf.ec50
            nH_change = -0.1*self.hf.nH
            limit = 1-10**-7

        # refine the regression equation of oxidation into plots of log-reduction 
        inactivation_y_values = list(eval(f'{self.hf.bottom} + ({self.hf.top}-{top_increment}-{self.hf.bottom})*xs**({self.hf.nH}+{nH_change}) / (({self.hf.ec50}+{ec50_change})**({self.hf.nH}+{nH_change}) + xs**({self.hf.nH}+{nH_change}))'))
        count = 0
        relative_to_1 = 'greater'
        for top_increment_change in increments:
            count, inactivation_y_values, top_increment = asymptote(xs, limit, inactivation_y_values, top_increment, top_increment_change, count, relative_to_1)
            print('refinement loop: ', count)
            if relative_to_1 == 'greater':
                relative_to_1 = 'lesser'
            else:
                relative_to_1 = 'greater'
                
        self.processed_data['inactivation'] = inactivation_y_values
        self.processed_data['log10-inactivation'] = -log10(1-array(inactivation_y_values)) 
        
        # create the figure
        pyplot.rcParams['figure.figsize'] = (11, 7)
        pyplot.rcParams['figure.dpi'] = 150
        
        self.figure, self.ax = pyplot.subplots()
        self.ax.plot(xs, self.processed_data['log10-inactivation'], label = 'Inactivation')
        if display_fa_oxidation:
            self.ax.plot(xs, self.processed_data['log10-oxidation'], label = 'Oxidation')
        if display_ps_excitation:
            sec_ax = self.ax.twinx()
            sec_ax.plot(xs, self.processed_data['excitation'], label = 'Excitation', color = 'g')
            sec_ax.set_ylabel('Photosensitizer excitation proportion', color = 'g')
            sec_ax.set_ylim(0.8,1.05)
            sec_ax.legend()
            
        self.ax.set_ylabel(y_label)
        self.ax.set_xlabel(index_label)
        if figure_title is None:
            figure_title = 'Cytoplasmic oxidation and inactivation of {} via PDI'.format(self.parameters['bacterial_specie'])
        self.ax.set_title(figure_title)
        self.ax.legend()    

        if self.verbose:
            message = f'The oxidation data was refined into inactivation data after {count} loops'
            print(message)
            self.messages.append(message)
        
        if export_contents:
            self._export(export_name, export_directory)
        
    def parse_data(self,
                     log_reduction: float = None, # the specified log_reduction that is achieved at the investigated time
                     target_hours: float = None    # the specified time at which an investigated log_reduction is achieved
                     ):
        value = unit = None
        if isnumber(log_reduction):
            if log_reduction > self.processed_data['log10-inactivation'].iloc[-1]:
                message = 'The inquired log-reduction is never reached in the simulation.'
                raise ValueError(message)
            else:
                for index, point in self.processed_data.iterrows():
                    if point['log10-inactivation'] >= log_reduction:
                        value = index
                        unit = 'hours'
                        message = f'{unit} to target: {value}'
                        break
        elif isnumber(target_hours):
            if target_hours > self.processed_data.index[-1]:
                message = '--> ERROR: The inquired time is never reached.'
                raise ValueError(message)
            else:
                for index, point in self.processed_data.iterrows():
                    if index >= target_hours:
                        value = point['log10-inactivation']
                        unit = 'log10-inactivation'
                        message = '{} at {} hours: {}'.format(unit, target_hours, value)
                        break
        else:
            message = '--> ERROR: Neither the target_time nor the target_reduction are parameterized as numbers.'
            raise TypeError(message)
            
        print(message)
        self.messages.append(message)    
            
        return value, unit