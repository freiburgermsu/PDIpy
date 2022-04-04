from scipy.constants import femto, pico, angstrom, nano, micro, milli, centi, liter, N_A, h, c, minute, hour
from numpy import array, log10, diff
from matplotlib import pyplot
#from hillfit import HillFit
from math import pi, exp
from chemw import ChemMW
from glob import glob
import tellurium
import pandas
import sigfig
import warnings, json, re, os

pandas.set_option('max_colwidth', None)

# define useful functions
def sigfigs_conversion(num, sigfigs_in = 2):
    return sigfig.round(num, sigfigs=sigfigs_in, notation = 'sci', warn = False)

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
                summation += float(num)
                total += 1
        if total > 0:
            return summation/total
        raise ValueError(f'The arguments {num_1} and {num_2} must be numbers or a list of numbers')
    elif isnumber(num_2):
        return num_2
    else:
        raise ValueError(f'The arguments {num_1} and {num_2} must be numbers or a list of numbers')
    
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
                 solution_dimensions: dict = {},   # defines dimensions of the simulated solution 
                 well_count: str = '24',           # The well count size of the simulated petri dish, if a petri dish will be simulated
                 verbose: bool = False,  
                 printing: bool = True,
                 jupyter: bool = False
                 ):
        # define base organizations
        self.parameters, self.variables, self.results, self.paths, self.defined_model = {}, {}, {}, {}, {}
        self.messages = []
        self.verbose = verbose
        self.jupyter = jupyter
        self.printing = printing
        self.chem_mw = ChemMW(printing = False)
        
        # initial parameters
        self.parameters['so_diffusion_m'] = 80 * nano       # Moan1990 
        self.parameters['timestep_s'] = 3 * minute
#         self.parameters['oxidation_angle'] = 5           # degrees
        self.parameters['root_path'] = os.path.dirname(__file__)
        
        # identify predefined parameter options
        self.bacteria = [re.search(r'(?<=bacteria\\)(.+)(?=.json)', x).group() for x in glob(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'*.json'))]
        with open(os.path.join(self.parameters['root_path'], 'parameters', 'light.json'), encoding = 'utf-8') as light:
            self.light_parameters = json.load(light)
        with open(os.path.join(self.parameters['root_path'], 'parameters', 'photosensitizers.json'), encoding = 'utf-8') as ps:
            self.photosensitizers = json.load(ps)        
                
        # define the solution 
        with open(os.path.join(self.parameters['root_path'], 'parameters', 'wells.json')) as solution:
            self.solution = json.load(solution)[well_count]
        operational_depth = 0.7  # proportion of well volume that is filled with solution
        self.solution['depth_cm'] = self.solution['depth_cm']*operational_depth
        well_count = str(well_count)
        if solution_dimensions != {}:
            for key, value in solution_dimensions.items():
                self.solution[key] = value
                
        self.parameters['solution_depth_m'] = self.solution['depth_cm']*centi
        self.parameters['solution_sqr_m'] = self.solution['area_sqr_cm']*centi**2
        self.parameters['solution_cub_m'] = self.parameters['solution_depth_m']*self.parameters['solution_sqr_m']
        
        # initiate a tracker of completing functions that define parameters
        self.defined_model = {
                'define_system':True,
                'define_bacterium': False,
                'define_photosensitizer':False,
                'define_light':False
                }
        
                
    def _define_bacterium(self, 
                         bacterial_specie: str = 'S_aureus',  # specifies one of the bacteria in the parameters directory of PDI to simulate
                         bacterial_characteristics: dict = {},      # passes a custom dictionary of characteristics of the simulated bacterium, which can refine characteristics from the bacterial_specie argument
                         bacterial_cfu_ml: float = 1E6,       # specifies the solution concentration of the simulated bacterium for solution simulations
                         biofilm: bool = False,               # specifies whether a biofilm simulation will be conducted
                         ):    
        def shell_volume(outer_radi, inner_radi_2):
            return (4*pi/3)*(outer_radi**3 - inner_radi_2**3)  
        def shell_radius(volume):
            return pow((volume*3)/(4*pi), 1/3)
    
        self.parameters['bacterial_cfu_ml'] = bacterial_cfu_ml
        self.parameters['biofilm'] = biofilm
        
        # load the bacterial parameters
        self.parameters['bacterial_specie'] = bacterial_specie
        with open(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', f'{self.parameters["bacterial_specie"]}.json'),encoding='utf-8') as bacterium:
            self.bacterium = json.load(bacterium)
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
            total_proportion += self.bacterium['membrane_chemicals'][chemical]['proportion']['value']
        for chemical in self.bacterium['membrane_chemicals']:
            fa_density_proportion = self.bacterium['membrane_chemicals'][chemical]['density_gL']['value'] * self.bacterium['membrane_chemicals'][chemical]['proportion']['value']
            self.variables['fa_gL_conc'] += fa_density_proportion
            mw = average([self.chem_mw.mass(formula) for formula in self.bacterium['membrane_chemicals'][chemical]['formula']])
            weighed_mws.append(mw*(fa_density_proportion/total_proportion))
        self.variables['fa_molar'] = self.variables['fa_gL_conc'] / average(weighed_mws)
        
        # calculate the mass proportion of membrane fatty acids in the cell
        if self.parameters['biofilm']:
            self.variables['cell_radius_m'] = shell_radius(self.bacterium['cell_volume_fL']['value']*femto/liter)
            membrane_inner_radius = self.variables['cell_radius_m'] - self.bacterium['membrane_thickness_nm']['value']*nano
            self.variables['membrane_cub_m'] = shell_volume(self.variables['cell_radius_m'], membrane_inner_radius)
            
            fa_g = self.variables['fa_gL_conc']*self.variables['membrane_cub_m']/liter
            self.variables['fa_mass_proportion'] = fa_g/(self.bacterium['cell_mass_pg']['value']*pico)

        # completion of the function
        self.defined_model.update({'define_bacterium':True})
                

    def _define_photosensitizer(self, 
                               photosensitizer: str = 'A3B_4Zn',      # specify which photosensitizer from the predefined options will be simulated
                               photosensitizer_characteristics: dict = {},   # Custom specifications of the simulation photosensitizer, which can be used to refine the parameterized photosensitizer
                               absorbance_nm: dict = {},              # specify the absorbance of each wavelength for the simulated PS system
                               transmittance: dict = {},              # specify the transmittance of each wavelength for the simulated PS system
                               photosensitizer_molar: float = None,   # specify the photosensitizer molar concentration for solution simulations
                               surface_system: bool = False,          # specifies whether cross-linked photosensitizer will be simulated
                               photosensitizer_g: float = 90e-9,      # specify the mass of photosensitizer for cross-linked simulations
                               cross_linked_sqr_m: float = 0.0191134, # define the cross-linked square-meters for surface simulations
                               area_coverage: float = False           # the fraction of the cross-linked surface that is covered with the photosensitizer
                               ):
        # define the photosensitizer absorption properties
        self.absorbance_nm = {}
        for ab in absorbance_nm:
            if '-' in ab:
                freq_range = ab.split('-')
                for freq in range(int(freq_range[0]), int(freq_range[1])):
                    self.absorbance_nm[freq] = absorbance_nm[ab]
            else:
                self.absorbance_nm[ab] = absorbance_nm[ab]
        
        self.transmittance = {}
        for tr in transmittance:
            if '-' in tr:
                freq_range = tr.split('-')
                for freq in range(int(freq_range[0]), int(freq_range[1])):
                    transmittance[freq] = transmittance[tr]
            else:
                self.transmittance[tr] = transmittance[tr]

        # define the bottom area of the simulated solution
        self.area = self.parameters['solution_sqr_m']
        self.parameters['surface_system'] = surface_system
        if self.parameters['surface_system']:            
            self.area = self.parameters['cross_linked_surface_sqr_m'] = cross_linked_sqr_m
            self.parameters['photosensitizer_g_sqr_m'] = photosensitizer_g/cross_linked_sqr_m
            
        # load and refine photosensitizer parameters
        self.parameters['photosensitizer_selection'] = photosensitizer
        self.photosensitizer = self.photosensitizers[photosensitizer]             
        if photosensitizer_characteristics != {}:
            for key, value in photosensitizer_characteristics.items():
                if key == 'name':
                    self.parameters['photosensitizer_selection'] = value
                else:
                    if key not in self.photosensitizer:
                        self.photosensitizer[key] = {}
                    for key2, value2 in value.items():
                        self.photosensitizer[key][key2] = value2
                        
        # classify the photonic absorbance characteristics of the photosensitizer
        self.chem_mw.mass(self.photosensitizer['formula']['value'])
        photosensitizer_mw = self.chem_mw.raw_mw
        self.parameters['soret_m'] = {
                'upper': self.photosensitizer['excitation_nm']['value'][0][1]*nano, 
                'lower': self.photosensitizer['excitation_nm']['value'][0][0]*nano
                }
        self.parameters['q_m'] = {
                'upper': self.photosensitizer['excitation_nm']['value'][1][1]*nano, 
                'lower': self.photosensitizer['excitation_nm']['value'][1][0]*nano
                }
        self.parameters['excitation_range_m'] = self.parameters['q_m']['upper'] - self.parameters['soret_m']['lower']
        
        self.relative_soret_excitation = 9
        soret_spectra, q_spectra = [], []
        if self.absorbance_nm or self.transmittance != {}:
            # define the photosensitizer concentration
            if photosensitizer_molar is None:
                error = '--> ERROR: A molar photosensitizer concentration must be defined for systems of dissolved photosensitizers.'
                self.messages.append(error)
                raise TypeError(error)
            self.variables['photosensitizer_molar'] = photosensitizer_molar
            self.variables['photosensitizers'] = (self.variables['photosensitizer_molar']*N_A) * (self.parameters['solution_cub_m']/liter)
            
            if self.absorbance_nm != {}:
                for key, val in self.absorbance_nm.items():
                    key = int(key)*nano
                    if self.parameters['soret_m']['lower'] <= key <= self.parameters['soret_m']['upper']:
                        soret_spectra.append(val)
                    elif self.parameters['q_m']['lower'] <= key <= self.parameters['q_m']['upper']:
                        q_spectra.append(val)
                average_absorbance = (self.relative_soret_excitation/10)*average(soret_spectra)+ (10-self.relative_soret_excitation)/10*average(q_spectra)
                self.collided_photons_fraction = 1-10**(-average_absorbance)
            elif self.transmittance != {}:
                for key, val in self.transmittance.items():
                    key = int(key)*nano
                    if self.parameters['soret_m']['lower'] <= key <= self.parameters['soret_m']['upper']:
                        soret_spectra.append(val)
                    elif self.parameters['q_m']['lower'] <= key <= self.parameters['soret_m']['upper']:
                        q_spectra.append(val)
                average_absorbance = (self.relative_soret_excitation/10)*average(soret_spectra)+ (10-self.relative_soret_excitation)/10*average(q_spectra)        
                self.collided_photons_fraction = 1-average(average_absorbance)
        else:
            # calculate physical dimensions of the photosensitizer
            photosensitizer_length = self.photosensitizer['dimensions']['length_A']*angstrom
            photosensitizer_width = self.photosensitizer['dimensions']['width_A']*angstrom
            photosensitizer_depth = self.photosensitizer['dimensions']['depth_A']*angstrom
            photosensitizer_shape = self.photosensitizer['dimensions']['shape']

            # calculate the volume proportion that is constituted by the photosensitizer, in the volume where the photosensitizer is present        
            if photosensitizer_shape == 'disc':
                self.variables['molecular_volume_cub_m'] = (photosensitizer_length/2)**2*pi*photosensitizer_depth
            elif photosensitizer_shape == 'rect':
                self.variables['molecular_volume_cub_m'] = photosensitizer_width * photosensitizer_length * photosensitizer_depth
            else:
                error = f'--> ERROR: The volume formula for the < {photosensitizer_shape} > shape is not available.'
                self.messages.append(error)
                raise ValueError(error)
                
            # define PS concentrations and quantities
            if not self.parameters['surface_system']:
                if photosensitizer_molar is None:
                    error = '--> ERROR: A molar photosensitizer concentration must be defined for systems of dissolved photosensitizers.'
                    self.messages.append(error)
                    raise TypeError(error)
                self.variables['photosensitizer_molar'] = photosensitizer_molar
                self.variables['photosensitizers'] = (self.variables['photosensitizer_molar']*N_A) * (self.parameters['solution_cub_m']/liter)
            else:
                self.variables['photosensitizers'] = self.parameters['photosensitizer_g_sqr_m'] / (photosensitizer_mw/N_A) * self.parameters['cross_linked_surface_sqr_m']
                
                # determine the fraction of the bottom surface that is covered with photosensitizers
                if area_coverage:  
                    self.parameters['area_coverage'] = area_coverage
                    photosensitizers_area = self.parameters['cross_linked_surface_sqr_m'] * self.parameters['area_coverage']
                    self.parameters['photosensitizer_area'] = average(photosensitizer_width,photosensitizer_length)*average(photosensitizer_depth, max(photosensitizer_width,photosensitizer_length)/2)
                    self.variables['photosensitizers'] = photosensitizers_area/self.parameters['photosensitizer_area']
                    
                # determine the effective concentration of photosensitizers in the volume where these cross-linked photosensitizers reside
                self.parameters['solution_cub_m'] = self.parameters['cross_linked_surface_sqr_m']*photosensitizer_length
                self.variables['photosensitizer_molar'] = self.variables['photosensitizers']/N_A / (self.parameters['solution_cub_m']/liter)  
                
            molecules_volume = self.variables['photosensitizers'] * self.variables['molecular_volume_cub_m']
            self.variables['volume_proportion'] = molecules_volume / self.parameters['solution_cub_m']
            
            # print calculated content
            if self.verbose:
                messages = [
                        'The photosensitizer dimensions as a {} = {} m x {} m x {} m.'.format(
                                photosensitizer_shape, photosensitizer_length, photosensitizer_width, photosensitizer_depth), 
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
                     light_source: str = 'LED',             # specifies a light source from the predefined options
                     light_characteristics: dict = {},      # specifies custom characteristics of the light source, as a complement or supplement to content from the parameter files
                     measurement: dict = None, # provides the measurement, in the proper respective units, for the photonic intensity of the light source
                     ):
        self.parameters['visible_m'] = [390*nano, 780*nano]
        
        # define properties of the light source
        self.parameters['light_source'] = light_source
        self.light = self.light_parameters[self.parameters['light_source']]            
        if light_characteristics != {}:
            for key, value in light_characteristics.items():                 
                self.light[key] = value
    
        # define the available light
        lumens_per_watt = self.light['lumens_per_watt']['value']
        if 'irradiance' in measurement: # mW / cm^2
            self.parameters['watts'] = (measurement['irradiance']*milli/centi**2) * self.area
        elif 'lux' in measurement: # lumen / m^2
            self.parameters['watts'] = measurement['lux'] / lumens_per_watt * self.area
        elif 'lumens' in measurement: # lumen
            self.parameters['watts'] = measurement['lumens'] / lumens_per_watt
        else:
            error = '--> ERROR: The light source has an unrecognized dimension.'
            self.messages.append(error)
            raise ValueError(error)
            
        # completion of the function
        self.defined_model.update({'define_light':True})
        
    def _singlet_oxygen_calculations(self,):
        # ensure that all definition functions have completed
        missing_functions = []
        for function in self.defined_model:
            if not self.defined_model[function]:
                missing_functions.append(function)
        if missing_functions != []:
            error = f'--> ERROR: The {function} function(s) must be defined before the simulation can execute.'
            self.messages.append(error)
            raise SystemError(error)

        # define the watts of light that resides in the excitation range of the photosensitizer
        excitation_visible_proportion = (self.parameters['excitation_range_m'] / diff(self.parameters['visible_m'])[0])
        visible_light_watts = self.parameters['watts'] * self.light['visible_proportion']['value']
        self.parameters['effective_excitation_watts'] = excitation_visible_proportion * visible_light_watts  #!!! homogeneous intesity throughout the visible range is assumed, yet this should be refined with measured spectra of the emission distribution
        weighted_average_excitation_wavelength = (self.parameters['q_m']['upper'] + self.parameters['soret_m']['lower']*self.relative_soret_excitation) / (1+self.relative_soret_excitation)
        joules_per_photon = (h*c) / weighted_average_excitation_wavelength
        self.variables['photon_moles_per_timestep'] = (self.parameters['effective_excitation_watts']/joules_per_photon)/N_A * self.parameters['timestep_s']
        if self.absorbance_nm == {} and self.transmittance == {}:   
            # calculate the quantity of photons that potentially reach the photosensitizers
            non_reflected_photons = 0.96    
            average_photosensitizer_depth_m = self.parameters['solution_depth_m']
            if not self.parameters['surface_system']:
                average_photosensitizer_depth_m /= 2  # the average depth for dissolved photosensitizers is presumed to be the center of the solution
            non_scattered_photons = exp(-self.solution['extinction_coefficient (1/m)']*average_photosensitizer_depth_m)
            self.variables['photon_moles_per_timestep'] *= (non_reflected_photons*non_scattered_photons)

        if self.verbose:
            messages = [
                    'photon moles per timestep:\t{}'.format(self.variables['photon_moles_per_timestep']), 
                    'effective excitation watts:\t{}'.format(sigfigs_conversion(self.parameters['effective_excitation_watts']))
                    ]
            
            self.messages.extend(messages)
            for message in messages:            
                print(message)
        
        # assign eps oxidation and attenuated bacteriual oxidation for biofilm simulations.
        if self.parameters['biofilm']:
            self.variables['eps_oxidation'] = self.bacterium['eps_oxidation_rate_constant']['value']         # for reference, 81% eps reduction versus 3-6-log bacterial reduction “Photodynamic Inactivation of Bacterial and Yeast Biofilms With a Cationic Porphyrin” by Beirao et al., 2014
            self.variables['bacterial_biofilm_mass_proportion'] = self.bacterium['cellular_dry_mass_proportion_biofilm']['value']
        
    def define_conditions(self,
                         bacterial_specie: str = 'S_aureus',    # specifies one of the bacteria in the parameters directory of PDI to simulate
                         bacterial_characteristics: dict = {},  # passes a custom dictionary of characteristics of the simulated bacterium, which can refine characteristics from the bacterial_specie argument
                         bacterial_cfu_ml: float = 1E6,         # specifies the solution concentration of the simulated bacterium for solution simulations
                         biofilm: bool = False,                 # specifies whether a biofilm simulation will be conducted
                         photosensitizer: str = 'A3B_4Zn',      # specify which photosensitizer from the predefined options will be simulated
                         photosensitizer_characteristics: dict = {},   # Custom specifications of the simulation photosensitizer, which can be used to refine the parameterized photosensitizer
                         absorbance_nm: dict = {},              # specify the absorbance of each wavelength for the simulated PS system
                         transmittance: dict = {},              # specify the transmittance of each wavelength for the simulated PS system
                         photosensitizer_molar: float = None,   # specify the photosensitizer molar concentration for solution simulations
                         surface_system: bool = False,   # specifies whether cross-linked photosensitizer will be simulated
                         photosensitizer_g: float = 90e-9,      # specify the mass of photosensitizer for cross-linked simulations
                         cross_linked_sqr_m: float = 0.0191134, # define the cross-linked square-meters for surface simulations
                         area_coverage: float = False,          # the fraction of the cross-linked surface that is covered with the photosensitizer
                         light_source: str = 'LED',             # specifies a light source from the predefined options
                         light_characteristics: dict = {},      # specifies custom characteristics of the light source, in addition to or substitute of a predefined option
                         measurement: dict = None,              # provides the measurement, in the proper respective units, for the photonic intensity of the light source
                          ):
        self._define_bacterium(bacterial_specie, bacterial_characteristics, bacterial_cfu_ml, biofilm)
        self._define_photosensitizer(photosensitizer, photosensitizer_characteristics, absorbance_nm, transmittance, photosensitizer_molar, surface_system, photosensitizer_g, cross_linked_sqr_m, area_coverage)                         
        self._define_light(light_source, light_characteristics, measurement)
        self._singlet_oxygen_calculations()
        
       
    def _kinetic_calculation(self,):      
        # ========== define bacterial doubling ===========
        k_2x = self.bacterium['doubling_rate_constant']['value']
        if self.parameters['biofilm']:    
            k_2x /= 10       #!!! verify this assumption of cellular reproduction in sessile versus planktonic states.

        # ========== define photosensitizer excitation ===========
        if ('so_specificity' and 'e_quantum_yield') in self.photosensitizer:
            so_conversion = self.photosensitizer['so_specificity']['value'] 
            qy_e = self.photosensitizer['e_quantum_yield']['value']    
        elif 'so_quantum_yield' in self.photosensitizer:
            so_conversion = qy_e = self.photosensitizer['so_quantum_yield']**0.5
        else:
            message = 'Either the so_quantum_yield or the e_quantum_yield and the so_specificity quantum yields must be defined for the simulated PS.'
            self.messages.append(message)
            raise ValueError(message)
            
        ps = self.variables['photosensitizer_molar']
        if self.absorbance_nm == {} and self.transmittance == {}:
            self.collided_photons_fraction = self.variables['volume_proportion']
        excited_ps_fraction = min(
                1, ((self.variables['photon_moles_per_timestep']*self.collided_photons_fraction) / (self.variables['photosensitizers']/N_A))
                )
        k_ps_rlx = 1/(self.photosensitizer['ps_decay (ns)']['value']*nano)
        k_e_ps = 1/(self.photosensitizer['ps_rise (fs)']['value']*femto)
        photosensitizer = f'{k_e_ps}*{excited_ps_fraction}*{qy_e}*ps - {k_ps_rlx}*e_ps'
            
        # ============== define photobleaching ==============
        k_b_ps = self.photosensitizer['photobleaching_constant (cm2/(J*M))']['value'] * (self.parameters['watts']/(self.area/centi**2))
        self.variables['hv_photobleaching_s'] = 1/k_b_ps
        
        # ============== define singlet oxygen generation ==============
        mo = 9*milli/float(self.chem_mw.mass('O2'))   # mg/L -> M, ambient water quality criteria for Dissolved Oxygen, EPA   
        k_so = 1/(self.photosensitizer['ps_charge_transfer (ns)']['value']*nano)
        
        lifetime_logcfu_slope = (40-10)/(8-4) 
        self.variables['so_decay_time_s'] = max(
                3.5, lifetime_logcfu_slope*log10(self.parameters['bacterial_cfu_ml'])
                )*micro 
        k_rlx_so = 1/self.variables['so_decay_time_s']
        
        # ============== define fatty acid and EPS oxidation ==============
        self.variables['k_fa'] = k_fa = 240 * self.variables['fa_gL_conc']                   
        fa = self.variables['fa_molar']
        if not self.parameters['biofilm']:
            k_fa *= 13  #!!! The physical meaning of this empirical augmentation must be articulated
            k_fa *= (1E8/self.parameters['bacterial_cfu_ml'])**0.2      # empirical factor that imposes the intuitive inverse proportionality of extinction rate and bacterial CFU/mL
        
        biofilm = eps = ''
        if self.parameters['biofilm']:
            biofilm = 'so + eps => o_eps; {}*so*eps'.format(self.variables['eps_oxidation'])         
            eps = 'eps = {}'.format((fa/self.variables['fa_mass_proportion'])/self.bacterium['cellular_dry_mass_proportion_biofilm']['value'])

        # ============== SBML kinetic model ==============
        self.model = (f'''
          model pdipy_oxidation
            # kinetic expressions
            ps -> e_ps; {photosensitizer}
            e_ps + $mo => so + ps;  {so_conversion}*{k_so}*e_ps*mo
            ps + so => b_ps ; {k_b_ps}*ps*so
            so => $mo; {k_rlx_so}*so
            so + fa => o_fa; {k_fa}*so*fa
            {biofilm}
             => fa; {k_2x}*fa

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
        total_points = int(self.parameters['total_time_s'] / self.parameters['timestep_s'])  # the HillFit method required 3 minute timesteps & 720 minutes
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
            self.raw_data.columns = ['[ps]', '[e_ps]', '[so]', '[b_ps]', '[fa]', '[ofa]', '[eps]', '[oeps]']
        else:
            self.raw_data.columns = ['[ps]', '[e_ps]', '[so]', '[b_ps]', '[fa]', '[ofa]']
        
        if self.verbose:
            messages = [
                    tellurium_model.getCurrentAntimony(), 
                    f'\nCurrent integrator:\n{tellurium_model.integrator}', 
                    ]  
            self.messages.extend(messages)
            for message in messages:            
                print(message)
                           
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
        
        # export the simulation results 
        self.raw_data.to_csv(os.path.join(self.paths['export_path'], 'raw_data.csv'))
        self.processed_data.to_csv(os.path.join(self.paths['export_path'], 'processed_data.csv'))
#        self.hf.export(self.paths['export_path'], 'hillfit-regression')
        self.figure.savefig(os.path.join(self.paths['export_path'], 'inactivation.svg'))
        
        # define and export tables of the simulation parameters and variables
        self.parameters['simulation_path'] = self.variables['simulation_path'] = self.paths['export_path']
        parameters = {'parameter':[], 'value':[]}
        for parameter, value in self.parameters.items():
            parameters['parameter'].append(parameter)
            if isnumber(value):
                parameters['value'].append(sigfigs_conversion(value, 5))
            else:
                parameters['value'].append(value)      
        parameters_table = pandas.DataFrame(parameters)
        parameters_table.to_csv(os.path.join(self.paths['export_path'], 'parameters.csv'))
                
        variables = {'variable':[], 'value':[]}
        for variable, value in self.variables.items():
            variables['variable'].append(variable)
            if isnumber(value):
                variables['value'].append(sigfigs_conversion(value, 5))
            else:
                variables['value'].append(value)
        variables_table = pandas.DataFrame(variables)   
        variables_table.to_csv(os.path.join(self.paths['export_path'], 'variables.csv'))
        
        if self.printing:
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
                 total_time: float = 720,              # the total simulation time in minutes
                 timestep: float = 3,                  # the simulation timestep in minutes
                 experimental_data: dict = {
                         'x':[], 'y':[]
                         },                            # The x- and y-values of experimental data that will be plotted with the predicted simulation results
                 display_fa_oxidation: bool = False,   # optionally overlaying the fatty acid oxidation proportion in the figure  
                 display_ps_excitation: bool = False,
                 display_inactivation: bool = True,
                 export_contents: bool = True
                 ):
#        def asymptote(xs, limit, top_increment, top_increment_change, count, relative_to_limit):
#            if self.processed_data['inactivation'].tolist()[-1] < limit:
#                relative_to_limit = 'lesser'
#            if relative_to_limit == 'lesser':
#                while self.processed_data['inactivation'].tolist()[-1] < limit:
#                    self.processed_data['inactivation'] = list(eval(f'{self.hf.bottom} + ({self.hf.top}-{top_increment}-{self.hf.bottom})*xs**({self.hf.nH}+{nH_change}) / (({self.hf.ec50}+{ec50_change})**({self.hf.nH}+{nH_change}) + xs**({self.hf.nH}+{nH_change}))'))
#                    top_increment -= top_increment_change
#                    count += 1
#            elif relative_to_limit == 'greater':
#                while self.processed_data['inactivation'].tolist()[-1] > limit:
#                    self.processed_data['inactivation'] = list(eval(f'{self.hf.bottom} + ({self.hf.top}-{top_increment}-{self.hf.bottom})*xs**({self.hf.nH}+{nH_change}) / (({self.hf.ec50}+{ec50_change})**({self.hf.nH}+{nH_change}) + xs**({self.hf.nH}+{nH_change}))'))
#                    top_increment += top_increment_change
#                    count += 1
#            return count, top_increment
        
        # calculate the kinetics of the simulation 
        self.parameters['total_time_s'] = total_time*minute
        self.parameters['timestep_s'] = timestep*minute
        self._kinetic_calculation()
                    
        # parse the simulation results
        x_values, oxidation_ys, excitation_ys, inactivation_ys = [], [], [], []
        first = True
        for index, point in self.raw_data.iterrows():   
            if first: # The inital point of 0,0 crashes the regression of HillFit, and thus it is skipped
                first = False
                continue
            
            # calculate the oxidation_proportion of the cytoplasmic fatty acids
            oxidation_proportion = point['[ofa]'] / (point['[ofa]']+point['[fa]']) 
            if oxidation_proportion > 1:
                break
            if oxidation_ys == []:
                first_oxidation_proportion = oxidation_proportion
            if first_oxidation_proportion > oxidation_proportion:
                warnings.warn('The simulation conditions are insufficient to inactivate the bacteria.')
            oxidation_ys.append(oxidation_proportion)
            
            # determine the x-axis values
            x_value = index/hour

            x_values.append(x_value)
            
            # calculate the PS excitation proportion 
            excitation_proportion = point['[e_ps]'] / (point['[ps]'] + point['[e_ps]'] + point['[b_ps]'])
            excitation_ys.append(excitation_proportion)

        # create the DataFrame of processed data  
        xs = array(x_values)   
        index_label = 'time (hr)'
        if exposure_axis:  # J/cm^2
            xs *= self.parameters['watts']*hour/(self.area/centi**2)
            index_label = 'exposure (J/cm\N{superscript two})'
        self.processed_data = pandas.DataFrame(index = xs)
        self.processed_data.index.name = index_label
        
        # populate the DataFrame of data
        self.processed_data['oxidation'] = oxidation_ys
        self.processed_data['excitation'] = excitation_ys
        self.processed_data['log10-oxidation'] = -log10(1-array(oxidation_ys))
        self.processed_data['log10-excitation'] = -log10(1-array(excitation_ys))
        
        # define the figure
        pyplot.rcParams['figure.figsize'] = (11, 7)
        pyplot.rcParams['figure.dpi'] = 150
        
        self.figure, self.ax = pyplot.subplots()
        if display_fa_oxidation:
            self.ax.plot(xs, self.processed_data['log10-oxidation'], label = 'Oxidation')
        if experimental_data['x'] != []:
            self.ax.scatter(experimental_data['x'], experimental_data['y'], label = 'Experimental Inactivation')
        if display_ps_excitation:
            if not display_fa_oxidation and not display_inactivation:
                self.ax.set_ylabel('Photosensitizer excitation proportion')
                self.ax.set_xlabel(index_label)
                self.ax.plot(xs, self.processed_data['excitation'], label = 'Excitation', color = 'g')
                self.ax.set_ylim(
                    min(self.processed_data['excitation'])-.05,
                    min(1,max(self.processed_data['excitation']))+.05
                    )
            else:
                sec_ax = self.ax.twinx()
                sec_ax.plot(xs, self.processed_data['excitation'], label = 'Excitation', color = 'g')
                sec_ax.set_ylabel('Photosensitizer excitation proportion', color = 'g')
                sec_ax.set_ylim(
                    min(self.processed_data['excitation'])-.05,
                    min(1,max(self.processed_data['excitation']))+.05
                    )
                sec_ax.legend(loc = 'lower right')            
        if figure_title is None:
            figure_title = 'Cytoplasmic oxidation and inactivation of {}'.format(self.parameters['bacterial_specie'])
        self.ax.set_title(figure_title)
        
        if display_inactivation or display_fa_oxidation:
            self.ax.set_ylabel(y_label)
            self.ax.set_xlabel(index_label)
            self.ax.legend(loc = 'lower center')    

        if self.printing:
            if self.jupyter:
                display(self.processed_data)
            else:
                print(self.processed_data)   
                        
        # determine the regression equation for the fitted curve via the Hill equation                
#        ys = array(oxidation_ys)
#        if ys[0] >= ys[-1]:
#            raise ValueError(f'The rate of inactivation is less than the rate of reproduction. Change the simulation conditions and attempt another simulation.')  
#        self.hf = HillFit(xs, ys)
#        self.hf.fitting(x_label = index_label, y_label = 'oxidation proportion', view_figure = False)
        
        # define and refine the fitted Hill equation parameters
#        for y in reversed(oxidation_ys):
#            if y<1:
#                final_y = y
#                break

        # extrapolate inactivation from oxidation
        self.processed_data['log10-inactivation'] = -log10(1-array(oxidation_ys)) + -log10(
                self.bacterium['biofilm_oxidation_fraction_lysis']['value']
                )
        if not self.parameters['biofilm']:
            self.processed_data['inactivation'] = array(oxidation_ys)
            
            # Hill-equation method
#            count =1
#            # define parameter changes
#            num_increments = 8
#            increments = logspace(-1,-(num_increments+4),ceil(num_increments))*self.hf.top
#            top_increment = 0.01*self.hf.top
#            ec50_change = .5*self.hf.ec50
#            nH_change = 3.5*self.hf.nH # +self.parameters['watts']
#            limit = 1-10**-(self.parameters['watts']**(0.2)-log10(1-final_y))
#            
#            # refine the regression equation of oxidation into plots of log-reduction 
#            self.processed_data['inactivation'] = eval(f'{self.hf.bottom} + ({self.hf.top}-{top_increment}-{self.hf.bottom})*xs**({self.hf.nH}+{nH_change}) / (({self.hf.ec50}+{ec50_change})**({self.hf.nH}+{nH_change}) + xs**({self.hf.nH}+{nH_change}))')
#            count = 0
#            relative_to_limit = 'greater'
#            for top_increment_change in increments:
#                count, top_increment = asymptote(xs, limit, top_increment, top_increment_change, count, relative_to_limit)
#                if self.printing:
#                    print('refinement loop: ', count)
#                if relative_to_limit == 'greater':
#                    relative_to_limit = 'lesser'
#                else:
#                    relative_to_limit = 'greater'
            
            # geometric translation method
            extrapolation = 4.5 + log10(self.parameters['effective_excitation_watts']/0.004314663461538462)   # the denominator is an empirical calibration, which imposes the intuitive proportionality between light intensity and inactivation
            if self.parameters['surface_system']:
                extrapolation = 2
            self.processed_data['log10-inactivation'] = -log10(1-array(self.processed_data['inactivation'])) + extrapolation

        if display_inactivation:
            self.ax.plot(xs, self.processed_data['log10-inactivation'], label = 'Inactivation')
        self.ax.legend(loc = 'lower center')

        if self.verbose:
            # for the Hill-equation method
#            if not self.parameters['biofilm']:
#                message = f'The oxidation data was refined into inactivation data after {count} loops'
#                print(message)
#                self.messages.append(message)
            
            if self.jupyter:
                display(self.raw_data)
            else:
                print(self.raw_data)
                self.figure.show()
        
        if export_contents:
            self._export(export_name, export_directory)
        
    def parse_data(self,
                     log_reduction: float = None,  # the specified log-reduction that is achieved at the investigated time
                     target_hours: float = None    # the specified time at which an investigated log_reduction is achieved
                     ):
        value = unit = None
        if isnumber(log_reduction):
            if log_reduction > max(self.processed_data['log10-inactivation'].squeeze()):
                message = 'The inquired log-reduction is never reached. Change the simulation conditions.'
                raise ValueError(message)
            else:
                for index, point in self.processed_data.iterrows():
                    if point['log10-inactivation'] >= log_reduction:
                        value = index
                        unit = 'hours'
                        message = f'{unit} to target: {value}'
                        break
                print(message)
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