from math import pow, pi
import pdipy


def test_cytoplasmic_membrane():
    # define calculation functions     
    def sector_volume(r, h):
        volume = (2*pi/3)*r**2*h
        return volume       
    
    if self.bacterium['shape']['value'] == "sphere":            
        # calculate the cellular dimensions
        

        #calculate the cellular and oxidation areas
        oxidized_cap_area = 2 * pi * self.variables['cell_radius (m)'] * outer_h
        cell_area = 4 * pi * self.variables['cell_radius (m)'] ** 2
        self.variables['oxidized_area_ratio'] = oxidized_cap_area / cell_area

        message2 = 'volume:area consistency', round(self.variables['oxidized_area_ratio'],7) == round(self.variables['oxidized_membrane_volume_ratio'],7)

        # calculate the cellular and oxidation volumes
        membrane_volume = shell_volume(self.variables['cell_radius (m)'], membrane_inner_radius) # M^3
        outer_h = self.variables['cell_radius (m)'] * (1 - cos(radians(self.parameters['oxidation_angle'])))
        inner_h = membrane_inner_radius * (1 - cos(radians(self.parameters['oxidation_angle'])))
        shell_sector_volume = sector_volume(self.variables['cell_radius (m)'], outer_h) - sector_volume(membrane_inner_radius, inner_h)
        self.variables['oxidized_membrane_volume_ratio'] = shell_sector_volume / membrane_volume

        oxidized_volume = membrane_volume * self.variables['oxidized_membrane_volume_ratio']

    if self.verbose:
        message1 = 'oxidized volume proportion: ', self.variables['oxidized_membrane_volume_ratio']
        self.messages.extend([message1, message2])

        print(message1)
        print(message2)        
        
        
def test_photosensitizer_dimensions():
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
    
    if self.parameters['photosensitizer_selection'] == 'A3B_4Zn':
        # calculate the volume of conjugated region of the photosensitizer
        self.variables['center_porphyrin_length'] = (2*(chemical_dimensions['bond']['c-c']*(2*cos(radians(chemical_dimensions['angle']['sp2']-90))+cos(radians(180-chemical_dimensions['angle']['sp2']))))) * angstrom
        self.variables['sp2_extension'] = (chemical_dimensions['bond']['c-c'] * (2 + cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-n'] * cos(radians(180-chemical_dimensions['angle']['sp2'])) + chemical_dimensions['bond']['c-n']) * angstrom
        self.variables['sp3_diazirine'] = (chemical_dimensions['bond']['c-c']*cos(radians(chemical_dimensions['angle']['sp3']-90)) + 2*chemical_dimensions['bond']['c-c']*(1+cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-c']*sin(radians(chemical_dimensions['angle']['sp3']-90)) + chemical_dimensions['bond']['c-f']) * angstrom
        self.variables['linked_sp3_diazirine'] = (chemical_dimensions['bond']['c-c']*cos(radians(chemical_dimensions['angle']['sp3']-90)) + 2*chemical_dimensions['bond']['c-c']*(1+cos(radians(180-chemical_dimensions['angle']['sp2']))) + chemical_dimensions['bond']['c-c']*sin(radians(chemical_dimensions['angle']['sp3']-90))) * angstrom
        
        conjugated_length = self.variables['center_porphyrin_length'] + 2*self.variables['sp2_extension']

        total_length = self.variables['center_porphyrin_length'] + 2*(self.variables['sp2_extension'] + self.variables['sp3_diazirine'])
        linked_total_length = self.variables['center_porphyrin_length'] + 2* self.variables['sp2_extension'] + self.variables['sp3_diazirine'] + self.variables['linked_sp3_diazirine']
        
        message1 = 'The tetrapyrrole length is {} meters'.format(sigfigs_conversion(self.variables['center_porphyrin_length']))
        message2 = 'The benzyl extension is {} meters'.format(sigfigs_conversion(self.variables['sp2_extension']))
        message3 = 'The diazirine is {} meters'.format(sigfigs_conversion(self.variables['sp3_diazirine']))
        
        

        message6 = ''
        if not self.surface_system:
            tilted_height = total_length*float(sin(radians(45)))
            layers = 1
            if not self.surface_system:
                layers = self.parameters['solution_depth (m)'] / tilted_height
            photosensitizers_per_layer = self.variables['photosensitizers'] / layers
            orthogonal_area = pi*(total_length)**2
            parallel_area = (total_length)*self.parameters['atomic_thickness']
            photosensitizers_layer_area = photosensitizers_per_layer*average(orthogonal_area, parallel_area)
            self.variables['area_proportion'] = photosensitizers_layer_area/self.parameters['surface_area (m^2)'] 
            message6 = 'The photosensitizer area proportion is {}'.format(sigfigs_conversion(self.variables['area_proportion']))