import pdipy


def test_cytoplasmic_membrane():
    # define calculation functions
    def shell_volume(radi_1, radi_2, coeff = 1):
        volume = (4*pi/3) * coeff * (radi_1**3 - radi_2**3)
        return volume        
    def sector_volume(r, h):
        volume = (2*pi/3)*r**2*h
        return volume       
    
    if self.bacterium['shape']['value'] == "sphere":            
        # calculate the cellular dimensions
        self.variables['cell_radius (m)'] = pow((self.bacterium['cell_volume (pL)']['value']*pico/liter*3) / (4*pi), 1/3)
        membrane_inner_radius = self.variables['cell_radius (m)'] - self.bacterium['membrane_thickness (nm)']['value'] * nano

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