from pprint import pprint
from math import pi
from pdipy import PDI
from sigfig import round
from chemw import ChemMW
import shutil, os, re
import pandas
import matplotlib

# define simulation conditions
bacterial_characteristics = {
        "cellular_dry_mass_proportion_biofilm": {
                'value':0.3            
        },
        "membrane_thickness_nm": {
                'value': 8
        }
}
photosensitizer_characteristics = {
        "q_nm": {
                'value':[ 530, 580 ]            
        },
        "photobleaching_constant (cm^2/J)": {
                'value': 2e-16
        }
}
light_characteristics = {
        "visible_proportion": {
                'value':0.8            
        },
        "lumens_per_watt": {
                'value': 13
        }
}

def isnumber(num):
    remainder = re.sub('[0-9.e-]', '', str(num))
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

def test_init():
    # define the PDI instance
    pdi = PDI(total_time = 300)
    
    # assert the presence of content
    for bol in ['surface_system', ]:
        assert type(pdi.parameters[bol]) is bool 
    for dic in [pdi.light_parameters, pdi.photosensitizers, pdi.defined_model]:
        assert type(dic) is dict
    for quant in ['so_diffusion_m', 'total_time_s', 'timestep_s', 'solution_depth_m', 'solution_sqr_m', 'solution_cub_m']:
        assert isnumber(pdi.parameters[quant])
        
    assert type(pdi.bacteria) is list
    assert os.path.exists(pdi.parameters['root_path'])
    
def test_define_conditions():
    # define the PDI conditions
    pdi = PDI(total_time = 300)
    pdi.define_conditions('S_aureus', bacterial_characteristics, 1e5, False, 'A3B_4Zn', photosensitizer_characteristics, 18e-9, False, False, 'LED', light_characteristics, {'irradiance':8})
    
    # assert qualities of the simulation
    
    for bol in ['biofilm', ]:
        assert type(pdi.parameters[bol]) is bool 
    for quant in [
            pdi.area, 
            pdi.parameters['watts'],
            pdi.variables['fa_gL_conc'], 
            pdi.variables['fa_molar'], 
            pdi.variables['membrane_cub_m'], 
            pdi.variables['fa_mass_proportion'],
            pdi.variables['membrane_cub_m'], 
            pdi.parameters['excitation_range_m'],
            pdi.variables['photosensitizers'],
            pdi.variables['photosensitizer_molar'],
            pdi.parameters['solution_cub_m'],
            pdi.variables['molecular_volume_cub_m'],
            pdi.variables['volume_proportion']
            ]:
        assert isnumber(quant)
        
    for dic in [
            pdi.bacterium, 
            pdi.photosensitizer, 
            pdi.light,
            pdi.parameters['soret_m'],
            pdi.parameters['q_m'],
            ]:
        assert type(dic) is dict

    for param in pdi.light:
        if param in light_characteristics:
            assert pdi.light[param]['value'] == light_characteristics[param]['value']
        
    pprint(pdi.bacterium)
    for characteristic in pdi.bacterium:
        if characteristic == 'membrane_chemicals':
            for chem in pdi.bacterium[characteristic]:
                assert isnumber(pdi.bacterium[characteristic][chem]['density_gL']['value'])
                assert isnumber(pdi.bacterium[characteristic][chem]['proportion']['value'])
                for form in pdi.bacterium[characteristic][chem]['formula']:
                    assert isnumber(ChemMW().mass(form))
        if quant in bacterial_characteristics:
            assert pdi.bacterium[characteristic]['value'] == bacterial_characteristics[characteristic]['value']

    for param in pdi.photosensitizer:
        if param in photosensitizer_characteristics:
            assert pdi.photosensitizer[param]['value'] == photosensitizer_characteristics[param]['value']
            
    for string in [
            pdi.parameters['photosensitizer_selection'],
            ]:
        assert type(string) is str
    
    assert type(pdi.parameters['visible_m']) is list
    assert pdi.parameters['light_source'] == 'LED'
        
def test_simulate():            
    # execute the simulation
    pdi = PDI(total_time = 300)
    pdi.define_conditions('S_aureus', bacterial_characteristics, 1e5, False, 'A3B_4Zn', photosensitizer_characteristics, 18e-9, False, False, 'LED', light_characteristics, {'irradiance':8})
    pdi.simulate('test-PDIpy', exposure_axis = True, display_fa_oxidation = True, display_ps_excitation = True)
    
    # assert qualities of the simulation
    assert type(pdi.figure) is matplotlib.figure.Figure
    
    assert set(pdi.processed_data.columns) == set(['oxidation', 'inactivation', 'excitation', 'log10-oxidation', 'log10-inactivation', 'log10-excitation'])
    assert pdi.processed_data.index.name == 'exposure (J/cm\N{superscript two})'
    for df in [pdi.raw_data, pdi.processed_data]:
        assert type(df) is pandas.core.frame.DataFrame
        
    for path in [pdi.paths['export_path']]:
        assert os.path.exists(path)
        for file in [
                'inactivation.svg',
                'parameters.csv', 
                'pdipy_input.omex',
                'processed-data.csv',
                'raw_data.csv',
                'variables.csv'
                ]:
            assert os.path.exists(os.path.join(pdi.paths['export_path'], path))
        assert os.path.exists(os.path.join(pdi.paths['export_path'], 'hillfit-regression'))
        
    shutil.rmtree(pdi.paths['export_path'])
        
def test_parse_data():
    # execute the simulation
    pdi = PDI(total_time = 300)
    pdi.define_conditions('S_aureus', bacterial_characteristics, 1e5, False, 'A3B_4Zn', photosensitizer_characteristics, 18e-9, False, False, 'LED', light_characteristics, {'irradiance':8})
    pdi.simulate('test-PDIpy', exposure_axis = True, display_fa_oxidation = True, display_ps_excitation = True, export_contents = False)
    
    # assert that the parsed values are expected
    value, unit = pdi.parse_data(log_reduction = 4)
    assert round(value, 4) == 2.909
    assert unit == 'hours'
    
    value, unit = pdi.parse_data(target_hours = 6)
    assert round(value, 4) == 6.773
    assert unit == 'log10-inactivation'