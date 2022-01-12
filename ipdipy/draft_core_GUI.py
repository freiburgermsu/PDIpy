# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 13:24:14 2021

@author: Ethan
"""

# -*- coding: utf-8 -*-
"""
@authors: Ethan Chan, Matthew Freiburger
"""

#Import Statements
import os #File Managment
from tkinter import *
from tkinter import ttk
import scipy.constants
import math
import re
from datetime import date
import os.path
from os import path

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from tkinter.filedialog import asksaveasfile
from tkinter.filedialog import askopenfile

#Normal Program
bacterial_species = ""
cell_shape = ""
cellular_vitality = ""
light_source = ""
plt.figure(figsize=(3,3))

csv_header_flag = True

oxidation_graph_data = []
oxidation_graph_timesteps = []
graphed_oxidation_molecules = []
plot_out_path = ""

log_file = open("readme.txt", "w")
csv_file = open("readme.txt", "w")

def append_to_console_box(message):
    global log_file
    log_file.write(message + "\n")
    console_box.configure(state='normal')
    console_box.insert("end", message + "\n")
    console_box.configure(state='disabled')
    
    
def run_calculations(timestep_entry, species_entry, porphyrin_entry, porphyrin_concentration_entry, light_spectrum_entry, luminous_intensity_entry, light_source_entry, light_watts_entry):
    t = threading.Thread(taget=initial_conditions, args=(timestep_entry, species_entry, porphyrin_entry, porphyrin_concentration_entry, light_spectrum_entry, luminous_intensity_entry, light_source_entry, light_watts_entry))
    t.start()
    
def import_environment(timestep_entry, species_entry, porphyrin_entry, porphyrin_concentration_entry, light_spectrum_entry, luminous_intensity_entry, light_source_entry, light_watts_entry):
    file = askopenfile(mode ='r', filetypes =[('Text Files', '*.txt')])
    if file is not None:
        content = file.read()
        if "timestep,species,photosensitizer,porphyrin_concentration,q_soret,luminous,light_source,light_wattage\n" in content:
            content = content.replace("timestep,species,photosensitizer,porphyrin_concentration,q_soret,luminous,light_source,light_wattage", "").replace("\n", "").replace("\r", "")
            fields = content.split(",")
            
            
            try:
                temp = float(fields[0])
                timestep_entry.delete(0, END)
                timestep_entry.insert(0, temp)
            except:
                pass
            
            try:
                temp = str(fields[1])
                if temp in species_choices:
                    species_entry.delete(0, END)
                    species_entry.insert(0, temp)
            except:
                pass
            
            
            try:
                temp = str(fields[2])
                if temp in photosensitizer_choices:
                    porphyrin_entry.delete(0, END)
                    porphyrin_entry.insert(0, temp)
            except:
                pass
            
            try:
                temp = float(fields[3])
                porphyrin_concentration_entry.delete(0, END)
                porphyrin_concentration_entry.insert(0, temp)
            except:
                pass
            
            try:
                temp = str(fields[4])
                if temp in excite_choices:
                    light_spectrum_entry.delete(0, END)
                    light_spectrum_entry.insert(0, temp)
            except:
                pass
            
            try:
                temp = float(fields[5])
                luminous_intensity_entry.delete(0, END)
                luminous_intensity_entry.insert(0, temp)
            except:
                pass
            
            try:
                temp = str(fields[6])
                if temp in source_choices:
                    light_source_entry.delete(0, END)
                    light_source_entry.insert(0, temp)
            except:
                pass
            
            try:
                temp = float(fields[7])
                light_watts_entry.delete(0, END)
                light_watts_entry.insert(0, temp)
            except:
                pass
            

def export_environment(timestep_entry, species_entry, porphyrin_entry, porphyrin_concentration_entry, light_spectrum_entry, luminous_intensity_entry, light_source_entry, light_watts_entry):
    f = asksaveasfile(mode='w', defaultextension=".txt")
    if f is None:
        return
    f.write("timestep,species,photosensitizer,porphyrin_concentration,q_soret,luminous,light_source,light_wattage\n")
    f.write(str(timestep_entry.get()) + "," + str(species_entry.get()) + "," + str(porphyrin_entry.get()) + "," + str(porphyrin_concentration_entry.get()) + "," + str(light_spectrum_entry.get()) + "," + str(luminous_intensity_entry.get()) + "," + str(light_source_entry.get()) + "," + str(light_watts_entry.get()))
    f.close()

    
def initial_conditions(timestep_entry, species_entry, porphyrin_entry, porphyrin_concentration_entry, light_spectrum_entry, luminous_intensity_entry, light_source_entry, light_watts_entry):
    """
        Description:
            Parameterizes and calculates the initial conditions of the PDI system

        Used:
            main()
    """
    global singlet_oxygen_concentration
    global oxidized_death_threshold
    global stop_biofilm_threshold
    global membrane_chemicals
    global cellular_vitality
    global cell_shape
    global timestep
    global bacterial_species
    global oxidation_graph_data
    global oxidation_graph_timesteps
    global graphed_oxidation_molecules
    global log_file
    global csv_file
    global plot_out_path
    
    console_box.configure(state='normal')
    console_box.delete(1.0, "end")
    console_box.configure(state='disabled')
    
    return_flag = False
    
    timestep = 0
    try:
        timestep = float(timestep_entry.get())
    except:
        console_box.configure(state='normal')
        console_box.insert("end", "INPUT ERROR: Time step should be a number.\n")
        console_box.configure(state='disabled')
        return_flag = True
    
    porphyrin_concentration_mgl = 0
    try:
        porphyrin_concentration_mgl = float(porphyrin_concentration_entry.get())
    except:
        console_box.configure(state='normal')
        console_box.insert("end", "INPUT ERROR: Porphyrin concentration should be a number.\n")
        console_box.configure(state='disabled')
        return_flag = True
    
    luminous_intensity = 0
    try:
        luminous_intensity = float(luminous_intensity_entry.get())
    except:
        console_box.configure(state='normal')
        console_box.insert("end", "INPUT ERROR: Luminous intensity should be a number.\n")
        console_box.configure(state='disabled')
        return_flag = True
    
    light_watts = 0
    try:
        light_watts = float(light_watts_entry.get())
    except:
        console_box.configure(state='normal')
        console_box.insert("end", "INPUT ERROR: Light watts should be a number.\n")
        console_box.configure(state='disabled')
        return_flag = True
    
    if return_flag:
        return
    
    today = date.today()
    directory_name = "Membrane-Simulation-" + today.strftime("%b-%d-%Y")
    
    index_file = 0
    
    while True:
        directory_path = directory_name + "-" + str(index_file)
        
        if not path.exists(directory_path):
            os.mkdir(directory_name + "-" + str(index_file))
            log_file = open(directory_path + "/console-log-" + str(index_file) + ".txt", "w")
            csv_file = open(directory_path + "/oxidation-data-" + str(index_file) + ".csv", "w")
            plot_out_path = directory_path + "/oxidation-graph-" + str(index_file) + ".png"
            break
        index_file += 1
    
    oxidation_graph_data = []
    oxidation_graph_timesteps = []
    graphed_oxidation_molecules = []


    
    




    
   
   
    # initial biological conditions
    cellular_vitality = "alive"
    biofilm_growth = 'yes'
    bacterial_species = "S aureus"
    bacterial_species = species_entry.get()
    if bacterial_species == 'S aureus':
        membrane_chemicals = {'anteiso_C17:0':{'density': 1,
                                               'chemical formula': 'C_18_H_35_O_2',  # determined from the structure
                                               'molecular weight': 283.476,  # determined from the structure
                                               'proportion': 0.2},
                              'oleic acid':{'density': 0.895, # g / mL
                                            'chemical formula': 'C_18_H_34_O_2',  # determined from the structure
                                            'molecular weight': 282.468,  # determined from the structure
                                            'proportion': 0.2},
                              'anteiso_C15:0':{'density': 1,
                                               'chemical formula': 'C_16_H_31_O_2',  # determined from the structure
                                               'molecular weight': 255.422,  # determined from the structure
                                               'proportion': 0.2}
                         }
        cell_shape = "sphere"
        stop_biofilm_threshold = 0.05
        oxidized_death_threshold = 0.1

    # system constants
    mw_molecular_oxygen = 32 / scipy.constants.milli          #mg / mole
    solution_volume = 100 * scipy.constants.micro                          #liters
    dissolved_oxygen_concentration = 9              # mg / L, ambient water quality criteria for DO, EPA
    molecules_dissolved_oxygen = dissolved_oxygen_concentration * solution_volume / mw_molecular_oxygen * scipy.constants.N_A


    # organometallic photosensitizer system
    porphyrins_and_corroles = {'A3B-tetracationic zinc porphyrin': {'quantum yield':0.6,
                                                                    'so specificity':0.8,
                                                                    'chemical formula': 'C_76_H_48_N_16_F_12_Zn',  # determined from the structure
                                                                    'molecular weight': 1478.688},  # counted from the atomic proportions of the structure
                               'A3B-tetracationic copper porphyrin': {'quantum yield':0.5,
                                                                    'so specificity': 0.9,
                                                                    'chemical formula': 'C_76_H_48_N_16_F_12_Zn',  # determined from the structure
                                                                    'molecular weight': 1478.688},  # counted from the atomic proportions of the structure
                               'A3B-dicationic zinc porphyrin':{'quantum yield':0.6,
                                                                'so specificity': 1,
                                                                'chemical formula': 'C_60_H_28_N_10_F_16_Zn',  # determined from the structure
                                                                'molecular weight': 1258.302},  # counted from the atomic proportions of the structure
                               'A3-monocationic gallium corrole':{'quantum yield':0.5,
                                                                  'so specificity':0.9,
                                                                  'chemical formula': 'C_45_H_18_N_7_F_13_Ga',  # determined from the structure
                                                                  'molecular weight': 973.385}  # counted from the atomic proportions of the image
                              }
   


           
    porphyrin_selection = porphyrin_entry.get()
    
    
    
    
    
    porphyrin_quantity_moles = porphyrin_concentration_mgl * solution_volume * scipy.constants.milli / porphyrins_and_corroles[porphyrin_selection]['molecular weight']
    if re.search('(porphyrin)', porphyrin_selection):
        soret_band_upper_bound = 430 * scipy.constants.nano
        soret_band_lower_bound = 400 * scipy.constants.nano
        q_band_upper_bound = 625 * scipy.constants.nano
        q_band_lower_bound = 530 * scipy.constants.nano
   
    light_spectrum = light_spectrum_entry.get()

    
    
   
   
    light_source = light_source_entry.get()

   
    visible_proportion = 0.1

    if light_source == 'Incandescent':   #Macisaac et al., 1999
        visible_proportion = 0.1
    elif light_source == 'LED':    #~80% more energy efficient and thus 5x more visible light per energy relative to incandescent bulbs
        visible_proportion = 0.5
       
    
           
    effective_visible_light_watts = light_watts * visible_proportion
    visible_region = 780E-9 - 390E-9
    excitation_visible_proportion = ((q_band_upper_bound - q_band_lower_bound) + (soret_band_upper_bound - soret_band_lower_bound)) / visible_region
    effective_excitation_light_watts = excitation_visible_proportion * effective_visible_light_watts  # homogeneous light intesity throughout the visible spectrum is assumed

    # photonic calculations
    average_excitation_wavelength = (q_band_upper_bound + soret_band_lower_bound) / 2
    joules_per_photon = (scipy.constants.h * scipy.constants.c) / average_excitation_wavelength
    photons_per_second = effective_excitation_light_watts / joules_per_photon

    # limiting reagent calculations
    so_from_light = photons_per_second * molecules_dissolved_oxygen
   
    so_quantum_yield = porphyrins_and_corroles[porphyrin_selection]['quantum yield'] * porphyrins_and_corroles[porphyrin_selection]['so specificity']
    molecules_photosensitizer = porphyrin_quantity_moles * scipy.constants.N_A
    so_from_photosensitizer = molecules_photosensitizer * so_quantum_yield * molecules_dissolved_oxygen
   
    if so_from_photosensitizer > so_from_light:  # the excitation collision is assumed to be in the first order of both the porphyrin and light
        limiting_component = 'light source'
        percent_difference = (so_from_photosensitizer - so_from_light) / (so_from_photosensitizer) * 100

       
        append_to_console_box('The limiting excitation agent is %s, which is %s %% less than the photosensitizer concentration' %(limiting_component, percent_difference))

        singlet_oxygen_amount = so_from_light
       
    elif so_from_photosensitizer < so_from_light:
        limiting_component = 'photosensitizer concentration'
        percent_difference = (so_from_light - so_from_photosensitizer) / so_from_light * 100
        append_to_console_box('The limiting excitation agent is the %s, which is %s %% less than the light intensity' %(limiting_component, percent_difference))
        singlet_oxygen_amount =  so_from_photosensitizer

    append_to_console_box('Dissolved oxygen molecules: ' + str(molecules_dissolved_oxygen))
    if singlet_oxygen_amount > molecules_dissolved_oxygen:
        singlet_oxygen_amount = molecules_dissolved_oxygen
       
    singlet_oxygen_concentration = singlet_oxygen_amount / solution_volume
    append_to_console_box('Generated singlet oxygen molecules: ' + str(singlet_oxygen_amount))
   
    main()

   
def cell_geometry():
    """
        Description:
            Calculates geometric properties of the cell such as surface area, absorption volume, and cell_radius

        Used:
            main()
    """
    global membrane_solution_interface_volume
    global oxidized_volume_proportion
    global cell_surface_area
    global membrane_volume
    global cell_radius
    global cell_volume
    global cell_mass
    global bacterial_species
   



    if bacterial_species == 'S aureus':
        membrane_thickness = 50 * scipy.constants.nano # ____ citation needed
        cell_mass = 1 * scipy.constants.pico # ____ citation needed
        cell_volume = 1 * scipy.constants.femto # ____ citation needed
       
    if cell_shape == "sphere":  
        cell_radius = math.pow((cell_volume * 3) / (4 * math.pi), 1/3)
        membrane_outer_radius = cell_radius
        membrane_inner_radius = cell_radius - membrane_thickness
        membrane_volume = (4 * math.pi / 3) * (membrane_outer_radius ** 3 - membrane_inner_radius ** 3) # Liters
   
        #calculate the region of membrane oxidation
        oxidation_angle = 5
        outer_component_height = cell_radius * math.cos(math.radians(oxidation_angle))
        outer_h = cell_radius - outer_component_height
        outer_oxidized_cap_area = 2 * math.pi * cell_radius * outer_h
        inner_component_height = (cell_radius - membrane_thickness) * math.cos(math.radians(oxidation_angle))
        inner_h = (cell_radius - membrane_thickness) * inner_component_height
        inner_oxidized_cap_area = 2 * math.pi * (cell_radius - membrane_thickness) * inner_h
       
        # the membrane area of oxidation  
        outer_cell_area = 4 * math.pi * membrane_outer_radius ** 2
        outer_area_ratio = outer_oxidized_cap_area / outer_cell_area
        inner_cell_area = 4 * math.pi * membrane_inner_radius ** 2
        inner_area_ratio = inner_oxidized_cap_area / inner_cell_area
        consistency = (outer_area_ratio) == (inner_area_ratio)
        append_to_console_box('inner and outer ratio equivalency: ' + str(consistency))
        if not consistency:
            digit_consistency = 3

            append_to_console_box('inner and outer ratio equivalency to ' + str(digit_consistency) + ' digits: ' + str(round(outer_area_ratio, digit_consistency) == round(inner_area_ratio, digit_consistency)))
            if round(outer_area_ratio, digit_consistency) != round(inner_area_ratio, digit_consistency):
                append_to_console_box('membrane outer area ratio: ' + str(outer_area_ratio))
                append_to_console_box('membrane inner area ratio: ' + str(inner_area_ratio) + '\n\n')
       
        cap_volume = (2 * math.pi / 3) * (1 - math.cos(math.radians(oxidation_angle))) * (membrane_outer_radius ** 3 - membrane_inner_radius ** 3)

        oxidized_volume_proportion = (cap_volume / membrane_volume)

        singlet_oxygen_diffusion_distance = 75 * scipy.constants.nano # ____ citation needed
        singlet_oxygen_interaction_radius = cell_radius + singlet_oxygen_diffusion_distance
        membrane_solution_interface_volume = (2 * math.pi / 3) * (1 - math.cos(math.radians(oxidation_angle))) * (singlet_oxygen_interaction_radius ** 3 - membrane_outer_radius ** 3)


def reactions():
    """
        Description:
            Executes the simulated reactions

        Used:
            main()
    """  
    global membrane_chemical_concentrations
    global oxidized_chemical_quantities
    global bacterial_species
    global csv_header_flag
   
    if bacterial_species == 'S aureus':
        c17_oxidized_volume = membrane_chemicals['anteiso_C17:0']['proportion'] * membrane_volume * oxidized_volume_proportion  # Liters
        anteiso_c17_0_gl = c17_oxidized_volume * membrane_chemicals['anteiso_C17:0']['density'] * scipy.constants.milli / membrane_solution_interface_volume
       
        oleic_oxidizd_volume = membrane_chemicals['oleic acid']['proportion'] * membrane_volume * oxidized_volume_proportion  # Liters
        oleic_acid_gl = oleic_oxidizd_volume * membrane_chemicals['oleic acid']['density']  * scipy.constants.milli / membrane_solution_interface_volume
       
        c15_oxidized_volume = membrane_chemicals['anteiso_C15:0']['proportion'] * membrane_volume * oxidized_volume_proportion  # Liters
        anteiso_c15_0_gl = c15_oxidized_volume * membrane_chemicals['anteiso_C15:0']['density'] * scipy.constants.milli / membrane_solution_interface_volume
               
        if time == 0:
            membrane_chemical_concentrations = {'anteiso_C17:0': anteiso_c17_0_gl  / membrane_chemicals['anteiso_C17:0']['molecular weight'],
                                                'oleic acid': oleic_acid_gl / membrane_chemicals['oleic acid']['molecular weight'],
                                                'anteiso_C15:0': anteiso_c15_0_gl / membrane_chemicals['anteiso_C15:0']['molecular weight']}
        if csv_header_flag:
            csv_header_string = "time"
            for molecule in membrane_chemical_concentrations:
                csv_header_string = csv_header_string + "," + molecule
            csv_file.write(csv_header_string)
            csv_header_flag = False
         
        # arbitrary example kinetic expressions
        '''The general form of the  kinetics reactions:
        (vmax * [fatty acid] * [SO]) / (((km1 + 1) * [fatty acid]) + ((km2 + 1) * [SO]))'''
       
        membrane_reactions = {'anteiso_C17:0':{'vmax': .000010,
                                            'km1': .2,
                                            'km2': 1},
                             'oleic acid':{'vmax': .000010,
                                            'km1': .2,
                                            'km2': 1},
                             'anteiso_C15:0':{'vmax': .000010,
                                            'km1': .2,
                                            'km2': 1}
                             }
        c17_michaelis_numerator = membrane_chemical_concentrations['anteiso_C17:0'] * singlet_oxygen_concentration * membrane_reactions['anteiso_C17:0']['vmax']
        c17_michaelis_denominator = ((membrane_reactions['anteiso_C17:0']['km1'] + 1) * membrane_chemical_concentrations['anteiso_C17:0'] + (membrane_reactions['anteiso_C17:0']['km2'] + 1) * singlet_oxygen_concentration)
       
        oleic_michaelis_numerator = membrane_chemical_concentrations['oleic acid'] * singlet_oxygen_concentration * membrane_reactions['oleic acid']['vmax']
        oleic_michaelis_denominator = ((membrane_reactions['oleic acid']['km1'] + 1) * membrane_chemical_concentrations['oleic acid'] + (membrane_reactions['oleic acid']['km2'] + 1) * singlet_oxygen_concentration)
       
        c15_michaelis_numerator = membrane_chemical_concentrations['anteiso_C15:0'] * singlet_oxygen_concentration * membrane_reactions['anteiso_C15:0']['vmax']
        c15_michaelis_denominator = ((membrane_reactions['anteiso_C15:0']['km1'] + 1) * membrane_chemical_concentrations['anteiso_C15:0'] + (membrane_reactions['anteiso_C15:0']['km2'] + 1) * singlet_oxygen_concentration)
       
        membrane_kinetics_oxidation = {'anteiso_C17:0': c17_michaelis_numerator / c17_michaelis_denominator,
                                       'oleic acid': oleic_michaelis_numerator / oleic_michaelis_denominator,
                                       'anteiso_C15:0': c15_michaelis_numerator / c15_michaelis_denominator
                                      }
       
       
        # oxidation healing kinetics
        lipid_recovery = 5000 # molecules / second
        recover_timestep_delay = 10 # the number of timesteps after which lipid reduction can occur
       
    if time == 0:
        oxidized_chemical_quantities = {}

    for molecule in membrane_kinetics_oxidation:
        if time == 0:
            oxidized_chemical_quantities['oxidized %s' %(molecule)] = 0

        else:    
            oxidized_molecule_quantity = membrane_kinetics_oxidation[molecule] * timestep * membrane_solution_interface_volume / scipy.constants.N_A
            if timestep >= recover_timestep_delay and oxidized_chemical_quantities['oxidized %s' %(molecule)] > 0:
                reduced_molecule_quantity = lipid_recovery * timestep
            else:
                reduced_molecule_quantity = 0

            net_molecule_oxidation = oxidized_molecule_quantity - reduced_molecule_quantity
            oxidized_chemical_quantities['oxidized %s' %(molecule)] += net_molecule_oxidation
            membrane_chemical_concentrations[molecule] += (net_molecule_oxidation / membrane_solution_interface_volume) / scipy.constants.N_A


def oxidation():  
    """
        Description:
            Calculations for cellular death

        Used:
            main()
    """  
    global broken_biofilm_agents
    global dead_bacterium_agents
    global cellular_vitality
    global biofilm_growth
   
    oxidation_vitality = {}
    total_oxidized_quantity = 0
    broken_biofilm = []
    dead_bacterium = []
    biofilm_growth = 'yes'
    csv_out_string = ""
    for molecule in membrane_chemical_concentrations:
        oxidation_vitality[molecule] = oxidized_chemical_quantities['oxidized %s' %(molecule)] * 100 / (oxidized_chemical_quantities['oxidized %s' %(molecule)] + membrane_chemical_concentrations[molecule])
        append_to_console_box('Oxidation proportion for %s:\t %s %%' %(molecule, oxidation_vitality[molecule]))
        csv_out_string = csv_out_string + "," + str(oxidation_vitality[molecule])

        
        if molecule not in graphed_oxidation_molecules:
            oxidation_graph_data.append([molecule, [oxidation_vitality[molecule]]])
            graphed_oxidation_molecules.append(molecule)
        else:
            index = 0
            for entry in oxidation_graph_data:
                if entry[0] == molecule:         
                    oxidation_graph_data[index][1].append(oxidation_vitality[molecule])
                index += 1

        if oxidation_vitality[molecule] >= stop_biofilm_threshold:
            biofilm_growth = 'no'
            broken_biofilm.append(molecule)
           
        if oxidation_vitality[molecule] >= oxidized_death_threshold:
            cellular_vitality = "dead"
            dead_bacterium.append(molecule)    
            
    csv_file.write("\n" + str(time) + csv_out_string)

    broken_biofilm_agents = ', '.join(biofilm_growth)
    dead_bacterium_agents = ', '.join(dead_bacterium)
    
    oxidation_graph_timesteps.append(time)

    
    
def final():    

    if cellular_vitality == "dead":  
        append_to_console_box('The < ' + str(dead_bacterium_agents) + ' > membrane component(s) caused bacterial death.')
        append_to_console_box("The bacterium has died at " + str(time) + " minutes.\nThe simulation is over.\n\n\n")
       
    else:                                                                                        
        if biofilm_growth == 'yes':
            append_to_console_box('The biofilm is still growing at ' + str(time) + ' minutes.\n')

        elif biofilm_growth == 'no':
            append_to_console_box('The < ' + broken_biofilm_agents + ' > membrane component(s) compromised biofilm maintenance.')
            append_to_console_box('The biofilm is has stopped growing at ' + str(time) + ' minutes.\n')
    
    plot_oxidation(oxidation_graph_timesteps, oxidation_graph_data)
    

       
def main():
    """
        Description:
            Execute the simulation by sequentially calling the aforementioned functions

        Used:
            base code
    """
    global time
    global log_file
    global csv_file
    global plot_out_path
   
    time = 0
    cell_geometry()
   
    while cellular_vitality == 'alive':
        root.update()
        reactions()
        oxidation()
        final()
        time += timestep
   
    log_file.close()
    csv_file.close()
    plt.savefig(plot_out_path, dpi=500) 
   
main()





#UI

window_width = 1280
window_height = 720
window_title = "Membrane Model"
background = "#FFFFFF"
primary = "#1c1c1c"
primary_variant = "#3700B3"
secondary = "#8a8a8a"
secondary_variant = "#1c1c1c"
header_frame_height = 50
body_frame_height = window_height - header_frame_height

root = Tk()
root.title(window_title)
#root.resizable(False, False)
root.geometry(str(window_width) + "x" + str(window_height))
root.configure(background=background)

root_frame = Frame(root, bg=background)
root_frame.pack()


#LEFT FRAME START -------------------------

left_frame = Frame(root_frame, width = window_width/5, height = window_height, bg= secondary)
left_frame.pack_propagate(0)
left_frame.grid(row = 0, column = 0)

left_frame_header = Label(left_frame, text="Control Panel", font="fixedsys 22", fg=primary, bg=secondary)
left_frame_header.pack(pady=8)





#LEFT FRAME END -------------------------


#RIGHT FRAME START ------------------------


right_frame = Frame(root_frame, width = window_width/5, height = window_height, bg= secondary)
right_frame.pack_propagate(0)
right_frame.grid(row = 0, column = 2)

right_frame_header = Label(right_frame, text="Environment", font="fixedsys 22", fg=primary, bg=secondary)
right_frame_header.pack(pady=8)


time_scale_label = Label(right_frame, text="Time step:", font="fixedsys 16", fg=primary, bg=secondary)
time_scale_label.place(x=10, y=50)

time_scale_entry = Entry(right_frame)
time_scale_entry.place(x=10, y=72)
time_scale_unit_label = Label(right_frame, text="minutes", font="fixedsys 16", fg=primary, bg=secondary)
time_scale_unit_label.place(x=120, y=72)


species_label = Label(right_frame, text="Species:", font="fixedsys 16", fg=primary, bg=secondary)
species_label.place(x=10, y=90+10)
species_choices = ['S aureus']
species_variable = StringVar(right_frame)
species_variable.set('S aureus')
species_choices_box = OptionMenu(right_frame, species_variable, *species_choices)
species_choices_box.place(x=10, y=115+10)


photosensitizer_label = Label(right_frame, text="Photosensitizer:", font="fixedsys 16", fg=primary, bg=secondary)
photosensitizer_label.place(x=10, y=130+10+10+10)
photosensitizer_choices = ['A3B-tetracationic zinc porphyrin', 'A3B-tetracationic copper porphyrin', 'A3B-dicationic zinc porphyrin', 'A3-monocationic gallium corrole']
photosensitizer_variable = StringVar(right_frame)
photosensitizer_variable.set('A3B-tetracationic zinc porphyrin')
photosensitizer_choices_box = OptionMenu(right_frame, photosensitizer_variable, *photosensitizer_choices)
photosensitizer_choices_box.place(x=10, y=155+10+10+10)


porphyrin_concentration_label = Label(right_frame, text="Porphyrin concentration:", font="fixedsys 16", fg=primary, bg=secondary)
porphyrin_concentration_label.place(x=10, y=190+30)
porphyrin_concentration_entry = Entry(right_frame)
porphyrin_concentration_entry.place(x=10, y=212+30)
porphyrin_concentration_unit_label = Label(right_frame, text="mg/L", font="fixedsys 16", fg=primary, bg=secondary)
porphyrin_concentration_unit_label.place(x=120, y=212+30)


light_excite_label = Label(right_frame, text="Excites Q/Soret bands?", font="fixedsys 16", fg=primary, bg=secondary)
light_excite_label.place(x=10, y=240+30-5)
excite_choices = ['Yes', 'No']
light_excite_variable = StringVar(right_frame)
light_excite_variable.set('Yes')
light_excite_box = OptionMenu(right_frame, light_excite_variable, *excite_choices)
light_excite_box.place(x=10, y=265+30-5)



luminous_intensity_label = Label(right_frame, text="Luminous Intensity:", font="fixedsys 16", fg=primary, bg=secondary)
luminous_intensity_label.place(x=10, y=290+30+5)
luminous_intensity_entry = Entry(right_frame)
luminous_intensity_entry.place(x=10, y=312+30+5)
luminous_intensity_unit_label = Label(right_frame, text="candela(s)", font="fixedsys 16", fg=primary, bg=secondary)
luminous_intensity_unit_label.place(x=120, y=312+30+5)



source_type_label = Label(right_frame, text="Light source?", font="fixedsys 16", fg=primary, bg=secondary)
source_type_label.place(x=10, y=370)
source_choices = ['LED', 'Incandescent']
source_variable = StringVar(right_frame)
source_variable.set('LED')
source_type_box = OptionMenu(right_frame, source_variable, *source_choices)
source_type_box.place(x=10, y=395)



source_wattage_label = Label(right_frame, text="Light source wattage:", font="fixedsys 16", fg=primary, bg=secondary)
source_wattage_label.place(x=10, y=430)
source_wattage_entry = Entry(right_frame)
source_wattage_entry.place(x=10, y=452)
source_wattage_unit_label = Label(right_frame, text="watts", font="fixedsys 16", fg=primary, bg=secondary)
source_wattage_unit_label.place(x=120, y=452)


#IMPORT BUTTON
start_button = Button(right_frame, text="Import Setup", command = lambda: import_environment (time_scale_entry, species_variable, photosensitizer_variable, porphyrin_concentration_entry, light_excite_variable, luminous_intensity_entry, source_variable, source_wattage_entry))
start_button.place(x = 30, y = 670)

#EXPORT BUTTON
start_button = Button(right_frame, text="Export Setup", command = lambda: export_environment (time_scale_entry, species_variable, photosensitizer_variable, porphyrin_concentration_entry, light_excite_variable, luminous_intensity_entry, source_variable, source_wattage_entry))
start_button.place(x = 140, y = 670)


#RIGHT FRAME END ------------------------

#START BUTTON
start_button = Button(left_frame, text="Start Simulation", command = lambda: initial_conditions (time_scale_entry, species_variable, photosensitizer_variable, porphyrin_concentration_entry, light_excite_variable, luminous_intensity_entry, source_variable, source_wattage_entry))
start_button.place(x = 20, y = 50)


#MIDDLE FRAME START --------
middle_frame = Frame(root_frame, width = window_width*3/5, height = window_height, bg = primary_variant)
middle_frame.grid(row = 0, column = 1)
#MIDDLE FRAME END ---------

header_frame = Frame(middle_frame, bg = background, width = window_width*3/5, height = 50, highlightthickness=3, highlightbackground=secondary_variant)
header_frame.pack_propagate(0)
header_frame.pack(side = TOP)



header = Label(header_frame, text="Membrane Model", font="fixedsys 33", fg=primary, bg=background)
header.pack()


body_frame = Frame(middle_frame, width = window_width*3/5, bg=background, height = body_frame_height, highlightthickness=3, highlightbackground=secondary_variant)
body_frame.pack_propagate(0)
body_frame.pack(side = BOTTOM)







#middle top frame for console outputs
middle_top_frame = Frame(body_frame, width = window_width*3/5, height = body_frame_height/2, bg = primary_variant)
middle_top_frame.pack_propagate(0)
middle_top_frame.pack()


scrollbar = Scrollbar(middle_top_frame)
console_box = Text(middle_top_frame, height=100, width= int(window_width*3/5), yscrollcommand=scrollbar.set, font=("Helvetica", 10), borderwidth=1, relief="solid")
scrollbar.config(command=console_box.yview)
scrollbar.pack(side=RIGHT, fill=Y)
console_box.pack(side="left")
console_box.configure(state='disabled')


#middle bottom frame for matplot lib graphs
middle_bottom_frame = Frame(body_frame, width = window_width*3/5, height = body_frame_height/3, bg = background)
middle_top_frame.pack_propagate(0)
middle_bottom_frame.pack()


fig, ax = plt.subplots()
ax.set_title('Fatty acid oxidation porportion')
ax.set_xlabel("Time (minutes)")
ax.set_ylabel("% Oxidiation (oxidized / total)")

# creating the Tkinter canvas
# containing the Matplotlib figure
canvas = FigureCanvasTkAgg(fig, master = middle_bottom_frame)  
canvas.draw()
 
# placing the canvas on the Tkinter window
canvas.get_tk_widget().pack()
 
# creating the Matplotlib toolbar
toolbar = NavigationToolbar2Tk(canvas, middle_bottom_frame)
toolbar.update()
 
# placing the toolbar on the Tkinter window
canvas.get_tk_widget().pack()





def plot_oxidation(oxidation_graph_timesteps_in, oxidation_graph_data_in):

    oxidized_death_threshold = 10
    ax.clear()
  
    for entry in oxidation_graph_data_in:
        ax.plot(oxidation_graph_timesteps_in, entry[1], label = entry[0])
        
    ax.plot(oxidation_graph_timesteps_in, [oxidized_death_threshold for i in range(len(oxidation_graph_timesteps_in))], label = "Oxidized death threshold")
    ax.set_title('Fatty acid oxidation porportion')
    ax.set_xlabel("Time (minutes)")
    ax.set_ylabel("% Oxidiation (oxidized / total)")
    ax.legend(loc='best', fontsize = 'small')

    canvas.draw()
    


    
    
    #.savefig('%s.%s' %(export_file_name, export_format))













#MIDDLE FRAME END ---------

testlabel = Label(body_frame,  text="Awaiting Start", font="fixedsys 40", fg=primary, bg=background)

testlabel.pack(pady=280)







root.mainloop()