# libraries for ui
from tkinter import *
from tkinter import ttk
import os
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from tkinter.filedialog import asksaveasfile, askopenfile
import threading

# import the software
import pdipy
pdi = pdipy.PDIBacterialPkg()
          
# interface options
header_string = ""
define_system_choices = [ "Solution depth + Solution volume", "Solution depth + Surface area", "Solution volume + Surface area"]
light_choices = ["irradiance (mW / cm²)", "exposure (J / cm²)", "lux (lumen / cm²)", "lumens (lumen)"]
data_inquiry_choices = ["% reduction", "time elapsed"]
hide_labels = []
hide_entries = []
hide_choices = []

# clear the console text
def clear_console():
    console_box.configure(state='normal')
    console_box.delete(1.0, "end")
    console_box.configure(state='disabled')
        
#UI functions
def console_message(message):
    console_box.configure(state='normal')
    console_box.insert("end", message + "\n")
    console_box.configure(state='disabled')    

# develop the UI
window_title = "iPDIpy"
background = "#FFFFFF"
primary = "#1c1c1c"
primary_variant = "#3700B3"
secondary = "#8a8a8a"
secondary_variant = "#1c1c1c"
header_frame_height = 50

root = Tk()
root.title(window_title)

#Setup resolution
scale = 1
window_width = int(root.winfo_screenwidth())*scale
window_height = int(root.winfo_screenheight())*scale
body_frame_height = window_height - header_frame_height
    
root.geometry(str(window_width) + "x" + str(window_height))
root.configure(background=background)
root.resizable(True, True)

for r in range(1000):
    root.rowconfigure(r, weight=1)    
for c in range(1000):
    root.columnconfigure(c, weight=1)

#LEFT FRAME START -------------------------
# parameter column content
left_frame = Frame(root, bg=secondary)
left_frame.grid(row = 1, column = 0, rowspan=999, columnspan = 400, stick = W+E+N+S)
#left_frame_header = Label(left_frame, text="Parameters", font="fixedsys 22", fg=primary, bg=secondary, borderwidth = 2, relief="solid")
left_frame_header_frame = Frame(root, bg=secondary, highlightthickness=3, highlightbackground=secondary_variant)
left_frame_header_frame.grid(row = 0, column = 0, rowspan=1, columnspan = 400, stick = W+E+N+S)
left_frame_header = Label(left_frame_header_frame, text="Parameters", font="fixedsys 22", fg=primary, bg=background)
left_frame_header.pack(fill=BOTH)

# header row content
left_top_header_frame = Frame(root, bg=secondary, highlightthickness=3, highlightbackground=secondary_variant)
left_top_header_frame.grid(row = 0, column = 400, rowspan = 1, columnspan = 250, stick = W+E+N+S)
header = Label(left_top_header_frame, text="Console", font="fixedsys 22", fg=primary, bg=background)
header.pack(fill=BOTH)

right_top_header_frame = Frame(root, bg=secondary, highlightthickness=3, highlightbackground=secondary_variant)
right_top_header_frame.grid(row = 0, column = 650, rowspan = 1, columnspan = 350, stick = W+E+N+S)
header = Label(right_top_header_frame, text="Graph", font="fixedsys 22", fg=primary, bg=background)
header.pack(fill=BOTH)

# console column content
console_frame = Frame(root, bg=primary_variant)
console_frame.grid(row = 1, column = 400, rowspan = 999, columnspan = 250, stick = W+E+N+S)
scrollbar = Scrollbar(console_frame)
console_box = Text(console_frame, width=1, height=1, yscrollcommand=scrollbar.set, font=("Helvetica", 10), borderwidth=1, relief="solid")
scrollbar.config(command=console_box.yview)
scrollbar.pack(side=RIGHT, fill=BOTH)
console_box.pack(fill=BOTH, expand=True)
console_box.configure(state='disabled')

# graph column content
graph_frame = Frame(root, bg=secondary)
graph_frame.grid(row = 1, column = 650, rowspan=999, columnspan = 350, stick = W+E+N+S)
graph_sub_frame = Frame(graph_frame, bg=secondary)
graph_sub_frame.pack_propagate(0)
graph_sub_frame.pack(expand=True, fill=BOTH)

fig, ax = plt.subplots()
ax.set_title('Oxidized proportion of cytoplasmic fatty acids')
ax.set_xlabel("Time (min)")
ax.set_ylabel("Oxidized proportion")

# creating the Tkinter canvas of the Matplotlib figure
canvas = FigureCanvasTkAgg(fig, master = graph_sub_frame)  
canvas.get_tk_widget().pack()
canvas.draw()

# creating the Matplotlib toolbar and placing it in the Tkinter window
toolbar = NavigationToolbar2Tk(canvas, graph_sub_frame)
toolbar.update()
canvas.get_tk_widget().pack(expand=True, fill=BOTH)


define_system_labels = []
define_system_entries = []
define_system_unit_labels = []

y_increment_before_define_system = 0
x_increment_hide = 10            
    
# define the data elements
def hide(choice):
    '''
    The options for either the surface system or the solution system are defined.
    '''
    global y_increment_before_define_system
    global data_elements
    global system_input_choice_label
    global system_input_option_menu
    global photosensitizer_concentration_moles_per_area_label
    global photosensitizer_concentration_moles_per_area_entry
    global photosensitizer_concentration_moles_per_area_unit_label
    global photosensitizer_concentration_moles_per_l_label
    global photosensitizer_concentration_moles_per_l_entry
    global photosensitizer_concentration_moles_per_l_unit_label
    global hide_labels
    global hide_entries
    global hide_choices
    global variate_layout_x
    global variate_layout_y
    
    def forget_place():
        global system_input_choice_label
        global system_input_option_menu

        for i in range(len(define_system_labels)):
            define_system_labels[i].place_forget()
            define_system_entries[i].place_forget()
            define_system_unit_labels[i].place_forget()
            
        for j in [hide_labels, hide_entries, hide_choices]:
            for i in range(len(j)):
                j[i].place_forget()
        
        system_input_option_menu.place_forget()
        system_input_choice_label.place_forget()
        photosensitizer_concentration_moles_per_area_label.place_forget()
        photosensitizer_concentration_moles_per_area_entry.place_forget()
        photosensitizer_concentration_moles_per_area_unit_label.place_forget()
        photosensitizer_concentration_moles_per_l_label.place_forget()
        photosensitizer_concentration_moles_per_l_entry.place_forget()
        photosensitizer_concentration_moles_per_l_unit_label.place_forget()
        
    
    x = x_increment_hide
    y_increment_entry = 45
    y_increment_choice = 60
    y = y_increment_before_define_system
    variate_layout_flag = False
    
    if x == variate_layout_x and y == variate_layout_y:
        variate_layout_flag = True
   
    forget_place()
    
    first_entry = 0
    second_entry = 0
    
    system_choice = ""
    system_input_choice = ""
    
    for element in data_elements:
        if element["variable_name"] == "System input choice":
            system_input_choice = str(element["get_element"].get())
        elif element["variable_name"] == "System choice":
            system_choice = str(element["get_element"].get())

    
    if system_choice == "Solution":
        
        system_input_choice_label.place(x=x, y=y)
        system_input_option_menu.place(x=x, y=y+25)
        
        y += y_increment_choice
        
        if system_input_choice == define_system_choices[0]:
            first_entry = 1
            second_entry = 2
        elif system_input_choice == define_system_choices[1]:
            first_entry = 1
            second_entry = 0
        elif system_input_choice == define_system_choices[2]:
            first_entry = 2
            second_entry = 0
    
        define_system_labels[first_entry].place(x=x, y=y)
        define_system_entries[first_entry].place(x=x, y=y+22)
        define_system_unit_labels[first_entry].place(x=x+110, y=y+22)
    
        y += y_increment_entry
        
        define_system_labels[second_entry].place(x=x, y=y)
        define_system_entries[second_entry].place(x=x, y=y+22)
        define_system_unit_labels[second_entry].place(x=x+110, y=y+22)
        
        y += y_increment_entry
        
        photosensitizer_concentration_moles_per_l_label.place(x=x, y=y)
        photosensitizer_concentration_moles_per_l_entry.place(x=x, y=y+22)
        photosensitizer_concentration_moles_per_l_unit_label.place(x=x+125, y=y+22)
        
        y += y_increment_entry
        
        
        
    elif system_choice == "Surface":
        
        define_system_labels[0].place(x=x, y=y)
        define_system_entries[0].place(x=x, y=y+22)
        define_system_unit_labels[0].place(x=x+110, y=y+22)
        
        y += y_increment_entry
        
        photosensitizer_concentration_moles_per_area_label.place(x=x, y=y)
        photosensitizer_concentration_moles_per_area_entry.place(x=x, y=y+22)
        photosensitizer_concentration_moles_per_area_unit_label.place(x=x+125, y=y+22)
        
        y += y_increment_entry
    
    if not variate_layout_flag:
        hide_labels[0].place(x=x, y=y)
        hide_entries[0].place(x=x, y=y+22)
        hide_choices[0].place(x=x+130, y=y+17)
    else:
        hide_labels[0].place(anchor=CENTER, relx=0.5, y=y+40-10, x=-30-20)
        hide_entries[0].place(anchor=CENTER, relx=0.5, y=y+40+22-10, x=-30-10)
        hide_choices[0].place(anchor=CENTER, relx=0.5, y=y+40+20-10, x=90-10)
        
def pack_inputs():
    global y_increment_before_define_system
    global system_input_choice_label
    global system_input_option_menu
    global photosensitizer_concentration_moles_per_area_label
    global photosensitizer_concentration_moles_per_area_entry
    global photosensitizer_concentration_moles_per_area_unit_label
    global photosensitizer_concentration_moles_per_l_label
    global photosensitizer_concentration_moles_per_l_entry
    global photosensitizer_concentration_moles_per_l_unit_label
    global hide_labels
    global hide_entries
    global hide_choices
    
    
    def ui_choice(frame, choices, default, x, y, command=None):
        choice_variable = StringVar(frame)
        choice_variable.set(default)
        if command:
            OptionMenu(frame, choice_variable, *choices, command=command).place(x=x, y=y)
        else:
            OptionMenu(frame, choice_variable, *choices).place(x=x, y=y)
        return choice_variable
    
    # coordinate values
    x = 10
    y = 10
    y_increment_choice = 60
    y_increment_entry = 45
    
    # parameter types
    data_type = ["num", "string", "string", "num", "string", "num", "num", "string", 'string', "num", "num", "num"]
    
    # ========== the first set of parameters ========== 
    label_text = ["Simulation time:", "Bacterial species:", "Photosensitizer:", "Photosensitizer concentration:", "Light source:", "Molecular proportion:"]
    variable_name = ["Simulation time", "Bacterial species", "Photosensitizer", "Photosensitizer concentration", "Light source", "Molecular proportion"]
    input_type = ["entry", "choice", "choice", "entry", "choice", "entry"]
    choice_increment = 0
    choices = [['Saureus'],
               ['A3B_4Zn', 'A3B-tetracationic zinc porphyrin', 'A3B-tetracationic copper porphyrin', 'A3B-dicationic zinc porphyrin', 'A3-monocationic gallium corrole'],
               ['LED', 'incandescent', 'fluorescent']]
    unit_label = ["minutes", "n/a", "n/a", "molar", "n/a", "n/a"]
    default_choices = ['Saureus', 'A3B_4Zn', 'LED']
    default_values = ["", 'Saureus', 'A3B_4Zn', "", 'LED', ""]
    
    data_elements = []
    for i in range(len(label_text)):
        data_element = {"input_type":input_type[i], "variable_name":variable_name[i], "data_type":data_type[i], "get_element":"", "default_value":default_values[i]}
        temp_label = Label(left_frame, text=label_text[i], font="fixedsys 16", fg=primary, bg=secondary)
        temp_label.place(x=x, y=y)
        if input_type[i] == "entry":
            temp_entry = Entry(left_frame)
            temp_entry.place(x=x, y=y+22)
            if unit_label[i] != "n/a":
                temp_unit_label = Label(left_frame, text=unit_label[i], font="fixedsys 16", fg=primary, bg=secondary)
                temp_unit_label.place(x=x+110, y=y+22)
            y += y_increment_entry
            data_element["get_element"] = temp_entry
        elif input_type[i] == "choice":
            temp_variable = ui_choice(left_frame, choices[choice_increment], default_choices[choice_increment], x, y+25)
            data_element["get_element"] = temp_variable
            y += y_increment_choice
            choice_increment += 1
        data_elements.append(data_element)
        
    # ========== Light input ========== 
    light_input_label = Label(left_frame, text="Light magnitude:", font="fixedsys 16", fg=primary, bg=secondary)
    light_input_label.place(x=x, y=y)
    light_input_entry = Entry(left_frame)
    light_input_entry.place(x=x, y=y+22)
    temp_variable = ui_choice(left_frame, light_choices, light_choices[0], x+130, y+17)
    data_element = {"input_type":"choice", "variable_name":"Light choice", "data_type":"string", "get_element":temp_variable, "default_value":light_choices[0]}
    data_elements.append(data_element)
    data_element = {"input_type":"entry", "variable_name":"Light value", "data_type":"num", "get_element":light_input_entry, "default_value":""}
    data_elements.append(data_element)
    
    
    

    # ========== System input choice ========== 
    system_choices = ["Solution", "Surface"]
    y += y_increment_entry
    
    system_choice_label = Label(left_frame, text="System type:", font="fixedsys 16", fg=primary, bg=secondary)
    system_choice_label.place(x=x, y=y)
    temp_variable = ui_choice(left_frame, system_choices, system_choices[1], x, y+25, command=hide)
    data_element = {"input_type":"choice", "variable_name":"System choice", "data_type":"string", "get_element":temp_variable, "default_value":system_choices[1]}
    data_elements.append(data_element)
        
    # System inputs per the system choice
    y += y_increment_choice
    y_increment_before_define_system = y
    
    system_input_choice_label = Label(left_frame, text="System inputs:", font="fixedsys 16", fg=primary, bg=secondary)
    #system_input_choice_label.place(x=x, y=y)
    temp_variable = StringVar(left_frame)
    temp_variable.set(define_system_choices[0])
    system_input_option_menu = OptionMenu(left_frame, temp_variable, *define_system_choices, command=hide)
    #system_input_option_menu.place(x=x, y=y+25)
    data_element = {"input_type":"choice", "variable_name":"System input choice", "data_type":"string", "get_element":temp_variable, "default_value":define_system_choices[0]}
    data_elements.append(data_element)
        
    #y += y_increment_choice
    
    temp_labels = ["Surface area:", "Solution depth:", "Solution volume:"]
    temp_variable_names = ["Surface area", "Solution depth", "Solution volume"]
    for label in temp_labels:
        define_system_labels.append(Label(left_frame, text=label, font="fixedsys 16", fg=primary, bg=secondary))
    
    temp_entry_count = 3
    for entry in range(temp_entry_count):
        temp_entry = Entry(left_frame)
        define_system_entries.append(temp_entry)
        data_element = {"input_type":"entry", "variable_name":temp_variable_names[entry], "data_type":"num", "get_element":temp_entry , "default_value":""}
        data_elements.append(data_element)
    
    temp_unit_labels = ["m²", "m", "m³"]
    for label in temp_unit_labels:
        define_system_unit_labels.append(Label(left_frame, text=label, font="fixedsys 16", fg=primary, bg=secondary))
    
    """
    define_system_labels[1].place(x=x, y=y)
    define_system_entries[1].place(x=x, y=y+22)
    define_system_unit_labels[1].place(x=x+110, y=y+22)

    y += y_increment_entry
    
    define_system_labels[2].place(x=x, y=y)
    define_system_entries[2].place(x=x, y=y+22)
    define_system_unit_labels[2].place(x=x+110, y=y+22)
    """
    
    photosensitizer_concentration_moles_per_area_label = Label(left_frame, text="Photosensitizer concentration:", font="fixedsys 16", fg=primary, bg=secondary)
    photosensitizer_concentration_moles_per_area_unit_label = Label(left_frame, text="moles per m²", font="fixedsys 16", fg=primary, bg=secondary)
    photosensitizer_concentration_moles_per_area_entry = Entry(left_frame)
    
    data_element = {"input_type":"entry", "variable_name":"Photosensitizer concentration moles", "data_type":"num", "get_element":photosensitizer_concentration_moles_per_area_entry, "default_value":""}
    data_elements.append(data_element)
    
    photosensitizer_concentration_moles_per_l_label = Label(left_frame, text="Photosensitizer concentration:", font="fixedsys 16", fg=primary, bg=secondary)
    photosensitizer_concentration_moles_per_l_unit_label = Label(left_frame, text="moles per liter", font="fixedsys 16", fg=primary, bg=secondary)
    photosensitizer_concentration_moles_per_l_entry = Entry(left_frame)
    
    """
    y += y_increment_entry
    
    photosensitizer_concentration_moles_per_l_label.place(x=x, y=y)
    photosensitizer_concentration_moles_per_l_entry.place(x=x, y=y+22)
    photosensitizer_concentration_moles_per_l_unit_label.place(x=x+125, y=y+22)
    """
    data_element = {"input_type":"entry", "variable_name":"Photosensitizer concentration moles per liter", "data_type":"num", "get_element":photosensitizer_concentration_moles_per_l_entry, "default_value":""}
    data_elements.append(data_element)
    
    define_system_labels[0].place(x=x, y=y)
    define_system_entries[0].place(x=x, y=y+22)
    define_system_unit_labels[0].place(x=x+110, y=y+22)
    
    y += y_increment_entry
    
    photosensitizer_concentration_moles_per_area_label.place(x=x, y=y)
    photosensitizer_concentration_moles_per_area_entry.place(x=x, y=y+22)
    photosensitizer_concentration_moles_per_area_unit_label.place(x=x+125, y=y+22)
    
    # ========== Data inquiry input ========== 
    y += y_increment_entry
    data_type_label = Label(left_frame, text="Data inquiry:", font="fixedsys 16", fg=primary, bg=secondary)
    data_type_label.place(x=x, y=y)
    hide_labels.append(data_type_label)
    data_type_entry = Entry(left_frame)
    data_type_entry.place(x=x, y=y+22)
    hide_entries.append(data_type_entry)
    
    #temp_variable = ui_choice(left_frame, data_inquiry_choices, data_inquiry_choices[0], x+130, y+17)
    
    
    temp_variable = StringVar(left_frame)
    temp_variable.set(data_inquiry_choices[0])
    data_inquiry_option_menu = OptionMenu(left_frame, temp_variable, *data_inquiry_choices, command=hide)
    data_inquiry_option_menu.place(x=x+130, y=y+17)
    hide_choices.append(data_inquiry_option_menu)
    
    data_element = {"input_type":"choice", "variable_name":"Data inquiry choice", "data_type":"string", "get_element":temp_variable, "default_value":data_inquiry_choices[0]}
    data_elements.append(data_element)
    data_element = {"input_type":"entry", "variable_name":"Data inquiry value", "data_type":"num", "get_element":data_type_entry, "default_value":""}
    data_elements.append(data_element)
 
    
    
    return data_elements 

def create_header_string(data_elements):
    global header_string
    
    first_element = True
    for element in data_elements:
        if first_element:
            header_string = element["variable_name"]
            first_element = False
        else:
            header_string = header_string + "," + element["variable_name"]

    header_string = header_string + "\n"
    
    
data_elements = pack_inputs()    
create_header_string(data_elements)


# IMPORT BUTTON
def import_environment(data_elements):
    f = askopenfile(mode = 'r', filetypes =[("Text Files", "*.txt")])
    if f is not None:
        content = f.read()
        try:
            if header_string in content:
                content = content.replace(header_string, "").replace("\n", "").replace("\r", "")
                fields = content.split(",")
                print(fields)
                counter = 0
                
                for element in data_elements:
                    print(element["variable_name"])
                    print(str(fields[counter]))
                    if element["input_type"] == "entry":
                        element["get_element"].delete(0, END)
                        element["get_element"].insert(0, str(fields[counter]))
                               
                    elif element["input_type"] == "choice":
                        element["get_element"].set(str(fields[counter]))
                    counter += 1    
                console_message("Parameters successfully imported.")
        except:
            #console_message("ERROR: Failed to import: " + str(os.path.basename(f.name)) + ". File malformed.")
            console_message("ERROR: The parameter file is unable to be imported.")   


import_button = Button(left_frame, text="Import parameters", command = lambda: import_environment(data_elements))
import_button.place(relx = 1, x = -140, y = 105+65-40-65+25-40, anchor = "nw")

            
# EXPORT BUTTON
def export_environment(data_elements):
    global header_string
    f = asksaveasfile(mode='w', defaultextension=".txt")
    if f is None:
        return
 
    f.write(header_string)
    data_string = ""
    first_element = True
    for element in data_elements:
        if first_element:
            data_string = str(element["get_element"].get())
            first_element = False
        else:
            data_string = data_string + ","  + str(element["get_element"].get())

    f.write(data_string)
    f.close()
    console_message("Parameters successfully exported.")

    
export_button = Button(left_frame, text="Export parameters", command = lambda: export_environment(data_elements))
export_button.place(relx = 1, x = -140, y = 140+65-40-65+25-40, anchor = "nw")


# START BUTTON
def start_simulation(data_elements):
    #verify the ui inputs
    def verify_ui_inputs(values, expected_types, input_names):
        def verify_ui_input(value, expected_type, input_name):
            if expected_type == "string":
                try:
                    temp = str(value)
                    return True
                except:
                    console_message("INPUT ERROR: " + input_name + " should be a number.")
                    return False
            elif expected_type == "num":
                try:
                    temp = float(value)
                    return True
                except:
                    console_message("INPUT ERROR: " + input_name + " should be a number.")
                    return False

        inputs_good = True

        for i in range(len(values)):    
            input_is_good = verify_ui_input(values[i], expected_types[i], input_names[i])
            if not input_is_good:
                inputs_good = False

        return inputs_good


    clear_console()    
    system_choice = ""
    system_choice_input_flag = -1
    
    for element in data_elements:
        if element["variable_name"] == "System input choice":
            for j in range(len(define_system_choices)):
                if element["get_element"].get() == define_system_choices[j]:
                    system_choice_input_flag = j
        elif element["variable_name"] == "System choice":
            system_choice = str(element["get_element"].get())
        
    get_list = []
    get_executed_list = []
    data_type_list = []
    variable_name_list = []
    
    for element in data_elements:
        #if element["variable_name"] == define_system_choices[choice_flag]:
        
        if system_choice == "Solution":
            if element["variable_name"] == "Surface area" and (system_choice_input_flag == 0):
                continue
            elif element["variable_name"] == "Solution depth" and (system_choice_input_flag == 2):
                continue
            elif element["variable_name"] == "Solution volume" and (system_choice_input_flag == 1):
                continue
            elif element["variable_name"] == "Photosensitizer concentration moles":
                continue
        elif system_choice == "Surface":
            if element["variable_name"] == "Solution depth":
                continue
            elif element["variable_name"] == "Solution volume":
                continue
            elif element["variable_name"] == "Photosensitizer concentration moles per liter":
                continue
            
        get_list.append(element["get_element"])
        get_executed_list.append(element["get_element"].get())
        data_type_list.append(element["data_type"])
        variable_name_list.append(element["variable_name"])
    
    inputs_good = verify_ui_inputs(get_executed_list, data_type_list, variable_name_list)
    
    if inputs_good:
        print("inputs were good")
        proccessed_data_list = {}
        for i in range(len(get_list)):
            if data_type_list[i] == "num":
                proccessed_data_list[variable_name_list[i]] = float(get_list[i].get())
            elif data_type_list[i] == "string":
                proccessed_data_list[variable_name_list[i]] = get_list[i].get()
            if proccessed_data_list[variable_name_list[i]] == "":
                proccessed_data_list[variable_name_list[i]] = None
        
        #t = threading.Thread(target=run_simulation, args=(proccessed_data_list["Simulation time"], proccessed_data_list["Bacterial species"], proccessed_data_list["Photosensitizer"], proccessed_data_list["Photosensitizer concentration"], proccessed_data_list["Light source"], proccessed_data_list["Surface area"], proccessed_data_list["Solution volume"], proccessed_data_list["Solution depth"], proccessed_data_list["Photosensitizer surface area"], proccessed_data_list["Molecular proportion"], proccessed_data_list["Photosensitizer moles per square cm"]))
        #t.start()
        t = threading.Thread(target=run_simulation, args=(proccessed_data_list,))
        t.start()
        
def run_simulation(proccessed_data_list):    
    system_choice_flag = -1
    for j in range(len(define_system_choices)):
        if proccessed_data_list["System input choice"] == define_system_choices[j]:
            system_choice_flag = j    
    
    # define the system conditions depending upon the system perspective
    if proccessed_data_list["System choice"] == "Solution":
        if system_choice_flag == 0:
            pdi.define_system(solution_depth = proccessed_data_list["Solution depth"], solution_volume=proccessed_data_list["Solution volume"])
        elif system_choice_flag == 1:
            pdi.define_system(solution_depth = proccessed_data_list["Solution depth"], surface_area = proccessed_data_list["Surface area"])
        elif system_choice_flag == 2:
            pdi.define_system(solution_volume = proccessed_data_list["Solution volume"], surface_area = proccessed_data_list["Surface area"])
    elif proccessed_data_list["System choice"] == "Surface":
        pdi.define_system(surface_area = proccessed_data_list["Surface area"], photosensitizer_surface_area=proccessed_data_list["Photosensitizer concentration moles"])
    
    # define the bacterium
    pdi.define_bacterium(bacterial_specie=proccessed_data_list["Bacterial species"])
    light_choice_flag = -1
    for j in range(len(light_choices)):
        if proccessed_data_list["Light choice"] == light_choices[j]:
            light_choice_flag = j
        
    # define the photosensitizer
    pdi.define_photosensitizer(photosensitizer_molar=proccessed_data_list["Photosensitizer concentration"], photosensitizer=proccessed_data_list["Photosensitizer"])
    
    # define the light
    if light_choice_flag == 0:
        pdi.define_light(light_source = proccessed_data_list["Light source"], simulation_time=proccessed_data_list["Simulation time"], irradiance=proccessed_data_list["Light entry"])
    elif light_choice_flag == 1:
        pdi.define_light(light_source = proccessed_data_list["Light source"], simulation_time=proccessed_data_list["Simulation time"], exposure=proccessed_data_list["Light entry"])
    elif light_choice_flag == 2:
        pdi.define_light(light_source = proccessed_data_list["Light source"], simulation_time=proccessed_data_list["Simulation time"], lumen=proccessed_data_list["Light entry"])
    elif light_choice_flag == 3:
        pdi.define_light(light_source = proccessed_data_list["Light source"], simulation_time=proccessed_data_list["Simulation time"], lux=proccessed_data_list["Light entry"])
    
    # conduct the simulation and export the results
    pdi.singlet_oxygen_calculations()
    pdi.kinetic_calculation()
    pdi.export()
            

button_place_y_pos = 70      
buton_y_pos_increment = 35  
        
start_button = Button(left_frame, text="Start Simulation", command = lambda: start_simulation(data_elements))
start_button.place(relx = 1, x = -140, y = button_place_y_pos-40+25-40, anchor = "nw")
button_place_y_pos += buton_y_pos_increment

# CLEAR BUTTON
def clear_environment(data_elements):
    clear_console()
    counter = 0
    for element in data_elements:
        if element["input_type"] == "entry":
            element["get_element"].delete(0, END)
            element["get_element"].insert(0, str(element["default_value"]))
            
        elif element["input_type"] == "choice":
            element["get_element"].set(str(element["default_value"]))
            counter += 1
    
    console_message("Parameters successfully cleared.")


clear_button = Button(left_frame, text="Clear parameters", command = lambda: clear_environment(data_elements))
clear_button.place(relx = 1, x = -140, y = 175+65-40-65+25-40, anchor = "nw")


# PARSE BUTTON
def parse_data(data_elements):
    target_reduction = 0.98
    if target_reduction > pdi.processed_data.iloc[-1]:
        print('never reached')

    for index, point in pdi.processed_data.iterrows():
        if index > 6:
            print(point['oxidation_proportion'])
        if point['oxidation_proportion'] >= target_reduction:
            print('time to target (hr): ', index)
            break
    

parse_button = Button(left_frame, text="Parse data", command = lambda: parse_data(data_elements))
parse_button.place(relx = 1, x = -140, y = 210+65-40-65+25-40, anchor = "nw")

def export_console():
    f = asksaveasfile(mode='w', defaultextension=".txt")
    if f is None:
        return
 
    out = console_box.get("1.0", "end")
    print(out)
    f.write(out)
    
    f.close()

export_console_button = Button(left_frame, text="Export console", command = lambda: export_console())
export_console_button.place(relx = 1, x = -140, y = 245+65-40-65+25-40, anchor = "nw")

saved_y_increment_before_define_system = y_increment_before_define_system

variate_layout_x = 310
variate_layout_y = 320

def variate_layout():
    global saved_y_increment_before_define_system
    global y_increment_before_define_system
    global x_increment_hide
    global variate_layout_x
    
    if x_increment_hide == 10:
        x_increment_hide = variate_layout_x
        saved_y_increment_before_define_system = y_increment_before_define_system
        y_increment_before_define_system = variate_layout_y
    else:
        y_increment_before_define_system = saved_y_increment_before_define_system
        x_increment_hide = 10
    hide(data_elements)

variation_button = Button(left_frame, text="Change Layout", command = lambda: variate_layout())
variation_button.place(relx = 1, x = -140, y = 280+65-40-65+25-40, anchor = "nw")

root.mainloop()