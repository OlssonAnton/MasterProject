import numpy as np

def read_save_file(file, line_index):
    file = open(file, 'r')
    data = file.readlines()[line_index]
    fixed_data = ''.join(data.split('[', 1))
    fixed_data = ''.join(fixed_data.split(']', 1))
    fixed_data = fixed_data.replace(',', '')
    list_of_floats = list(fixed_data.split(' '))
    for i in range(len(list_of_floats)):
        list_of_floats[i] = float(list_of_floats[i])
    return list_of_floats

def edit_line(file_path, line_index, new_line):
    opened_script_file = open(file_path)
    string_list = opened_script_file.readlines()
    opened_script_file.close()
    string_list[line_index] = new_line + '\n'
    opened_script_file = open(file_path, 'w')
    new_file_contents = ''.join(string_list)
    opened_script_file.write(new_file_contents)
    opened_script_file.close()

def inclusive_range(start, end, step):
    return np.arange(start, end+step, step).tolist()

def asimov_significance(luminosity, signal, background): #signal/background in either cross sections or number of events
    #return luminosity**(0.5) * signal / ((signal + background)**(0.5))
    return signal / (signal + background)**0.5
