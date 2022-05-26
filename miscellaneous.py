import numpy as np
import math

#Help function that reads lists of data in string format and return lists of floats
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

#Help function that replaces the contents of a save file at file_path with new_line
def edit_line(file_path, line_index, new_line):
    opened_script_file = open(file_path)
    string_list = opened_script_file.readlines()
    opened_script_file.close()
    string_list[line_index] = new_line + '\n'
    opened_script_file = open(file_path, 'w')
    new_file_contents = ''.join(string_list)
    opened_script_file.write(new_file_contents)
    opened_script_file.close()

#Redefined range() function that includes the final element (end) in the output
def inclusive_range(start, end, step):
    return np.arange(start, end+step, step).tolist()

# Basic way of computing significance of signal compared to background
def asimov_significance(S, B, sigma_b): #(signal and background in nr. of events)
    return math.sqrt(2 * ((S + B) * math.log(((S+B)*(B + sigma_b**2))/(B**2 + (S + B)*sigma_b**2))
    - (B**2/sigma_b**2) * math.log(1 + sigma_b**2 * S / (B*(B + sigma_b**2)))))
    #return signal / (signal + background)**0.5
