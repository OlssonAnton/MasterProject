import numpy as np

def number_of_points(mvd_list, mtad_list):
    nr_of_points = 0
    for mtad in mtad_list:
        for mvd in mvd_list:
            if (mvd < mtad):
                nr_of_points += 1
    return nr_of_points

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

def read_and_write(file, read, write):
    new_file_content = ""
    reading_file = open(file, "r")
    for line in reading_file:
        stripped_line = line.strip()
        new_line = stripped_line.replace(read, write)
        new_file_content += new_line + "\n"
    reading_file.close()
    writing_file = open(file, "w")
    writing_file.write(new_file_content)
    writing_file.close()

def inclusive_range(start, end, step):
    return np.arange(start, end+step, step).tolist()

def asimov_significance(luminosity, signal, background): #signal/background in either cross sections or number of events
    #return luminosity**(0.5) * signal / ((signal + background)**(0.5))
    return signal / (signal + background)**0.5

def generate_mvd_mtad_lists(mvd_min, mvd_step, mtad_list):
    mvd_list = []
    mtad_long_list = []
    mvd = mvd_min
    for mtad in mtad_list:
        while(mvd < mtad):
            mvd_list.append(mvd)
            mvd += mvd_step
            mtad_long_list.append(mtad)
        mvd = mvd_min
    return [mvd_list, mtad_long_list]

def z_list(L,sigma_s, sigma_b): #Significance (z = 2 corresponds to 95% C.L.)
    z_list = []
    for i in range(len(sigma_s)):
        z_list.append(L**(0.5) * sigma_s[i] / ((sigma_s[i] + sigma_b)**(0.5)))
    return z_list
