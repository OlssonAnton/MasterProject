import os
import subprocess
import numpy as np
import math
from general_stuff import read_save_file, edit_line

#Function that calculates the relic densities for all points of a simulation
#with output stored at scan_path
def relic_density(scan_path, dark_coupling):
    saved_scan_data_path = scan_path + '/saved_scan_data'
    saved_scan_data_file = open(saved_scan_data_path, 'r')
    saved_data = saved_scan_data_file.readlines()
    mvd_properties = read_save_file(saved_scan_data_path, 1)
    mtad_properties = read_save_file(saved_scan_data_path, 2)
    saved_scan_data_file.close()
    mvd_min = math.floor(mvd_properties[0]) #convert float to int
    mvd_max = math.floor(mvd_properties[1])
    mvd_step = math.floor(mvd_properties[2])
    mtad_min = math.floor(mtad_properties[0])
    mtad_max = math.floor(mtad_properties[1])
    mtad_step = math.floor(mtad_properties[2])

    mvd_length = 1 + math.floor((mvd_max - mvd_min) / mvd_step)
    mtad_length = 1 + math.floor((mtad_max - mtad_min) / mtad_step)
    relic_density_matrix = np.zeros((mvd_length, mtad_length), float) #Initialize relic_density matrix with zeros

    saved_relic_densities = open(saved_scan_data_path).readlines()[6]
    relic_saved = False
    if not saved_relic_densities == '[]\n': #Check if the relic densities were already retrieved and saved to file
        relic_saved = True
        print("Relic densities are saved, use old ones")
        #Modify the saved string to build the numpy array
        saved_relic_densities = saved_relic_densities.replace('] [', '][') #Remove spaces
        saved_relic_densities = saved_relic_densities.replace('[[', '[') #Remove double brackets at start
        saved_relic_densities = saved_relic_densities.replace(']]', ']') #Remove double brackets at end
        saved_relic_densities = saved_relic_densities.split(']') #Create columns by spliting the rows
        saved_relic_densities.pop() #Removes final linebreak (last element after split())
        for i in range(len(saved_relic_densities)):
            line = saved_relic_densities[i] #Retrieve a row with relic densities
            line = line.replace('[', '') #Remove the left bracket from the row
            line = line.replace('0.     ', '0. ') #Remove spaces from np.zeros() initializer
            while '  ' in line:
                line = line.replace('  ', ' ')
            line = line.split(' ') #Data in string separated by spaces
            while(line[-1] == ''):
                line.pop() #Remove all white spaces at the end of lines
            for j in range(len(line)):
                relic_density = line[j]
                relic_density = float(relic_density) #Convert from string-number to float-number
                relic_density_matrix[i][j] = relic_density #Add relic density to correct place in relic density 2d array
    if relic_saved == False:
        data_file_path = './micromegas_5.2.13/masterproject_FG/work/models/vars1.mdl' #From madgraph/various directory
        edit_line(data_file_path, 3, 'gD             |'+str(0.2)+'         |SU2D coupling constant      ')
        edit_line(data_file_path, 36, 'Mtap           |800          |Mass of tap.')
        edit_line(data_file_path, 38, 'Mh2            |300          |Mass of h2.')

        mtad = mtad_min
        mvd = mvd_min
        while(mtad <= mtad_max):
            while(mvd < mtad and mvd <= mvd_max):
                mvd_index = math.floor((mvd - mvd_min) / mvd_step)
                mtad_index = math.floor((mtad - mtad_min) / mtad_step)
                edit_line(data_file_path, 9, 'MVD            |'+str(mvd)+'          |Mass of GDP.')
                edit_line(data_file_path, 33, 'MtaD           |'+str(mtad)+'          |Mass of taD.')
                os.system('make -C ./micromegas_5.2.13/masterproject_FG main=main.cpp')
                output = subprocess.run(['micromegas_5.2.13/masterproject_FG/main', 'micromegas_5.2.13/masterproject_FG/data.par'], capture_output = True).stdout
                output = str(output)
                channel_index_first = output.index('Channels which contribute') - 4
                omega_index_last = output.index('Omega=') + 6
                output = output[omega_index_last:channel_index_first] #Slice everything but the float
                output = float(output) #Convert string with exponential representation to float number
                relic_density_matrix[mvd_index][mtad_index] = output
                mvd += mvd_step
            mvd = mvd_min
            mtad += mtad_step


        #Save relic densities to save file
        numpy_string = np.array2string(relic_density_matrix).replace(']\n', ']')
        numpy_string = numpy_string.replace('\n', '')
        numpy_string = numpy_string.replace('  ', ' ')
        edit_line(scan_path + '/saved_scan_data', 6, numpy_string)

    return relic_density_matrix
