import math
import numpy as np
import os
from os.path import exists
import lhe_parser
import gzip
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from general_stuff import read_save_file, inclusive_range, edit_line, asimov_significance
from significance_script import significance_calculator
from relic_density import relic_density

def retrieve_data(scan_path):
    saved_scan_data_path = scan_path + '/saved_scan_data' # Define path to save file for scan
    mvd_properties = read_save_file(saved_scan_data_path, 1) #Retrieve mass data from the save file
    mtad_properties = read_save_file(saved_scan_data_path, 2)
    #Masses saved as floats are converted to ints with floor()
    mvd_min = math.floor(mvd_properties[0])
    mvd_max = math.floor(mvd_properties[1])
    mvd_step = math.floor(mvd_properties[2])
    mtad_min = math.floor(mtad_properties[0])
    mtad_max = math.floor(mtad_properties[1])
    mtad_step = math.floor(mtad_properties[2])

    mvd_length = 1 + math.floor((mvd_max - mvd_min) / mvd_step) #Define length of y-axis as number of mvd points
    mtad_length = 1 + math.floor((mtad_max - mtad_min) / mtad_step) #Define length of x-axis as number of mtad points
    #Initialize correctly sized zero arrays for cross sections, detector efficiencies, dark tau odd width/mass ratios
    all_cross_sections = np.zeros((mvd_length,mtad_length), float)
    all_low_efficiencies = np.zeros((mvd_length,mtad_length), float)
    all_high_efficiencies = np.zeros((mvd_length,mtad_length), float)
    all_tau_width_mass_ratios = np.zeros((mvd_length,mtad_length), float)

    #Retrieve data from save file, gives '[]\n' if nothing saved
    saved_cross_sections = open(saved_scan_data_path).readlines()[0]
    saved_low_efficiencies = open(saved_scan_data_path).readlines()[3]
    saved_high_efficiencies = open(saved_scan_data_path).readlines()[4]
    saved_width_mass_ratios = open(saved_scan_data_path).readlines()[5]
    cross_saved = False
    if not saved_cross_sections == '[]\n': #Check if the cross sections were already retrieved and saved to file
        cross_saved = True
        print("Cross sections are saved, use old ones")
        #Modify the saved string to build the numpy array
        saved_cross_sections = saved_cross_sections.replace('] [', '][') #Remove spaces
        saved_cross_sections = saved_cross_sections.replace('[[', '[') #Remove double brackets at start
        saved_cross_sections = saved_cross_sections.replace(']]', ']') #Remove double brackets at end
        saved_cross_sections = saved_cross_sections.split(']') #Create columns by spliting the row
        saved_cross_sections.pop()
        for i in range(len(saved_cross_sections)):
            line = saved_cross_sections[i] #Retrieve a row with cross sections
            line = line.replace('[', '') #Remove the left bracket from the row
            line = line.replace('0.     ', '0. ') #Remove double spaces between rows
            line = line.replace(' ]', ']') # Remove space after final element
            line = line.split(' ') #Data in string separated by spaces
            while(line[-1] == ''):
                line.pop() #Remove all white spaces at the end of lines
            for j in range(len(line)):
                cross_section = line[j]
                cross_section = float(cross_section) #Convert from string-number to float-number
                all_cross_sections[i][j] = cross_section #Add cross section to correct place in cross section 2d array

    low_eff_saved = False
    if not saved_low_efficiencies == '[]\n': #Check if the cross sections were already retrieved and saved to file
        print("Low efficiencies are saved, use old ones")
        low_eff_saved = True
        saved_low_efficiencies = saved_low_efficiencies.replace('] [', '][') #Modify the saved string to recreate the numpy array
        saved_low_efficiencies = saved_low_efficiencies.replace('[[', '[') #Remove double brackets at start
        saved_low_efficiencies = saved_low_efficiencies.replace(']]', ']') #Remove double brackets at end
        saved_low_efficiencies = saved_low_efficiencies.split(']') #Create columns by spliting the rows
        saved_low_efficiencies.pop() #Removes final linebreak (last element after split())
        for i in range(len(saved_low_efficiencies)):
            line = saved_low_efficiencies[i] #Retrieve a row with low efficiencies
            line = line.replace('[', '') #Remove the left bracket from the row
            line = line.replace('  ', ' ') #Remove double spaces between rows
            line = line.split(' ') #Data in string separated by spaces
            while(line[-1] == ''):
                line.pop() #Remove all white spaces at the end of lines
            for j in range(len(line)):
                low_eff = line[j]
                low_eff = float(low_eff) #Convert from string-number to float
                all_low_efficiencies[i][j] = low_eff #Add low efficiency to correct place in low efficiency 2d array

    high_eff_saved = False
    if not saved_high_efficiencies == '[]\n': #Check if the cross sections were already retrieved and saved to file
        print("High efficiencies are saved, use old ones")
        high_eff_saved = True
        saved_high_efficiencies = saved_high_efficiencies.replace('] [', '][') #Modify the saved string to recreate the numpy array
        saved_high_efficiencies = saved_high_efficiencies.replace('[[', '[') #Remove double brackets at start
        saved_high_efficiencies = saved_high_efficiencies.replace(']]', ']') #Remove double brackets at end
        saved_high_efficiencies = saved_high_efficiencies.split(']') #Removes final linebreak (last element after split())
        saved_high_efficiencies.pop()
        for i in range(len(saved_high_efficiencies)):
            line = saved_high_efficiencies[i] #Retrieve a row with high efficiencies
            line = line.replace('[', '') #Remove the left bracket from the row
            line = line.replace('  ', ' ') #Remove double spaces between rows
            line = line.split(' ') #Data in string separated by spaces
            while(line[-1] == ''):
                line.pop() #Remove all white spaces at the end of lines
            for j in range(len(line)):
                high_eff = line[j]
                high_eff = float(high_eff) #Convert from string-number to float
                all_high_efficiencies[i][j] = high_eff #Add high efficiency to correct place in high efficiency 2d array

    width_mass_ratios_saved = False
    if not saved_width_mass_ratios == '[]\n': #Check if the cross sections were already retrieved and saved to file
        print("Width/Mass ratios are saved, use old ones")
        width_mass_ratios_saved = True
        saved_width_mass_ratios = saved_width_mass_ratios.replace('] [', '][') #Modify the saved string to recreate the numpy array
        saved_width_mass_ratios = saved_width_mass_ratios.replace('[[', '[') #Remove double brackets at start
        saved_width_mass_ratios = saved_width_mass_ratios.replace(']]', ']') #Remove double brackets at end
        saved_width_mass_ratios = saved_width_mass_ratios.split(']') #Removes final linebreak (last element after split())
        saved_width_mass_ratios.pop()
        for i in range(len(saved_width_mass_ratios)):
            line = saved_width_mass_ratios[i] #Retrieve a row with width/mass ratios
            line = line.replace('[', '') #Remove the left bracket from the row
            line = line.replace('  ', ' ') #Remove double spaces between rows
            line = line.replace(' ]', ']') # Remove space after final element
            line = line.split(' ') #Data in string separated by spaces
            for j in range(len(line)):
                width_mass_ratio = line[j]
                width_mass_ratio = float(width_mass_ratio) #Convert from string-number to float
                all_tau_width_mass_ratios[i][j] = width_mass_ratio #Add width/mass to correct place in width/mass 2d array

    for root, subdirectories, files in os.walk(scan_path):
        #Initialize cross section and mass variables
        cross_section = 0
        mvd_index = 0
        mtad_index = 0
        mvd = 0
        mtad = 0
        for file in files:
            if file == 'saved_data': #Find the data file to see which point (mvd, mtad) the loop is at
                save_file_path = os.path.join(root, file) # Attach full file path to '/saved_data'
                point = read_save_file(save_file_path, 0) # point = [mvd, mtad]
                mvd = point[0]
                mtad = point[1]
                mvd_index = math.floor((mvd - mvd_min) / mvd_step)
                mtad_index = math.floor((mtad - mtad_min) / mtad_step)

        if cross_saved == False: #Cross sections were not saved, retrieve from lhe-file
            for file in files:
                if file == 'unweighted_events.lhe': #Find the event file, which should be unzipped after run_simulation()
                    event_file_path = os.path.join(root, file) #Attach point directory path to file name
                    event_file = lhe_parser.EventFile(event_file_path)
                    sum_of_weights = 0
                    number_of_events = 0
                    print('mvd: ', mvd)
                    print('mtad: ', mtad)
                    for event in event_file: #Cross section as average of weights of all events
                        sum_of_weights += event.wgt
                        number_of_events += 1
                    cross_section = sum_of_weights / number_of_events
                    print('Cross section: ', cross_section)
                    all_cross_sections[mvd_index][mtad_index] = cross_section # 0 is inserted if there was no event file created by MadGraph (for not allowed, ignored or crashed points)

        if low_eff_saved == False: #Low efficiencies were not saved, retrieve from CLs file (madanalysis output)
            for file in files:
                if file == 'CLs_output_summary.dat': #Find the CLs file in point directory
                    CLs_file_path = os.path.join(root, file) #Attach point directory path to file name
                    data = open(CLs_file_path, "r").readlines()
                    first_line = data[1] #Containing low efficiencies
                    first_line = first_line[first_line.find('||'):] #Remove everything before efficiencies
                    first_line = first_line[2:] #Remove spaces
                    low_eff = float(first_line[0:9]) #Effiencies always given with same number of digits
                    all_low_efficiencies[mvd_index][mtad_index] = low_eff

        if high_eff_saved == False: #High efficiencies were not saved, retrieve from CLs file (madanalysis output)
            for file in files:
                if file == 'CLs_output_summary.dat': #Find the CLs file in point directory
                    CLs_file_path = os.path.join(root, file) #Attach point directory path to file name
                    data = open(CLs_file_path, "r").readlines()
                    second_line = data[2] #Containing high efficiencies
                    second_line = second_line[second_line.find('||'):] #Remove everything before efficiencies
                    second_line = second_line[2:] #Remove spaces
                    high_eff = float(second_line[0:9]) #Effiencies always given with same number of digits
                    all_high_efficiencies[mvd_index][mtad_index] = high_eff

        if width_mass_ratios_saved == False:
            for file in files:
                if file == 'mvd' + str(math.floor(mvd)) + '_mtad' + str(math.floor(mtad)) + '_tag_1_banner.txt': #Find the CLs file in point directory
                    banner_file_path = os.path.join(root, file) #Attach point directory path to file name
                    banner_file = open(banner_file_path, 'r')
                    tad_width_line = banner_file.readlines()[616] #Retrieve the line containing the dark tau odd width from the banner file
                    tad_width = tad_width_line[17:29] #Retrieve the value for dark tau odd width from the line
                    width_mass_ratio = float(tad_width) / mtad
                    all_tau_width_mass_ratios[mvd_index][mtad_index] = width_mass_ratio # 0 is inserted if there was no event file created by MadGraph (for not allowed, ignored or crashed points)

    #Save cross sections to save file
    numpy_string = np.array2string(all_cross_sections).replace(']\n', ']')
    numpy_string = numpy_string.replace('\n', '')
    numpy_string = numpy_string.replace('  ', ' ')
    edit_line(scan_path + '/saved_scan_data', 0, numpy_string)

    #Save low effiencies to save file
    numpy_string = np.array2string(all_low_efficiencies).replace(']\n', ']')
    numpy_string = numpy_string.replace('\n', '')
    numpy_string = numpy_string.replace('  ', ' ')
    edit_line(scan_path + '/saved_scan_data', 3, numpy_string)

    #Save high effiencies to save file
    numpy_string = np.array2string(all_high_efficiencies).replace(']\n', ']')
    numpy_string = numpy_string.replace('\n', '')
    numpy_string = numpy_string.replace('  ', ' ')
    edit_line(scan_path + '/saved_scan_data', 4, numpy_string) # Save the created numpy array as a numpy-array-string on one line in saved_scan_data

    #Save dark tau odd width/mass ratios to save file
    numpy_string = np.array2string(all_tau_width_mass_ratios).replace(']\n', ']')
    numpy_string = numpy_string.replace('\n', '')
    numpy_string = numpy_string.replace('  ', ' ')
    edit_line(scan_path + '/saved_scan_data', 5, numpy_string)

    return [all_cross_sections, all_low_efficiencies, all_high_efficiencies, all_tau_width_mass_ratios] #Data stored in a 2D numpy arrays

def exclusion_discovery_significance_contours(scan_paths, luminosity):
    significance_set = []
    tau_width_mass_ratios = []
    relic_densities = []
    dark_couplings = []
    X_axises = []
    Y_axises = []
    for i in range(len(scan_paths)):
        dark_couplings.append(read_save_file(scan_paths[i] + '/saved_scan_data', 7)[0])
        mvd_properties = read_save_file(scan_paths[i] + '/saved_scan_data', 1)
        mtad_properties = read_save_file(scan_paths[i] + '/saved_scan_data', 2)
        X_axises.append(np.array(inclusive_range(mtad_properties[0], mtad_properties[1], mtad_properties[2]))) #mtad on x-axis
        Y_axises.append(np.array(inclusive_range(mvd_properties[0], mvd_properties[1], mvd_properties[2]))) #mvd on y-axis
    for scan_path in scan_paths:
        dark_couplings.append(read_save_file(scan_paths[i] + '/saved_scan_data', 7)[0])
        cross_efficiencies = retrieve_data(scan_path)
        cross_sections = cross_efficiencies[0]
        low_efficiencies = cross_efficiencies[1]
        high_efficiencies = cross_efficiencies[2]
        tau_width_mass_ratios.append(cross_efficiencies[3])
        significance_set.append(significance_calculator(scan_path, luminosity, cross_sections, low_efficiencies, high_efficiencies, 'Zstats'))
        relic_densities.append(relic_density(scan_path, dark_couplings[i]))

    #Code supports plotting several scans in the same graph.
    #Separate the exclusion and discovery arrays from the significance output
    exclusion_significances = []
    discovery_significances = []
    for i in range(len(significance_set)):
        exclusion_significances.append(significance_set[i][1])
        discovery_significances.append(significance_set[i][2])

    print(exclusion_significances)
    print(discovery_significances)
    print('densities: ', relic_densities)
    print('ratios: ', tau_width_mass_ratios)

    discovery_levels = inclusive_range(5,5, 10) #Discovery at Z > 5
    exclusion_levels = inclusive_range(2, 2, 10) #Exclusion at Z > 2
    relic_levels = inclusive_range(0.12, 50000, 50000) # 'Everything' larger than 0.12 (Planck measurement for relic density)

    plt.figure()
    linestyles = ['solid', 'dashed', 'dashdot', 'dotted', 'solid', 'dashed', 'dashdot', 'dotted']
    relic_colors = ['purple', 'pink', 'cyan', 'indigo']
    exclusion_contours = []
    discovery_contours = []
    relic_density_regions = []
    for i in range(len(exclusion_significances)):
        linestyle = linestyles[i]
        exclusion_contours.append(plt.contour(X_axises[i], Y_axises[i], exclusion_significances[i], exclusion_levels, colors = 'green', linestyles = linestyles[i]))
        discovery_contours.append(plt.contour(X_axises[i], Y_axises[i], discovery_significances[i], discovery_levels, colors = 'blue', linestyles = linestyles[i]))
        relic_density_regions.append(plt.contourf(X_axises[i], Y_axises[i], relic_densities[i], relic_levels, colors = relic_colors[i], alpha = 0.3))
    NWA_safe = [True] * len(tau_width_mass_ratios)
    # Set NWA_safe to false if any ratios are larger than 0.01, which will plot the constant ratio lines if necessary
    for k in range(len(tau_width_mass_ratios)):
        for i in range(len(tau_width_mass_ratios[k])):
            line = tau_width_mass_ratios[k][i]
            for j in range(len(line)):
                ratio = float(line[j])
                if ratio >= 0.01:
                    NWA_safe[k] = False
                    break
    #Plot the constant ratio lines if even 'close' to breaking the NWA in any point
    if not NWA_safe:
        constant_ratios = []
        ratio_levels = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5]
        for i in range(len(constant_ratios)):
            constant_ratios.append(plt.contour(X_axises[i],Y_axises[i], tau_width_mass_ratios[i], ratio_levels, colors = 'black', alpha = 0.7, linestyles = linestyles[i]))
            plt.clabel(constant_ratios[i], inline_spacing = -12, fontsize = 14) #Label the lines with their ratios
    ax = plt.gca()
    ax.set_ylim([mvd_properties[0], mvd_properties[1]])
    ax.set_xlim([mtad_properties[0], mtad_properties[1]])
    x = np.arange(0, mtad_properties[1])
    plt.plot(x, x, 'k') # Kinematical limit is mvd = mtad
    ax.fill_between(x, 500, x, facecolor = 'red') #Fill the kinematically forbidden region with red
    plt.text(80, 205, 'm' + r'$_{W_{D \pm}}$ > m' + r'$_{\psi_D}$', fontsize = 20)
    plt.ylabel('m' + r'$_{W_{D \pm}}$ [GeV]',fontsize = 18)
    plt.xlabel('m' + r'$_{\psi_D}$ [GeV]', fontsize = 18)
    legends = []
    for i in range(len(exclusion_contours)):
        legends.append(mlines.Line2D([], [], color='blue', linestyle = linestyles[i], label='Z > 5, g' + r'$_D$ = ' +str(dark_couplings[i]), linewidth=3))
        legends.append(mlines.Line2D([], [], color='green', linestyle = linestyles[i], label='Z > 2, g' + r'$_D$ = ' +str(dark_couplings[i]), linewidth=3))
    ax.legend(handles = legends, loc = 'upper right', fontsize = 20) #loc = (0.43, 0.85)
    plt.title('Exclusion and discovery contours; No scalar mixing, m' + r'$_\psi$=800GeV, m' + r'$_{h2}$=300GeV, L='+str(luminosity)+'fb' + r'$^{-1}$', fontsize = 18) #g' + r'$_D$=0.01'
    plt.show()
