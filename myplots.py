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

def retrieve_cross_sections_efficiencies(scan_path):
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
    all_cross_sections = np.zeros((mvd_length,mtad_length), float) #Create empty array to insert cross sections from file
    all_low_efficiencies = np.zeros((mvd_length,mtad_length), float) #Create empty array to insert efficiencies from file
    all_high_efficiencies = np.zeros((mvd_length,mtad_length), float)

    saved_cross_sections = open(saved_scan_data_path).readlines()[0]
    saved_low_efficiencies = open(saved_scan_data_path).readlines()[3]
    saved_high_efficiencies = open(saved_scan_data_path).readlines()[4]
    cross_saved = False
    if not saved_cross_sections == '[]\n': #Check if the cross sections were already retrieved and saved to file
        cross_saved = True
        print("Cross sections are saved, use old ones")
        saved_cross_sections = saved_cross_sections.replace('] [', '][') #Modify the saved string to recreate the numpy array
        saved_cross_sections = saved_cross_sections.replace('[[', '[')
        saved_cross_sections = saved_cross_sections.replace(']]', ']')
        saved_cross_sections = saved_cross_sections.split(']')
        saved_cross_sections.pop()
        for i in range(len(saved_cross_sections)):
            line = saved_cross_sections[i]
            line = line.replace('[', '') #Final modifications of saved numpy-array-string
            line = line.split(' ')
            for j in range(len(line)):
                cross_section = line[j]
                #print(cross_section)
                cross_section = float(cross_section) #Convert from string-number to float
                all_cross_sections[i][j] = cross_section

    low_eff_saved = False
    if not saved_low_efficiencies == '[]\n': #Check if the cross sections were already retrieved and saved to file
        print("Low efficiencies are saved, use old ones")
        low_eff_saved = True
        saved_low_efficiencies = saved_low_efficiencies.replace('] [', '][') #Modify the saved string to recreate the numpy array
        saved_low_efficiencies = saved_low_efficiencies.replace('[[', '[')
        saved_low_efficiencies = saved_low_efficiencies.replace(']]', ']')
        saved_low_efficiencies = saved_low_efficiencies.split(']')
        saved_low_efficiencies.pop() #Remove linebreak at the end
        for i in range(len(saved_low_efficiencies)):
            line = saved_low_efficiencies[i]
            line = line.replace('[', '') #Final modifications of saved numpy-array-string
            line = line.replace('  ', ' ')
            line = line.split(' ')
            line.pop()
            print('line: ', line)
            for j in range(len(line)):
                low_eff = line[j]
                print('low eff: ', low_eff)
                low_eff = float(low_eff) #Convert from string-number to float
                all_low_efficiencies[i][j] = low_eff

    high_eff_saved = False
    if not saved_high_efficiencies == '[]\n': #Check if the cross sections were already retrieved and saved to file
        print("High efficiencies are saved, use old ones")
        high_eff_saved = True
        saved_high_efficiencies = saved_high_efficiencies.replace('] [', '][') #Modify the saved string to recreate the numpy array
        saved_high_efficiencies = saved_high_efficiencies.replace('[[', '[')
        saved_high_efficiencies = saved_high_efficiencies.replace(']]', ']')
        saved_high_efficiencies = saved_high_efficiencies.split(']')
        saved_high_efficiencies.pop()
        for i in range(len(saved_high_efficiencies)):
            line = saved_high_efficiencies[i]
            line = line.replace('[', '') #Final modifications of saved numpy-array-string
            line = line.replace('  ', ' ')
            line = line.split(' ')
            line.pop()
            for j in range(len(line)):
                high_eff = line[j]
                high_eff = float(high_eff) #Convert from string-number to float
                all_high_efficiencies[i][j] = high_eff

    for root, subdirectories, files in os.walk(scan_path):
        cross_section = 0 #0 = Default value for ignored, crashed or not allowed points
        mvd_index = 0
        mtad_index = 0
        for file in files:
            if file == 'saved_data': #Find the data file to see which point (mvd, mtad) the loop is at
                save_file_path = os.path.join(root, file) # Attach full file path to '/saved_data'
                point = read_save_file(save_file_path, 3) # point = [mvd, mtad]
                mvd = point[0]
                mtad = point[1]
                mvd_index = math.floor(mvd / mvd_step)
                mtad_index = math.floor(mtad / mtad_step)

        if cross_saved == False: #Cross sections were not saved, retrieve from lhe-file
            for file in files:
                if file == 'unweighted_events.lhe': #Find the event file, which should be unzipped after run_simulation()
                    event_file_path = os.path.join(root, file) #Attach point directory path to file name
                    event_file = lhe_parser.EventFile(event_file_path)
                    sum_of_weights = 0
                    number_of_events = 0
                    for event in event_file: #Cross section as average of all weights
                        sum_of_weights += event.wgt
                        number_of_events += 1
                    cross_section = sum_of_weights / number_of_events
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


    numpy_string = np.array2string(all_cross_sections).replace(']\n', ']')
    numpy_string = numpy_string.replace('\n', '')
    numpy_string = numpy_string.replace('  ', ' ')
    edit_line(scan_path + '/saved_scan_data', 0, numpy_string)

    numpy_string = np.array2string(all_low_efficiencies).replace(']\n', ']')
    numpy_string = numpy_string.replace('\n', '')
    numpy_string = numpy_string.replace('  ', ' ')
    edit_line(scan_path + '/saved_scan_data', 3, numpy_string)

    numpy_string = np.array2string(all_high_efficiencies).replace(']\n', ']')
    numpy_string = numpy_string.replace('\n', '')
    numpy_string = numpy_string.replace('  ', ' ')
    edit_line(scan_path + '/saved_scan_data', 4, numpy_string) # Save the created numpy array as a numpy-array-string on one line in saved_scan_data

    return [all_cross_sections, all_low_efficiencies, all_high_efficiencies] #Cross sections stored in a 2D numpy array

def asimov_significance_contour(scan_paths, luminosity):
    significance_set = []
    for scan_path in scan_paths:
        cross_efficencies = retrieve_cross_sections_efficiencies(scan_path)
        cross_sections = cross_efficencies[0]
        low_efficiencies = cross_efficencies[1]
        high_efficiencies = cross_efficencies[2]
        significance_set.append(significance_calculator(scan_path, luminosity, cross_sections, low_efficiencies, high_efficiencies, 'asimov'))


    mvd_properties = read_save_file(scan_paths[0] + '/saved_scan_data', 1)
    mtad_properties = read_save_file(scan_paths[0] + '/saved_scan_data', 2)
    X = np.array(inclusive_range(mtad_properties[0], mtad_properties[1], mtad_properties[2]))
    Y = np.array(inclusive_range(mvd_properties[0], mvd_properties[1], mvd_properties[2]))

    gD1exclusion = significance_set[0][0]
    gD1discovery = significance_set[0][0]
    x = np.arange(0, 500)
    ax = plt.gca()
    # ax.set_ylim([mvd[0], mvd[-1]])
    # ax.set_xlim([mtad[0], mtad[-1]])
    ax.set_ylim(0, 500)
    ax.set_xlim(0, 800)
    ax.fill_between(x, 500, x, facecolor = 'red')
    plt.plot(x, x, 'k')
    plt.text(145, 205, 'mvd > mtad', fontsize = 20, )
    plt.ylabel('Dark Matter Candidate Mass (mvd) GeV')
    plt.xlabel('Tau Odd Mass (mtad) GeV')
    exclusion_levels = inclusive_range(2,2,1)
    discovery_levels = inclusive_range(5,5,1)
    gD1_exclusion_contour = plt.contour(X, Y, gD1exclusion, exclusion_levels, colors = 'green')
    gD1_discovery_contour = plt.contour(X, Y, gD1discovery, discovery_levels, colors = 'blue')
    #gD3_exclusion_contour = plt.contour(X, Y, gD3, exclusion_levels, colors = 'green', linestyles = 'dashed')
    #gD3_discovery_contour = plt.contour(X, Y, gD3, discovery_levels, colors = 'blue', linestyles = 'dashed')
    discovery_gD1 = mlines.Line2D([], [], color='blue', label='Z > 5, gD = 1 (NWA)')
    #discovery_gD3 = mlines.Line2D([], [], color='blue', label='Z > 5, gD = 3', linestyle = 'dashed')
    exclusion_gD1 = mlines.Line2D([], [], color='green', label='Z > 2, gD = 1 (NWA)')
    #exclusion_gD3 = mlines.Line2D([], [], color='green', label='Z > 2, gD = 3', linestyle = 'dashed')
    ax.legend(handles=[discovery_gD1, exclusion_gD1], loc = 'lower right')
    plt.title('Asimov significance: Theta_s = 0, mtap=800GeV, mh2=300GeV, L=139fb^-1, gD=1')
    plt.show()

def asimov_significance_color(scan_paths, luminosity):
    data = create_significance_arrays(luminosity, signals, background, mvd, mtad, 'asimov', analysis_paths, event_paths)
    X = data[0]
    Y = data[1]
    Z = data[2][0]
    x = np.arange(0, 251)
    levels = inclusive_range(0, 6, 0.25)
    plt.contourf(X, Y, Z, levels) #extend='neither')
    ax = plt.gca()
    ax.set_ylim([mvd[0], mvd[-1]])
    ax.set_xlim([mtad[0], mtad[-1]])
    ax.fill_between(x, 300, x, facecolor = 'red')
    plt.plot(x, x, 'k')
    plt.text(145, 205, 'mvd > mtad', fontsize = 20, )
    plt.ylabel('Dark Matter Candidate Mass (mvd) GeV')
    plt.xlabel('Tau Odd Mass (mtad) GeV')
    plt.colorbar().set_label("Significance")  #"Cross Section (pb)"
    plt.show()

def exclusion_discovery_significance_contours(scan_paths, luminosity):
    significance_set = []
    for scan_path in scan_paths:
        cross_efficencies = retrieve_cross_sections_efficiencies(scan_path)
        cross_sections = cross_efficencies[0]
        low_efficiencies = cross_efficencies[1]
        high_efficiencies = cross_efficencies[2]
        significance_set.append(significance_calculator(scan_path, luminosity, cross_sections, low_efficiencies, high_efficiencies, 'Zstats'))

    data = open(scan_paths[0] + '/saved_scan_data', 'r')
    mvd_properties = data[1]
    mtad_properties = data[2]
    X = np.array(inclusive_range(mtad_properties[0], mtad_properties[1], mtad_properties[2]))
    Y = np.array(inclusive_range(mvd_properties[0], mvd_properties[1], mvd_properties[2]))

    gD1exclusion = significance_set[0][1]
    gD1discovery = significance_set[0][2]
    x = np.arange(0, 500)
    discovery_levels = inclusive_range(5,10,10)
    exclusion_levels = inclusive_range(2, 10, 10)
    plt.figure()
    # gD1_exclusion_contour = plt.contour(X, Y, gD1_exclusion, exclusion_levels, colors = 'green')
    # gD1_discovery_contour = plt.contour(X, Y, gD1_discovery, discovery_levels, colors = 'blue')
    # gD3_exclusion_contour = plt.contour(X, Y, gD3_exclusion, exclusion_levels, colors = 'green', linestyles = 'dashed')
    # gD3_discovery_contour = plt.contour(X, Y, gD3_discovery, discovery_levels, colors = 'blue', linestyles = 'dashed')
    gD_001_exclusion_contour = plt.contour(X, Y, gD1exclusion, exclusion_levels, colors = 'green')
    gD_001_discovery_contour = plt.contour(X, Y, gD1discovery, discovery_levels, colors = 'blue')
    ax = plt.gca()
    ax.set_ylim([mvd[0], mvd[-1]])
    ax.set_xlim([mtad[0], mtad[-1]])
    ax.fill_between(x, 500, x, facecolor = 'red')
    plt.plot(x, x, 'k')
    plt.text(100, 305, 'mvd > mtad', fontsize = 20, )
    plt.ylabel('Dark Matter Candidate Mass (mvd) GeV')
    plt.xlabel('Tau Odd Mass (mtad) GeV')
    # discovery_gD1 = mlines.Line2D([], [], color='blue', label='Z > 5, gD = 1 (NWA)')
    # discovery_gD3 = mlines.Line2D([], [], color='blue', label='Z > 5, gD = 3', linestyle = 'dashed')
    # exclusion_gD1 = mlines.Line2D([], [], color='green', label='Z > 2, gD = 1 (NWA)')
    # exclusion_gD3 = mlines.Line2D([], [], color='green', label='Z > 2, gD = 3', linestyle = 'dashed')
    gD1exclusion_legend = mlines.Line2D([], [], color='blue', label='Z > 5, gD = 0.01 (NWA)')
    gD1discovery_legend = mlines.Line2D([], [], color='green', label='Z > 2, gD = 0.01 (NWA)')
    #ax.legend(handles=[discovery_gD1, discovery_gD3, exclusion_gD1, exclusion_gD3], loc = 'lower right')
    ax.legend(handles=[gD1exclusion_legend, gD1discovery_legend], loc = 'lower right')
    plt.title('Zstats significance: Theta_s = 0, mtap=800GeV, mh2=300GeV, L=139fb^-1, gD=1')
    plt.show()
