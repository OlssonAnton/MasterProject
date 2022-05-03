import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from general_stuff import inclusive_range, z
from significance_script import significance_calculator

def create_significance_arrays(luminosity, signals, mvd, mtad, significance_type, analysis_paths, event_paths):
    mvd_arr = np.array(mvd)
    mtad_arr = np.array(mtad)
    X, Y = np.meshgrid(mtad_arr, mvd_arr)

    mvd_min = mvd[0]
    mtad_min = mtad[0]
    step = mvd[1] - mvd_min
    counter = 0
    while mvd_min < mtad_min: #Setup for counter
        counter += 1
        mvd_min += step

    significances = []
    discoveries = []
    exclusions = []

    for i in range(len(signals)):
        discovery_matrix = []
        exclusion_matrix = []
        significance_matrix = []

        if significance_type == 'Zstats':
            #Asymptotic significances not generated in this case
            significance_lists = significance_calculator(analysis_paths[i], signals[i], event_paths[i], 'Zstats')
            discovery_list = significance_lists[1]
            exclusion_list = significance_lists[2]
            current_counter = counter
            while(discovery_list):
                discovery_matrix.append(discovery_list[:current_counter])
                del discovery_list[:current_counter]
                exclusion_matrix.append(exclusion_list[:current_counter])
                del exclusion_list[:current_counter]
                current_counter += 1
            for i in range(len(discovery_matrix)):
                missing_length = len(mvd) - len(discovery_matrix[i])
                discovery_matrix[i] += [0] * missing_length
                exclusion_matrix[i] += [0] * missing_length

        if significance_type == 'asymptotic limit':
            #Zstats significances not generated in this case
            all_significance_lists = significance_calculator(analysis_paths[i], signals[i], event_paths[i], 'asymptotic limit')
            significance_list = all_significance_lists[0]
            current_counter = counter
            while(significance_list):
                significance_matrix.append(significance_list[:current_counter])
                del significance_list[:current_counter]
                current_counter += 1
            for i in range(len(significance_matrix)):
                missing_length = len(mvd) - len(significance_matrix[i])
                significance_matrix[i] += [0] * missing_length

        #Make arrays into numpy arrays for plotting
        significance_array = np.array(significance_matrix)
        exclusion_array = np.array(exclusion_matrix)
        discovery_array = np.array(discovery_matrix)

        #Putting DM candidate on y-axis, so z-values have to be transposed
        significance_array = np.transpose(significance_array)
        exclusion_array = np.transpose(exclusion_array)
        discovery_array = np.transpose(discovery_array)

        #If plotting several parameter-sets in one graph
        significances.append(significance_array)
        discoveries.append(discovery_array)
        exclusions.append(exclusion_array)

    return [X,Y, significances, exclusions, discoveries]

def asymptotic_limit_significance_contour(luminosity, signals, mvd, mtad, analysis_paths, event_paths):
    data = create_significance_arrays(luminosity, signals, mvd, mtad, 'asymptotic limit', analysis_paths, event_paths)
    X = data[0]
    Y = data[1]
    gD1 = data[2][0]
    gD3 = data[2][1]
    x = np.arange(0, 251)
    ax = plt.gca()
    ax.set_ylim([mvd[0], mvd[-1]])
    ax.set_xlim([mtad[0], mtad[-1]])
    ax.fill_between(x, 300, x, facecolor = 'red')
    plt.plot(x, x, 'k')
    plt.text(145, 205, 'mvd > mtad', fontsize = 20, )
    plt.ylabel('Dark Matter Candidate Mass (mvd) GeV')
    plt.xlabel('Tau Odd Mass (mtad) GeV')
    exclusion_levels = inclusive_range(2,2,1)
    discovery_levels = inclusive_range(5,5,1)
    gD1_exclusion_contour = plt.contour(X, Y, gD1, exclusion_levels, colors = 'green')
    gD1_discovery_contour = plt.contour(X, Y, gD1, discovery_levels, colors = 'blue')
    gD3_exclusion_contour = plt.contour(X, Y, gD3, exclusion_levels, colors = 'green', linestyles = 'dashed')
    gD3_discovery_contour = plt.contour(X, Y, gD3, discovery_levels, colors = 'blue', linestyles = 'dashed')
    discovery_gD1 = mlines.Line2D([], [], color='blue', label='Z > 5, gD = 1 (NWA)')
    discovery_gD3 = mlines.Line2D([], [], color='blue', label='Z > 5, gD = 3', linestyle = 'dashed')
    exclusion_gD1 = mlines.Line2D([], [], color='green', label='Z > 2, gD = 1 (NWA)')
    exclusion_gD3 = mlines.Line2D([], [], color='green', label='Z > 2, gD = 3', linestyle = 'dashed')
    ax.legend(handles=[discovery_gD1, discovery_gD3, exclusion_gD1, exclusion_gD3], loc = 'lower right')
    plt.title('Asymptotic limit significance: Theta_s = 0, mtap=800GeV, mh2=300GeV, L=139fb^-1')
    plt.show()

def asymptotic_limit_significance_color(luminosity, signals, background, mvd, mtad, analysis_paths, event_paths):
    data = create_significance_arrays(luminosity, signals, background, mvd, mtad, 'asymptotic limit', analysis_paths, event_paths)
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

def exclusion_discovery_significance_contours(luminosity, signals, mvd, mtad, analysis_paths, event_paths):
    data = create_significance_arrays(luminosity, signals, mvd, mtad, 'Zstats', analysis_paths, event_paths)
    X = data[0]
    Y = data[1]
    gD1_exclusion = data[3][0]
    gD3_exclusion = data[3][1]
    gD1_discovery = data[4][0]
    gD3_discovery = data[4][1]
    vf = [0, 250]
    x = np.arange(0, 251)
    discovery_levels = inclusive_range(5,10,10)
    exclusion_levels = inclusive_range(2, 10, 10)
    plt.figure()
    gD1_exclusion_contour = plt.contour(X, Y, gD1_exclusion, exclusion_levels, colors = 'green')
    gD1_discovery_contour = plt.contour(X, Y, gD1_discovery, discovery_levels, colors = 'blue')
    gD3_exclusion_contour = plt.contour(X, Y, gD3_exclusion, exclusion_levels, colors = 'green', linestyles = 'dashed')
    gD3_discovery_contour = plt.contour(X, Y, gD3_discovery, discovery_levels, colors = 'blue', linestyles = 'dashed')
    ax = plt.gca()
    ax.set_ylim([mvd[0], mvd[-1]])
    ax.set_xlim([mtad[0], mtad[-1]])
    ax.fill_between(x, 300, x, facecolor = 'red')
    plt.plot(x, x, 'k')
    plt.text(145, 205, 'mvd > mtad', fontsize = 20, )
    plt.ylabel('Dark Matter Candidate Mass (mvd) GeV')
    plt.xlabel('Tau Odd Mass (mtad) GeV')
    discovery_gD1 = mlines.Line2D([], [], color='blue', label='Z > 5, gD = 1 (NWA)')
    discovery_gD3 = mlines.Line2D([], [], color='blue', label='Z > 5, gD = 3', linestyle = 'dashed')
    exclusion_gD1 = mlines.Line2D([], [], color='green', label='Z > 2, gD = 1 (NWA)')
    exclusion_gD3 = mlines.Line2D([], [], color='green', label='Z > 2, gD = 3', linestyle = 'dashed')
    ax.legend(handles=[discovery_gD1, discovery_gD3, exclusion_gD1, exclusion_gD3], loc = 'lower right')
    plt.title('Zstats significance: Theta_s = 0, mtap=800GeV, mh2=300GeV, L=139fb^-1')
    plt.show()
