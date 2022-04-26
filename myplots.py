import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from general_stuff import inclusive_range, z
from significance_script import significance_calculator

def z_contour_plot(luminosity, sigma_s, sigma_b, mvd, mtad):
    mvd_arr = np.array(mvd)
    mtad_arr = np.array(mtad)
    X, Y = np.meshgrid(mtad_arr, mvd_arr)

    # print("X: ", X)
    # print("Y: ", Y)

    mvd_min = mvd[0]
    mtad_min = mtad[0]
    step = mvd[1] - mvd_min
    counter = 0
    while mvd_min < mtad_min: #Setup for counter
        counter += 1
        mvd_min += step

    current_counter = counter
    cross_section_list = sigma_s
    cross_section_matrix = []
    while(cross_section_list):
        cross_section_matrix.append(cross_section_list[:current_counter])
        del cross_section_list[:current_counter]
        current_counter += 1

    #significance_lists = significance_calculator("../../HEPTools/madanalysis5/madanalysis5/restructured_careful_scan2")
    significance_lists = significance_calculator('./gD1')
    discovery_list = significance_lists[0]
    exclusion_list = significance_lists[1]
    discovery_matrix = []
    exclusion_matrix = []
    current_counter = counter
    while(discovery_list):
        discovery_matrix.append(discovery_list[:current_counter])
        del discovery_list[:current_counter]
        exclusion_matrix.append(exclusion_list[:current_counter])
        del exclusion_list[:current_counter]
        current_counter += 1

    #fill out forbidden regions with zeros
    for i in range(len(cross_section_matrix)):
        missing_length = len(mvd) - len(cross_section_matrix[i])
        cross_section_matrix[i] += [0] * missing_length
        discovery_matrix[i] += [0] * missing_length
        exclusion_matrix[i] += [0] * missing_length

    print("DISCOVERY MATRIX: ", discovery_matrix)
    print("EXCLUSION MATRIX: ", exclusion_matrix)
    print("DISCOVERY LENGTH: ", len(discovery_matrix))
    print("EXCLUSION LENGTH: ", len(exclusion_matrix))

    significance_matrix = []
    for line in cross_section_matrix:
        significance_line = []
        for cross_section in line:
            significance_line.append(z(luminosity, cross_section, sigma_b))
        significance_matrix.append(significance_line)

    #print("CROSS SECTION MATRIX: ", cross_section_matrix)

    significance_array = np.array(significance_matrix)
    cross_section_array = np.array(cross_section_matrix)
    exclusion_array = np.array(exclusion_matrix)
    discovery_array = np.array(discovery_matrix)

    #Axises are switched
    significance_array = np.transpose(significance_array)
    cross_section_array = np.transpose(cross_section_array)
    exclusion_array = np.transpose(exclusion_array)
    discovery_array = np.transpose(discovery_array)

    vf = [0, 250]
    x = np.arange(0, 251)
    print("significance_array: ", significance_array)
    levels = inclusive_range(2, 5, 3)
    #plt.contour(X, Y, significance_array, levels, extend='neither', colors = 'green')
    discovery_levels = inclusive_range(4,5,1)
    exclusion_levels = inclusive_range(2, 5, 3)
    plt.figure()
    exclusion_contour = plt.contour(X, Y, significance_array, exclusion_levels, colors = 'green')
    discovery_contour = plt.contour(X, Y, significance_array, discovery_levels, colors = 'blue')
    ax = plt.gca()
    ax.set_ylim([mvd[0], mvd[-1]])
    ax.set_xlim([mtad[0], mtad[-1]])
    ax.fill_between(x, 300, x, facecolor = 'red')
    plt.plot(x, x, 'k')
    plt.text(145, 205, 'mvd > mtad', fontsize = 20, )
    plt.ylabel('Dark Matter Candidate Mass (mvd) GeV')
    plt.xlabel('Tau Odd Mass (mtad) GeV')
    discovery_gD1 = mlines.Line2D([], [], color='green', label='Z > 5, gD = 1 (NWA)')
    discovery_gD3 = mlines.Line2D([], [], color='green', label='Z > 5, gD = 3', linestyle = 'dashed')
    exclusion_gD1 = mlines.Line2D([], [], color='blue', label='Z > 2, gD = 1 (NWA)')
    exclusion_gD3 = mlines.Line2D([], [], color='blue', label='Z > 2, gD = 3', linestyle = 'dashed')
    ax.legend(handles=[discovery_gD1, discovery_gD3, exclusion_gD1, exclusion_gD3], loc = 'lower right')
    #plt.colorbar().set_label("Significance")  #"Cross Section (pb)"
    plt.title('All Diagrams: g_D = 0.5, theta_s = 0, m_tap = 800GeV, M_h2 = 300GeV, L=139fb^-1')
    plt.show()

def z_contour_plot_old(luminosity, sigma_s, sigma_b, mvd, mtad):
    mvd_arr = np.array(mvd)
    mtad_arr = np.array(mtad)
    X, Y = np.meshgrid(mvd_arr, mtad_arr)

    print("X: ", X)
    print("Y: ", Y)

    mvd_min = mvd[0]
    mtad_min = mtad[0]
    step = mvd[1] - mvd_min
    counter = 0
    while mvd_min < mtad_min:
        counter += 1
        mvd_min += step

    cross_section_list = sigma_s
    cross_section_matrix = []
    while(cross_section_list):
        cross_section_matrix.append(cross_section_list[:counter])
        del cross_section_list[:counter]
        counter += 1

    for i in range(len(cross_section_matrix)):
        missing_length = len(mvd) - len(cross_section_matrix[i])
        cross_section_matrix[i] += [0] * missing_length

    significance_matrix = []

    for line in cross_section_matrix:
        significance_line = []
        for cross_section in line:
            significance_line.append(z(luminosity, cross_section, sigma_b))
        significance_matrix.append(significance_line)

    print("CROSS SECTION MATRIX: ", cross_section_matrix)

    significance_array = np.array(significance_matrix)
    cross_section_array = np.array(cross_section_matrix)

    vf = [0, 250]
    x = np.arange(0, 251)

    levels = inclusive_range(0, 5, 0.25)
    plt.contourf(X, Y, significance_array, levels, extend='neither')
    ax = plt.gca()
    ax.set_xlim([80, 235])
    ax.set_ylim([135, 235])
    ax.fill_between(x, 0, x, facecolor = 'red')
    plt.plot(x, x, 'k')
    plt.text(180, 150, 'mvd > mtad', fontsize = 20, )
    plt.xlabel('Dark Matter Candidate Mass (mvd) GeV')
    plt.ylabel('Tau Odd Mass (mtad) GeV')
    plt.colorbar().set_label("Significance")  #"Cross Section (pb)"
    plt.title('All Diagrams: g_D = 0.5, theta_s = 0, m_tap = 800GeV, M_h2 = 300GeV, L=300fb^-1')
    plt.show()

def z_plot2(luminosity, cross_sections, background):
    axes = plt.axes()
    axes.set_ylim([0, 5]) #Not looking at significances > 5
    axes.set_xlim([0, 510])


    mtad_list = []
    mvd_list = []
    mvd_points = 5
    mtad = 100
    while mtad <= 500:
        mtad_list += [mtad] * mvd_points
        mvd_list += [1] + np.linspace(25,mtad,mvd_points - 1).tolist()
        mvd_points += 2
        mtad += 50

    print(mtad_list)
    print(mvd_list)

    mtad = 100
    groups = (500 - 100)/50 + 1
    start_point = 0
    end_point = 5
    for i in range(math.floor(groups)):
        plt.plot(mvd_list[start_point:end_point], z_list(luminosity, cross_sections[start_point:end_point],
                background), label = "mtad =" + str(mtad))
        print(mvd_list[start_point:end_point])
        print(z_list(luminosity, cross_sections[start_point:end_point], background))
        mtad += 50
        points = end_point - start_point
        start_point = end_point
        end_point += points + 2

        print("START POINT: ", start_point)
        print("END POINT: ", end_point)

    plt.ylabel('Significance')
    plt.xlabel('DM mass (GeV)')
    plt.title('g_D = 0.5, theta_s = 0, Mh1=1.255e2, Mh2 = 3.0e2')
    plt.legend(loc="upper left")
    #plt.xscale('log')
    plt.show()

def sigma_plot(cross_sections):
    axes = plt.axes()
    #axes.set_ylim([0, 5]) #Not looking at significances > 5
    axes.set_xlim([0, 510])


    mtad_list = []
    mvd_list = []
    mvd_points = 5
    mtad = 100
    while mtad <= 500:
        mtad_list += [mtad] * mvd_points
        mvd_list += [1] + np.linspace(25,mtad,mvd_points - 1).tolist()
        mvd_points += 2
        mtad += 50

    print(mtad_list)
    print(mvd_list)

    mtad = 100
    groups = (500 - 100)/50 + 1
    start_point = 0
    end_point = 5
    for i in range(math.floor(groups)):
        plt.plot(mvd_list[start_point:end_point], cross_sections[start_point:end_point],
                label = "mtad =" + str(mtad))
        print(mvd_list[start_point:end_point])
        print(cross_sections[start_point:end_point])
        mtad += 50
        points = end_point - start_point
        start_point = end_point
        end_point += points + 2

        print("START POINT: ", start_point)
        print("END POINT: ", end_point)

    plt.ylabel('Cross Section (pb)')
    plt.xlabel('DM mass (GeV)')
    plt.title('g_D = 0.5, theta_s = 0, Mh1=1.255e2, Mh2 = 3.0e2')
    plt.legend(loc="upper left")
    #plt.xscale('log')
    plt.show()

def z_plot(luminosity, cross_sections, background):
    axes = plt.axes()
    axes.set_ylim([0, 5]) #Not looking at significances > 5
    for i in range(len(cross_sections)):
        plt.plot(luminosity, z(luminosity, cross_sections[i], background), label = 'BP' + str(i+1))
    plt.ylabel('Significance')
    plt.xlabel('Luminosity (fb^-1)')
    plt.title('g_D = 0.5, theta_s = 0, Mh1=1.255e2, Mh2 = 3.0e2')
    plt.legend(loc="upper left")
    plt.xscale('log')
    plt.show()
