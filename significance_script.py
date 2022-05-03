import os
from Zstats import Zdisc, Zexcl
from general_stuff import z

def retrieve_mad_data(analysis_path, event_path):
    high_region_mad_significances = []
    low_region_mad_significances = []
    low_region_efficiencies = []
    high_region_efficiencies = []
    points = len([name for name in os.listdir(analysis_path)])
    for i in range(points):
        current_folder = "ANALYSIS_" + str(i)
        data_file = analysis_path + "/" + current_folder + "/Output/SAF/CLs_output_summary.dat"

        reading_file = open(data_file, "r")

        first_line = open(data_file).readlines()[1]
        second_line = open(data_file).readlines()[2]

        first_line = first_line.replace(" ", "")
        first_line = first_line[first_line.find('SRlow'):]
        first_line = first_line[5:]

        second_line = second_line.replace(" ", "")
        second_line = second_line[second_line.find('SRhigh'):]
        second_line = second_line[6:]

        low_exp = 0
        high_exp = 0

        if first_line[0] == "-":
            low_exp = -1
        else:
            low_exp = float(first_line[0:9])

        if second_line[0] == "-":
            high_exp = -1
        else:
            high_exp = float(second_line[0:9])

        high_region_mad_significances.append(high_exp) #expected significances
        low_region_mad_significances.append(low_exp)

        first_line = first_line[first_line.find('||'):]
        first_line = first_line[2:]

        second_line = second_line[second_line.find('||'):]
        second_line = second_line[2:]

        low_eff = float(first_line[0:9])
        high_eff = float(second_line[0:9])

        high_region_efficiencies.append(high_eff)
        low_region_efficiencies.append(low_eff)

    save_file = open(event_path + '/saved_data')
    saved_data = save_file.readlines()
    save_file.close()
    saved_data[1] = str(low_region_efficiencies) + '\n'
    saved_data[2] = str(high_region_efficiencies) + '\n'
    new_data = ''.join(saved_data)
    save_file = open(event_path + '/saved_data', 'w')
    save_file.write(new_data)
    save_file.close()

    return [low_region_mad_significances, high_region_mad_significances, low_region_efficiencies, high_region_efficiencies]

def significance_calculator(analysis_path, cross_sections, event_path, calculator):
    SM_irreducible_background = 0.4312 + 0.02746 * 2 #2 neutrino final state processes (SM background)
    low_region_background = 6.0 #nr of events
    low_region_background_uncertainty = 1.7
    high_region_background = 10.2
    high_region_background_uncertainty = 3.3
    luminosity = 139 #fb^-1

    mad_data = retrieve_mad_data(analysis_path, event_path)
    low_region_signals = mad_data[0]
    high_region_signals = mad_data[1]
    low_region_efficiencies = mad_data[2]
    high_region_efficiencies = mad_data[3]

    significance_list = []
    discovery_list = []
    exclusion_list = []

    for i in range(len(low_region_efficiencies)):
        cross_section = cross_sections[i] * 1000 #convert pb to fb
        low_region_signal = luminosity * cross_section * low_region_efficiencies[i] #luminosity in fb^-1, cross_section in fb (converted)
        high_region_signal = luminosity * cross_section * high_region_efficiencies[i]

        if calculator == 'asymptotic limit':
            z_low = z(luminosity, low_region_signal, low_region_background)
            z_high = z(luminosity, high_region_signal, high_region_background)

            if z_low > z_high:
                significance_list.append(z_low)
            else:
                significance_list.append(z_high)

        if calculator == 'Zstats':
            z_disc_low_region_signal = Zdisc(low_region_signal, low_region_background, low_region_background_uncertainty)
            z_disc_high_region_signal = Zdisc(high_region_signal, high_region_background, high_region_background_uncertainty)

            z_excl_low_region_signal = Zexcl(low_region_signal, low_region_background, low_region_background_uncertainty)
            z_excl_high_region_signal = Zexcl(high_region_signal, high_region_background, high_region_background_uncertainty)

            if z_disc_low_region_signal > z_disc_high_region_signal:
                discovery_list.append(z_disc_low_region_signal)
            else:
                discovery_list.append(z_disc_high_region_signal)

            if z_excl_low_region_signal > z_excl_high_region_signal:
                exclusion_list.append(z_excl_low_region_signal)
            else:
                exclusion_list.append(z_excl_high_region_signal)

    return [significance_list, discovery_list, exclusion_list]

def Zstats_significance_calculator(analysis_path, cross_sections, event_path):
    SM_irreducible_background = 0.4312 + 0.02746 * 2 #2 neutrino final state processes (SM background)
    low_region_background = 6.0 #nr of events
    low_region_background_uncertainty = 1.7
    high_region_background = 10.2
    high_region_background_uncertainty = 3.3
    luminosity = 139 #fb^-1

    mad_data = retrieve_mad_data(analysis_path, event_path)
    low_region_significance = mad_data[0]
    high_region_significance = mad_data[1]
    low_region_efficiencies = mad_data[2]
    high_region_efficiencies = mad_data[3]

    discovery_list = []
    exclusion_list = []

    cross_section = cross_sections[i] * 1000 #convert pb to fb
    for i in range(len(cross_sections)):
        low_region_signal = luminosity * cross_section * low_region_efficiencies[i] #luminosity in fb^-1, cross_section in fb (converted)
        high_region_signal = luminosity * cross_section * high_region_efficiencies[i]

        z_disc_low_sig = Zdisc(low_sig, background_low, background_low_uncertainty)
        z_disc_high_sig = Zdisc(high_sig, background_high, background_high_uncertainty)

        z_excl_low_sig = Zexcl(low_sig, background_low, background_low_uncertainty)
        z_excl_high_sig = Zexcl(high_sig, background_high, background_high_uncertainty)

        if z_disc_low_sig > z_disc_high_sig:
            discovery_list.append(z_disc_low_sig)
        else:
            discovery_list.append(z_disc_high_sig)

        if z_excl_low_sig > z_excl_high_sig:
            exclusion_list.append(z_excl_low_sig)
        else:
            exclusion_list.append(z_excl_high_sig)

    return [discovery_list, exclusion_list]

def asymptotic_limit_significance_calculator(analysis_path, cross_sections, event_path):
    SM_irreducible_background = 0.4312 + 0.02746 * 2 #2 neutrino final state processes (SM background)
    low_region_background = 6.0 #nr of events
    low_region_background_uncertainty = 1.7
    high_region_background = 10.2
    high_region_background_uncertainty = 3.3
    luminosity = 139 #fb^-1

    mad_data = retrieve_mad_data(analysis_path, event_path)
    low_region_signals = mad_data[0]
    high_region_signals = mad_data[1]
    low_region_efficiencies = mad_data[2]
    high_region_efficiencies = mad_data[3]

    significance_list = []
    cross_section = cross_sections[i] * 1000 #convert pb to fb
    for i in range(len(low_region_efficiencies)):
        low_region_signal = luminosity * cross_section * low_region_efficiencies[i] #luminosity in fb^-1, cross_section in fb (converted)
        high_region_signal = luminosity * cross_section * high_region_efficiencies[i]

        z_low = z(luminosity, low_region_signal, low_region_background)
        z_high = z(luminosity, high_region_signal, high_region_background)

        if z_low > z_high:
            significance_list.append(z_low)
        else:
            significance_list.append(z_high)

    return significance_list
