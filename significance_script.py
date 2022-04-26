import os
from Zstats import Zdisc, Zexcl

#scan_path = "./restructured_careful_scan"
def significance_calculator(scan_path):
    luminosity = 139 #fb^-1
    background_low = 6.0 #nr of events
    background_low_uncertainty = 1.7
    background_high = 10.2
    background_high_uncertainty = 3.3

    discovery_list = []
    exclusion_list = []
    points = len([name for name in os.listdir(scan_path)])
    for i in range(points):
        current_folder = "ANALYSIS_" + str(i)
        data_file = scan_path + "/" + current_folder + "/Output/SAF/CLs_output_summary.dat"

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
        low_eff = 0
        high_eff = 0

        if first_line[0] == "-":
            low_exp = -1
        else:
            low_exp = float(first_line[0:9])

        if second_line[0] == "-":
            high_exp = -1
        else:
            high_exp = float(second_line[0:9])

        first_line = first_line[first_line.find('||'):]
        first_line = first_line[2:]

        second_line = second_line[second_line.find('||'):]
        second_line = second_line[2:]

        low_eff = float(first_line[0:9])
        high_eff = float(second_line[0:9])

        low_sig = low_exp * luminosity * low_eff
        high_sig = high_exp * luminosity * high_eff

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

        # print("FIRST LINE: ", first_line)
        # print("SECOND LINE: ", second_line)
        #
        # print("Low sig: ", low_sig)
        # print("High sig: ", high_sig)
        #
        # print("Z-DISC: ", Z_disc)
        # print("Z-EXCL: ", Z_excl)

    return [discovery_list, exclusion_list]
