import os
from os.path import exists
import subprocess
import lhe_parser
import shutil
import gzip
from mpl_toolkits import mplot3d
from general_stuff import read_and_write, generate_mvd_mtad_lists, inclusive_range, number_of_points, edit_line
from myplots import asymptotic_limit_significance_color, asymptotic_limit_significance_contour, exclusion_discovery_significance_contours

def retrieve_cross_sections(event_path):
    if (exists(event_path + '/saved_data')):
        saved_cross_sections = open(event_path + '/saved_data').readlines()[0]
        if not saved_cross_sections == '[]\n':
            fixed_string = ''.join(saved_cross_sections.split('[', 1))
            fixed_string = ''.join(fixed_string.split(']', 1))
            fixed_string = fixed_string.replace(',', '')
            list_of_floats = list(fixed_string.split(' '))
            for i in range(len(list_of_floats)):
                list_of_floats[i] = float(list_of_floats[i])
            return list_of_floats
    else: #create and setup an empy data_file
        save_file = open(event_path + '/saved_data', 'w')
        save_file_setup = ['[]\n', '[]\n', '[]\n', '[]\n']
        save_file_setup = ''.join(save_file_setup)
        save_file.write(save_file_setup)
        save_file.close()

    number_of_scans = len(next(os.walk(event_path))[1]) #points  #Number of simulations to be plotted
    all_cross_sections = []
    for i in range(number_of_scans):
        event_exist = True #Safeguard for opening file, sets to false if LHE file missing
        run_path = event_path + '/run' + str(i + 1) #Path to folder containing all simulations
        if exists(run_path):
            event_path = run_path + '/unweighted_events.lhe'
            if not(exists(event_path)): #Is there an unzipped event file?
                zipped_event_path = event_path + '.gz'
                if exists(zipped_event_path): #Is there a zipped event file?
                    with gzip.open(zipped_event_path, 'rb') as f_in:
                        with open(event_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out) #Unzip event file
                else:
                    event_exist = False #There is no zipped or unzipped event file

            if event_exist:
                lhe = lhe_parser.EventFile(event_path)
                sum_of_weights = 0
                number_of_events = 0
                for event in lhe:
                    sum_of_weights += event.wgt
                    number_of_events += 1
                cross_section = sum_of_weights / number_of_events
                all_cross_sections.append(cross_section)
            else:
                all_cross_sections.append(0) #sigma = 0 indicates missing event file
                print('No zipped or unzipped event file in run' + str(i+1))
        else:
            print('No run' + str(i+1) + ' directory')
            all_cross_sections.append(0)
    save_file = open(event_path + '/saved_data')
    data = save_file.readlines()
    save_file.close()
    data[0] = str(all_cross_sections) + '\n'
    save_file = open(event_path + '/saved_data', 'w')
    new_data = ''.join(data)
    save_file.write(new_data)
    save_file.close()
    # with open(event_path + '/saved_data', 'w') as save_file:
    #     save_file.write(str(data))

    return all_cross_sections

def run_simulation(mvd_properties, mtad_properties, run_name, output_directory):
    mtad_min = mtad_properties[0]
    mtad_end = mtad_properties[1]
    mtad_step = mtad_properties[2]
    mvd_min = mvd_properties[0]
    mvd_end = mvd_properties[1]
    mvd_step = mvd_properties[2]
    mtad_list = inclusive_range(mtad_min, mtad_end, mtad_step)
    madevent_path = output_directory + '/bin/madevent'
    madanalysis_path = '../../HEPTools/madanalysis5/madanalysis5/bin/ma5'
    madevent_script_file = './madevent_script.txt'
    madanalysis_script_file = './madanalysis_script.txt'

    def setup_madevent_script():
        scan_file = open(madevent_script_file)
        string_list = scan_file.readlines()
        scan_file.close()
        if run_name:
            string_list[0] = 'launch ' + run_name[0] + '\n'
        else:
            string_list[0] = 'launch run1\n'
        string_list[6] = 'set mtad ' + str(mtad_min) + '\n'
        string_list[7] = 'set mvd ' + str(mvd_min) + '\n'
        scan_file = open(madevent_script_file, 'w')
        new_file_contents = ''.join(string_list)
        scan_file.write(new_file_contents)
        scan_file.close()
    def setup_madanalysis_script():
        opened_script_file = open(madanalysis_script_file)
        string_list = opened_script_file.readlines()
        opened_script_file.close()
        first_run_path = production_directory + '/Events/run1'
        first_hepmc_file = ''
        for file in os.listdir(first_run_path):
            if file.endswith('.hepmc'):
                first_hepmc_file = file
        string_list[3] = 'import ' + first_run_path + first_hepmc_file + '\n'
        opened_script_file = open(script_file, 'w')
        new_file_contents = ''.join(string_list)
        opened_script_file.write(new_file_contents)
        opened_script_file.close()
    setup_madevent_script()
    #setup_madanalysis_script()

    def hepmc_unzipper(run_path):
        current_hepmc_file = ''
        for file in os.listdir(run_path):
            if file.endswith('.hepmc.gz'):
                zip_path = run_path + '/' + file
                file_path = zip_path.replace('.gz', '')
                with gzip.open(zip_path, 'rb') as f_in:
                    with open(file_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(zip_path)

    run_number = 1
    mtad = mtad_min
    mvd = mvd_min
    while(mtad <= mtad_end):
        while(mvd < mtad and mvd <= mvd_end):
            subprocess.run([madevent_path, madevent_script_file])
            run_path = production_directory + '/Events/run' + str(run_number)
            hepmc_unzipper(run_path)
            current_hepmc_file = ''
            if (exists(run_path)):
                for file in os.listdir(run_path):
                    if file.endswith('.hepmc'):
                        current_hepmc_file = file
            edit_line(madanalysis_script_file, 3, 'import ' + run_path + '/' + current_hepmc_file)
            subprocess.run([madanalysis_path, '-R', '-s', madanalysis_script_file])

            read_and_write(madevent_script_file, 'launch run' + str(run_number), 'launch run' + str(run_number+1))
            run_number += 1
            read_and_write(madevent_script_file, 'set mvd ' + str(mvd), 'set mvd ' + str(mvd+mvd_step))
            mvd += mvd_step

            for file in os.listdir(run_path):
                if (not file.endswith('events.lhe.gz')) and (not file.endswith('banner.txt')) :
                    file = run_path + '/' + file
                    print("Trying to remove: ", file)
                    os.remove(file)

        read_and_write(madevent_script_file, 'set mtad ' + str(mtad), 'set mtad ' + str(mtad+mtad_step))
        read_and_write(madevent_script_file, 'set mvd ' + str(mvd), 'set mvd ' + str(mvd_min))
        mtad += mtad_step
        mvd = mvd_min

    mass_lists = generate_mvd_mtad_lists(mvd_min, mvd_step, mtad_list)
    mvd_data = mass_lists[0] + '\n'
    mtad_data = mass_lists[1] + '\n'
    save_file = open(event_path + '/saved_data', 'w')
    save_file_setup = ['[]\n', '[]\n', '[]\n', mvd_data, mtad_data]
    save_file_setup = ''.join(save_file_setup)
    save_file.write(save_file_setup)
    save_file.close()


def plot(analysis_paths, event_paths):
    def retrieve_mass_points():
        saved_data = []
        for path in event_paths:
            save_file = open(path + '/saved_data', 'r')
            saved_data.append(save_file.readlines())
        mass_lists = [saved_data[0][3], saved_data[0][4]]
        for i in range(len(mass_lists)):
            mass_lists[i] = ''.join(mass_lists[i].split('[', 1))
            mass_lists[i] = ''.join(mass_lists[i].split(']', 1))
            mass_lists[i] = mass_lists[i].replace(',', '')
            list_of_floats = list(mass_lists[i].split(' '))
            for j in range(len(list_of_floats)):
                list_of_floats[j] = float(list_of_floats[j])
            mass_lists[i] = list_of_floats
        mvd_long_list = mass_lists[0]
        mtad_long_list = mass_lists[1]
        mvd_step = 0
        if len(mvd_long_list) > 1:
            if not mvd_long_list[mvd_long_list.index(mvd_long_list[-1])] == mvd_long_list[mvd_long_list.index(mvd_long_list[-1]) - 1]:
                mvd_step = mvd_long_list[mvd_long_list.index(mvd_long_list[-1])] - mvd_long_list[mvd_long_list.index(mvd_long_list[-1]) - 1]
        mtad_step = 0
        if len(mtad_long_list) > 1:
            if not mtad_long_list[mtad_long_list.index(mtad_long_list[-1])] == mtad_long_list[mtad_long_list.index(mtad_long_list[-1]) - 1]:
                mtad_step = mtad_long_list[mtad_long_list.index(mtad_long_list[-1])] - mtad_long_list[mtad_long_list.index(mtad_long_list[-1]) - 1]
        mvd = inclusive_range(mvd_long_list[0], mvd_long_list[-1], mvd_step)
        mtad = inclusive_range(mtad_long_list[0], mtad_long_list[-1], mtad_step)
        return [mvd, mtad]
    mass_points = retrieve_mass_points()
    mvd = mass_points[0]
    mtad = mass_points[1]

    print('NUMBER OF POINTS: ', number_of_points(mvd, mtad))

    cross_section_set = [retrieve_cross_sections(gD1_event_path), retrieve_cross_sections(gD3_event_path)]
    # print("gD1 cross sections", retrieve_cross_sections(gD1_event_path))
    # print("gD3 cross sections", retrieve_cross_sections(gD3_event_path))

    #SM_background = 0.4312 + 0.02746 * 2 #2 neutrino final state processes (SM background)
    luminosity = 139 #Luminosity = 300 fb^-1, end of Run3 (L3)
    #asymptotic_limit_significance_color(L3, cross_section_set, mvd, mtad, ['./gD3'])
    #asymptotic_limit_significance_contour(luminosity, cross_section_set, mvd, mtad, ['./gD1', './gD3'])
    exclusion_discovery_significance_contours(luminosity, cross_section_set, mvd, mtad, analysis_paths, event_paths)

if '__main__' == __name__:
    #Running from madgraph/various directory
    gD1_event_path = '../../masterproject/tauDMproduction/newalldiagrams/Events/gD1'
    gD3_event_path = '../../masterproject/tauDMproduction/newalldiagrams/Events/gD3'
    gD1_analysis_path = './gD1'
    gD3_analysis_path = './gD3'
    event_paths = [gD1_event_path, gD3_event_path]
    analysis_paths = [gD1_analysis_path, gD3_analysis_path]
    plot(analysis_paths, event_paths)
    #run_simulation([5, 230, 25], [35, 235, 25], [], '../../masterproject/tauDMproduction/no_mixing')
