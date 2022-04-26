import os
from os.path import exists
import subprocess
import lhe_parser
import shutil
import gzip
from mpl_toolkits import mplot3d
from significance_script import significance_calculator
from general_stuff import read_and_write, generate_mvd_list, inclusive_range, number_of_points, edit_line
from myplots import z_plot, z_plot2, sigma_plot, z_contour_plot, z_contour_plot_old

def retrieve_cross_sections(scan_path):
    if(exists(scan_path + '/saved_cross_sections')):
        string_list = open(scan_path + '/saved_cross_sections').readlines()[0]
        fixed_string = ''.join(string_list.split('[', 1))
        fixed_string = ''.join(fixed_string.split(']', 1))
        fixed_string = fixed_string.replace(',', '')
        list_of_floats = list(fixed_string.split(' '))
        for i in range(len(list_of_floats)):
            list_of_floats[i] = float(list_of_floats[i])
        return list_of_floats

    number_of_scans = len(next(os.walk(scan_path))[1]) #points  #Number of simulations to be plotted
    all_cross_sections = []
    for i in range(number_of_scans):
        event_exist = True #Safeguard for opening file, sets to false if LHE file missing
        run_path = scan_path + '/run' + str(i + 1) #Path to folder containing all simulations
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

    with open(scan_path + '/saved_cross_sections', 'w') as save_file:
        save_file.write(str(all_cross_sections))

    return all_cross_sections

def run_simulation(mvd_properties, mtad_properties, run_name, production_directory):
    mtad_min = mtad_properties[0]
    mtad_end = mtad_properties[1]
    mtad_step = mtad_properties[2]
    mvd_min = mvd_properties[0]
    mvd_end = mvd_properties[1]
    mvd_step = mvd_properties[2]
    madevent_path = production_directory + '/bin/madevent'
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

def plot_stuff():
    mtad = inclusive_range(135, 235, 25)
    mvd = inclusive_range(80, 230, 25)
    print('NUMBER OF POINTS: ', number_of_points(mvd, mtad))

    #z_lists = significance_calculator('./gD1')
    #z_lists = significance_calculator('../../HEPTools/madanalysis5/madanalysis5/newmodel1point/')
    #print('Discovery Significance: ', z_lists[0])
    #print('Exclusion Significance: ', z_lists[1])

    alldiagrams_scan_path = '../../masterproject/tauDMproduction/alldiagrams/Events/careful_scan_pythia_delphes'
    resonance_scan_path = '../../masterproject/tauDMproduction/tadresonance/Events/careful_scan_pythia_delphes'
    newalldiagrams_scan_path = '../../masterproject/tauDMproduction/newalldiagrams/Events/newscan/'
    all_cross_sections = retrieve_cross_sections(newalldiagrams_scan_path)
    #print('CROSS SECTIONS: ', all_cross_sections)
    #print('CROSS SECTION LENGTH: ', len(all_cross_sections))

    SM_background = 0.4312 + 0.02746 * 2 #2 neutrino final state processes (SM background)
    L3 = 300 #Luminosity fb^-1, end of Run3
    #z_contour_plot_old(L3, all_cross_sections, SM_background, mvd, mtad)
    z_contour_plot(L3, all_cross_sections, SM_background, mvd, mtad)

if '__main__' == __name__:
    plot_stuff()
    #run_simulation([80, 230, 25], [135, 235, 25], [], '../../masterproject/tauDMproduction/newalldiagrams')
