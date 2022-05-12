import os
from os.path import exists
import subprocess
import shutil
import gzip
from general_stuff import inclusive_range, edit_line
from myplots import asimov_significance_color, asimov_significance_contour, exclusion_discovery_significance_contours

def run_simulation(mvd_properties, mtad_properties, dark_coupling, output_directory):
    mtad_min = mtad_properties[0]
    mtad_max = mtad_properties[1]
    mtad_step = mtad_properties[2]
    mvd_min = mvd_properties[0]
    mvd_max = mvd_properties[1]
    mvd_step = mvd_properties[2]
    mtad_list = inclusive_range(mtad_min, mtad_max, mtad_step)
    madevent_path = output_directory + '/bin/madevent'
    madanalysis_path = '../../HEPTools/madanalysis5/madanalysis5/bin/ma5'
    madevent_script_file = './madevent_script.txt'
    madanalysis_script_file = './madanalysis_script.txt'

    edit_line(madevent_script_file, 9, 'set gD ' + str(dark_coupling))

    event_directory = output_directory + '/Events'
    scan_name = 'mvd' + str(mvd_min) + '-' + str(mvd_max) + '_mtad' + str(mtad_min) + '-' + str(mtad_max)
    scan_directory = event_directory + '/mvd' + str(mvd_min) + '-' + str(mvd_max) + '_mtad' + str(mtad_min) + '-' + str(mtad_max) + '_gD' + str(dark_coupling)
    if not exists(scan_directory):
        os.mkdir(scan_directory) #Create directory where MadEvent and MadAnalysis outputs will be located
    print('Scan Directory:', scan_directory)

    save_file = open(scan_directory + '/saved_cross_sections', 'w')
    save_file_setup = ['[]\n', str(mvd_properties) + '\n', str(mtad_properties) + '\n', '[]\n', '[]\n']
    save_file_setup = ''.join(save_file_setup)
    save_file.write(save_file_setup)
    save_file.close()

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

    mtad = mtad_min
    mvd = mvd_min
    while(mtad <= mtad_max):
        while(mvd < mtad and mvd <= mvd_max):
            point_name = 'mvd' + str(mvd) + '_mtad' + str(mtad)
            point_directory = scan_directory + '/mvd' + str(mvd) + '_mtad' + str(mtad)
            if not exists(point_directory):
                os.mkdir(point_directory)
            print('Point Directory:', point_directory)
            edit_line(madevent_script_file, 7, 'set mvd ' + str(mvd))
            edit_line(madevent_script_file, 6, 'set mtad ' + str(mtad))
            edit_line(madevent_script_file, 0, 'launch ' + point_name)
            edit_line(madanalysis_script_file, 4, 'submit ' + point_name)
            ignored = 'False'
            #if mvd <= 300: #Doesnt seem to be stuff at more than mvd = 300 GeV, save time this way
            subprocess.run([madevent_path, madevent_script_file])
            if (exists(event_directory + '/' + point_name)):
                hepmc_unzipper(event_directory + '/' + point_name)
                current_hepmc_file = ''
                for file in os.listdir(event_directory + '/' + point_name):
                    if file.endswith('.hepmc'):
                        current_hepmc_file = file
                edit_line(madanalysis_script_file, 3, 'import ' + event_directory + '/' + point_name + '/' + current_hepmc_file)
                subprocess.run([madanalysis_path, '-R', '-s', madanalysis_script_file])

                for file in os.listdir(event_directory + '/' + point_name):
                    if file.endswith('events.lhe.gz'): #Find and unzip the event file
                        event_file_path = event_directory + '/' + point_name + '/unweighted_events.lhe'
                        zipped_event_file_path = event_file_path + '.gz'
                        with gzip.open(zipped_event_file_path, 'rb') as f_in:
                            with open(event_file_path, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out) #Unzip event file
                    if (not file.endswith('events.lhe.gz')) and (not file.endswith('events.lhe.gz')) and (not file.endswith('banner.txt')) and (not file == 'saved data') :
                        file_path = event_directory + '/' + point_name + '/' + file
                        print("Trying to remove: ", file_path)
                        if exists(file_path):
                            os.remove(file_path)
                for file in os.listdir(event_directory + '/' + point_name):
                    if exists(event_directory + '/' + point_name):
                        shutil.move(event_directory + '/' + point_name + '/' + file, point_directory + '/' + file)
                shutil.move('./' + point_name + '/Output/SAF/CLs_output_summary.dat', scan_directory + '/' + point_name + '/CLs_output_summary.dat')

            else: #madevent ran but a directory was not created, probably crashed (or ignored)
                os.mkdir(point_directory + '/madevent_ignored')
                os.mkdir(point_directory + '/madanalysis_ignored') #Add en empty madanalysis directory if MadGraph/MadAnalysis was not run
                ignored = 'True'
            ignored += '\n' #Prepare to write ignored status to file
            current_point_info = str([mvd, mtad]) + '\n'
            save_file = open(point_directory + '/saved_data', 'w')
            save_file_setup = ['[]\n', '[]\n', '[]\n', current_point_info, ignored] #cross section, low_eff, high_eff, mass_point, ignored?
            save_file_setup = ''.join(save_file_setup)
            save_file.write(save_file_setup)
            save_file.close()

            mvd += mvd_step
        edit_line(madevent_script_file, 7, 'set mvd ' + str(mvd_min))
        edit_line(madevent_script_file, 6, 'set mtad ' + str(mtad + mtad_step))
        mtad += mtad_step
        mvd = mvd_min

def plot(scan_paths, luminosity):
    #asimov_significance_color(scan_paths, luminosity)
    asimov_significance_contour(scan_paths, luminosity)
    #exclusion_discovery_significance_contours(scan_paths, luminosity)

if '__main__' == __name__:
    #Running from madgraph/various directory
    gD1_event_path = '../../masterproject/tauDMproduction/newalldiagrams/Events/gD1'
    gD3_event_path = '../../masterproject/tauDMproduction/newalldiagrams/Events/gD3'
    gD1_analysis_path = './gD1'
    gD3_analysis_path = './gD3'
    event_paths = [gD1_event_path, gD3_event_path]
    analysis_paths = [gD1_analysis_path, gD3_analysis_path]
    no_mixing_event_path = '../../masterproject/tauDMproduction/no_mixing/Events/gD_0.01'
    no_mixing_analysis_path = './gD_0.01'

    gD_001_path = '../../masterproject/tauDMproduction/no_mixing/Events/mvd5-455_mtad10-760_gD1'

    #SM_background = 0.4312 + 0.02746 * 2 #2 neutrino final state processes (SM background)
    luminosity = 139 #Luminosity = 300 fb^-1, end of Run3 (L3)

    #plot([gD_001_path], luminosity)
    run_simulation([5, 305, 50], [60, 560, 50], 0.01, '../../masterproject/tauDMproduction/no_mixing') #mvd, mtad, gD, output_directory
    #run_simulation([5, 5, 50], [100, 100, 50], [], 1, '../../masterproject/tauDMproduction/no_mixing') #mvd, mtad, gD, output_directory
