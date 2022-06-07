import os
from os.path import exists
import subprocess
import shutil
import gzip
from general_stuff import inclusive_range, edit_line

def run_simulation(mvd_properties, mtad_properties, dark_coupling, mtap, output_directory):
    mtad_min = mtad_properties[0] #Retrieve mass data for odd dark tau and DM candidate
    mtad_max = mtad_properties[1]
    mtad_step = mtad_properties[2]
    mvd_min = mvd_properties[0]
    mvd_max = mvd_properties[1]
    mvd_step = mvd_properties[2]
    madevent_path = output_directory + '/bin/madevent' #Define file paths to scripts and executables
    madanalysis_path = '../../HEPTools/madanalysis5/madanalysis5/bin/ma5'
    madevent_script_file = './madevent_script.txt'
    madanalysis_script_file = './madanalysis_script.txt'
    event_directory = output_directory + '/Events'

    edit_line(madevent_script_file, 9, 'set gD ' + str(dark_coupling)) #Set the value for the dark coupling from the input parameter
    edit_line(madevent_script_file, 10, 'set mtap ' + str(mtap))
    scan_name = 'mvd' + str(mvd_min) + '-' + str(mvd_max) + '_mtad' + str(mtad_min) + '-' + str(mtad_max) + '_gD' + str(dark_coupling)
    scan_directory = event_directory + '/' + scan_name
    if not exists(scan_directory):
        os.mkdir(scan_directory) #Create directory where MadEvent and MadAnalysis outputs will be located

    #Setup the save file which will contain cross sections, the scan region, low- and high signal region detector effiencies, relic_densities
    save_file = open(scan_directory + '/saved_scan_data', 'w')
    save_file_setup = ['[]\n', str(mvd_properties) + '\n', str(mtad_properties) + '\n', '[]\n', '[]\n', '[]\n', '[]\n', str(dark_coupling) + '\n', str(mtap) + '\n']
    save_file_setup = ''.join(save_file_setup)
    save_file.write(save_file_setup)
    save_file.close()

    def hepmc_unzipper(run_path): #A help function which unzips and removes the zipped version of all HepMC files at run_path
        current_hepmc_file = ''
        for file in os.listdir(run_path):
            if file.endswith('.hepmc.gz'):
                zip_path = run_path + '/' + file
                file_path = zip_path.replace('.gz', '')
                with gzip.open(zip_path, 'rb') as f_in:
                    with open(file_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(zip_path) #Don't keep the zipped version

    mtad = mtad_min
    mvd = mvd_min
    while(mtad <= mtad_max): # Loop for all masses of dark tau odd
        while(mvd < mtad and mvd <= mvd_max): #Loop for all masses of DM candidate below the kinematical limit
            point_name = 'mvd' + str(mvd) + '_mtad' + str(mtad) # Name of the mass point
            point_directory = scan_directory + '/mvd' + str(mvd) + '_mtad' + str(mtad) # Path to the directory of the mass point
            if not exists(point_directory):
                os.mkdir(point_directory) #Creaty empty directory for the mass point
            edit_line(madevent_script_file, 7, 'set mvd ' + str(mvd)) #Update the madevent and madanalysis scripts with the current points and output directories
            edit_line(madevent_script_file, 6, 'set mtad ' + str(mtad))
            edit_line(madevent_script_file, 0, 'launch ' + point_name)
            edit_line(madanalysis_script_file, 4, 'submit ' + point_name)
            subprocess.run([madevent_path, madevent_script_file])
            crashed = 'False' #Flag that will be set to 'True' should MadEvent fail to produce the output files/directory
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

            else: #madevent ran but a directory was not created, likely due to internal MadEvent crash
                crashed = 'True' #Flag will be saved to the mass point save file to indicate that MadEvent crashed at this point only
            crashed += '\n' #Prepare to write ignored status to file
            current_point_info = str([mvd, mtad]) + '\n'
            save_file = open(point_directory + '/saved_data', 'w')
            save_file_setup = [current_point_info, crashed] #mass_point, ignored?
            save_file_setup = ''.join(save_file_setup)
            save_file.write(save_file_setup)
            save_file.close()

            mvd += mvd_step #Prepare next point and loop for same mtad again if mvd + mvd_step<= mvd_max and mvd + mvd_step < mtad
        edit_line(madevent_script_file, 7, 'set mvd ' + str(mvd_min)) #Reset mvd if mvd_max or kinematical limit for current mtad was reached
        edit_line(madevent_script_file, 6, 'set mtad ' + str(mtad + mtad_step)) #Increase mtad
        mtad += mtad_step #Update loop parameters to match updated MadEvent script data
        mvd = mvd_min
