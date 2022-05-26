from myplots import asimov_significance_color, asimov_significance_contour, exclusion_discovery_significance_contours
from simulation import run_simulation

def main():
    luminosity = 300 #fb^-1
    #plot() decides what type of data is to be plotted
    def plot(scan_paths, luminosity):
        #asimov_significance_color(scan_paths, luminosity)
        #asimov_significance_contour(scan_paths, luminosity)
        exclusion_discovery_significance_contours(scan_paths, luminosity)

    gD_001_path = '../../masterproject/tauDMproduction/no_mixing/Events/mvd5-305_mtad10-560_gD0.01'
    gD_1_path = '../../masterproject/tauDMproduction/no_mixing/Events/mvd5-455_mtad10-760_gD1'

    #plot([gD_001_path], luminosity)
    plot([gD_1_path], luminosity)
    plot([gD_001_path, gD_1_path], luminosity)

    #run_simulation([5, 5, 50], [10, 10, 50], 0.01, '../../masterproject/tauDMproduction/no_mixing') #mvd, mtad, gD, output_directory

if '__main__' == __name__:
    main()
