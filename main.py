from simulation import run_simulation
from myplots import exclusion_discovery_significance_contours

def main():
    luminosity = [139] #fb^-1

    gD_001_path = '../../masterproject/tauDMproduction/no_mixing/Events/mvd5-305_mtad10-560_gD0.01'
    gD_1_path = '../../masterproject/tauDMproduction/no_mixing/Events/mvd5-455_mtad10-760_gD1'

    exclusion_discovery_significance_contours([gD_001_path], luminosity, decreased_background_uncertainty = False,
    relic_density_constraint = False, width_mass_ratio_lines = False)

    #run_simulation([405, 505, 50], [760, 760, 50], 0.01, 800, '../../masterproject/tauDMproduction/no_mixing') #mvd, mtad, gD, output_directory

if '__main__' == __name__:
    main()
