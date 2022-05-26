import numpy as np
from Zstats import Zdisc, Zexcl
from general_stuff import asimov_significance

def significance_calculator(scan_path, luminosity, cross_sections, low_region_efficiencies, high_region_efficiencies, calculator):
    SM_irreducible_background = 0.4312 + 0.02746 * 2 #2 neutrino final state processes (SM background)
    low_region_background = 6.0 #nr of events
    low_region_background_uncertainty = 1.7
    high_region_background = 10.2
    high_region_background_uncertainty = 3.3

    asimov_array = np.zeros((len(cross_sections), len(cross_sections[0])), float)
    discovery_array = np.zeros((len(cross_sections), len(cross_sections[0])), float)
    exclusion_array = np.zeros((len(cross_sections), len(cross_sections[0])), float)

    for i in range(len(exclusion_array)):
        for j in range(len(exclusion_array[i])):
            cross_section = 1000 * cross_sections[i][j] #Convert pb to fb
            low_region_signal = luminosity * cross_section * low_region_efficiencies[i][j] #luminosity in fb^-1, cross_section in fb (converted)
            high_region_signal = luminosity * cross_section * high_region_efficiencies[i][j]

            if calculator == 'Zstats':
                z_disc_low_region_signal = 0
                z_disc_high_region_signal = 0
                z_excl_low_region_signal = 0
                z_excl_high_region_signal = 0
                if (low_region_signal > 500) or (high_region_signal > 500): #Zstats can't handle too large number of events, compared to background
                    print('Too many signals for Zstats!')
                    z_disc_low_region_signal = asimov_significance(low_region_signal, low_region_background, low_region_background_uncertainty)
                    z_disc_high_region_signal = asimov_significance(high_region_signal, high_region_background, high_region_background_uncertainty)

                    z_excl_low_region_signal = asimov_significance(low_region_signal, low_region_background, low_region_background_uncertainty)
                    z_excl_high_region_signal = asimov_significance(high_region_signal, high_region_background, high_region_background_uncertainty)
                else:
                    z_disc_low_region_signal = Zdisc(low_region_signal, low_region_background, low_region_background_uncertainty)
                    z_disc_high_region_signal = Zdisc(high_region_signal, high_region_background, high_region_background_uncertainty)

                    z_excl_low_region_signal = Zexcl(low_region_signal, low_region_background, low_region_background_uncertainty)
                    z_excl_high_region_signal = Zexcl(high_region_signal, high_region_background, high_region_background_uncertainty)

                if z_disc_low_region_signal > z_disc_high_region_signal:
                    discovery_array[i][j] = z_disc_low_region_signal
                else:
                    discovery_array[i][j] = z_disc_high_region_signal

                if z_excl_low_region_signal > z_excl_high_region_signal:
                    exclusion_array[i][j] = z_excl_low_region_signal
                else:
                    exclusion_array[i][j] = z_excl_high_region_signal

    return [asimov_array, discovery_array, exclusion_array]
