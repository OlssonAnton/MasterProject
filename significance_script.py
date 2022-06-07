import numpy as np
from Zstats import Zdisc, Zexcl
from general_stuff import asimov_significance

def significance_calculator(scan_path, luminosity, cross_sections, low_region_efficiencies, high_region_efficiencies, calculator = 'Zstats', decreased_background_uncertainty = False):
    luminosity_scaling = luminosity / 139 #The background values are from a 139 fb^-1 search, assuming linear scaling of backgrounds and uncertainties
    low_uncertainty_ratio = 1.7 / 6.0
    high_uncertainty_ratio = 3.3 / 10.2
    low_region_background = 6.0 * luminosity_scaling #nr of events
    low_region_background_uncertainty = low_region_background * low_uncertainty_ratio
    high_region_background = 10.2 * luminosity_scaling
    high_region_background_uncertainty = high_region_background * high_uncertainty_ratio
    if decreased_background_uncertainty:
        high_region_background_uncertainty = high_region_background * 0.1
        low_region_background_uncertainty = low_region_background * 0.1

    discovery_array = np.zeros((len(cross_sections), len(cross_sections[0])), float)
    exclusion_array = np.zeros((len(cross_sections), len(cross_sections[0])), float)

    for i in range(len(exclusion_array)):
        for j in range(len(exclusion_array[i])):
            cross_section = 1000 * cross_sections[i][j] #Convert pb to fb
            low_region_signal = luminosity * cross_section * low_region_efficiencies[i][j] #luminosity in fb^-1, cross_section in fb (converted)
            high_region_signal = luminosity * cross_section * high_region_efficiencies[i][j]
            if calculator == 'asimov':
                z_low_region_sig = asimov_significance(low_region_signal, low_region_background, low_region_background_uncertainty).real
                z_high_region_sig = asimov_significance(high_region_signal, high_region_background, high_region_background_uncertainty).real
                if z_low_region_sig > z_high_region_sig:
                    discovery_array[i][j] = z_low_region_sig
                    exclusion_array[i][j] = z_low_region_sig
                else:
                    discovery_array[i][j] = z_high_region_sig
                    exclusion_array[i][j] = z_high_region_sig
            if calculator == 'Zstats':
                z_disc_low_region_signal = 0
                z_disc_high_region_signal = 0
                z_excl_low_region_signal = 0
                z_excl_high_region_signal = 0
                if (low_region_signal > 1000) or (high_region_signal > 1000): #Zstats can't handle too large number of events, compared to background
                    print('Too many signals for Zstats!')
                    z_disc_low_region_signal = asimov_significance(low_region_signal, low_region_background, low_region_background_uncertainty).real
                    z_disc_high_region_signal = asimov_significance(high_region_signal, high_region_background, high_region_background_uncertainty).real

                    z_excl_low_region_signal = asimov_significance(low_region_signal, low_region_background, low_region_background_uncertainty).real
                    z_excl_high_region_signal = asimov_significance(high_region_signal, high_region_background, high_region_background_uncertainty).real
                else:
                    z_disc_low_region_signal = Zdisc(low_region_signal, low_region_background, low_region_background_uncertainty)
                    z_disc_high_region_signal = Zdisc(high_region_signal, high_region_background, high_region_background_uncertainty)

                    z_excl_low_region_signal = Zexcl(low_region_signal, low_region_background, low_region_background_uncertainty)
                    z_excl_high_region_signal = Zexcl(high_region_signal, high_region_background, high_region_background_uncertainty)

                if z_disc_low_region_signal > z_disc_high_region_signal:
                    print(z_disc_low_region_signal)
                    discovery_array[i][j] = z_disc_low_region_signal
                else:
                    discovery_array[i][j] = z_disc_high_region_signal
                    print(z_disc_high_region_signal)

                if z_excl_low_region_signal > z_excl_high_region_signal:
                    print(z_excl_low_region_signal)
                    exclusion_array[i][j] = z_excl_low_region_signal
                else:
                    exclusion_array[i][j] = z_excl_high_region_signal
                    print(z_excl_high_region_signal)

    return [exclusion_array, discovery_array]
