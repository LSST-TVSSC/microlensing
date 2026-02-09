# Filter for microlensing events designed for the ANTARES broker

# Authors: Somayeh Khakpash, Natasha Abrams, Rachel Street, Atousa Kalantari

import antares
import antares.devkit as dk
from statsmodels.stats.weightstats import DescrStatsW
import numpy as np
from astropy.table import MaskedColumn
import warnings
import astropy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import skew
from KMTNET_Algorithm import run_kmtnet_fit

# Initialize development kit client
dk.init()

# Define a Paczyński microlensing model
def paczynski(t, t0, u0, tE, F_s):
    """
    Paczyński microlensing light curve model
    t0 : peak time
    u0 : impact parameter
    tE : Einstein crossing time
    F_s : source flux
    F_b : blended flux
    """
    u = np.sqrt(u0**2 + ((t - t0) / tE) ** 2)
    A = (u**2 + 2) / (u * np.sqrt(u**2 + 4))
    return F_s * (A) + (1-F_s)
    
def fit_paczynski(times, mags, flxs, flx_errs):
    """
    Fit the Paczyński microlensing model to flux data.
    Returns best-fit parameters and chi-squared value.
    """
    if len(times) < 4:
        return None, None  # Not enough data

    # initial guesses
    t0_guess = times[np.argmin(flxs)]
    u0_guess = 1.0 / (np.max(flxs))

    tE_guess = 20.0
    F0_guess = 0.5


    initial_guess = [t0_guess, u0_guess, tE_guess, F0_guess]

    bounds = (
                        [times.min() - 50, 0, 1.0, 0.0],
                        [times.max() + 50, np.inf, 500.0, 1]
                    )


    try:
        popt, _ = curve_fit(
            paczynski,
            times, flxs,
            p0=initial_guess,
            sigma=flx_errs,
            bounds=bounds,

            maxfev=5000
        )
        chi2 = np.sum(((flxs - paczynski(times, *popt)) / flx_errs) ** 2) / len(times)
        return popt, chi2
    except Exception as e:
        print(f"  Paczynski fitting error: {e}")
        return None, None


def mag_to_flux(mag, F0=1.0):
    """
    Convert magnitude to flux.

    Parameters:
    - mag : magnitude (float or array)
    - F0 : reference flux (zeropoint), default=1.0 for relative flux

    Returns:
    - flux : flux corresponding to the magnitude
    """
    flux = F0 * 10 ** (-0.4 * mag)
    flux = flux / np.min(flux)
    return flux


def magerr_to_fluxerr(mag, mag_err, F0=1.0):
    """
    Convert magnitude uncertainty to flux uncertainty.

    Parameters:
    - mag : magnitude value or array
    - mag_err : magnitude uncertainty value or array
    - F0 : zeropoint flux (default=1.0 for relative flux)

    Returns:
    - flux_err : flux uncertainty
    """
    flux = mag_to_flux(mag, F0)
    flux_err = 0.4 * np.log(10) * flux * mag_err
    return flux_err


class microlensing(dk.Filter):    
    INPUT_LOCUS_PROPERTIES = [
        'ztf_object_id',
    ]

    REQUIRED_TAGS = ['lc_feature_extractor']

    OUTPUT_TAGS = [
        {
            'name': 'microlensing_candidate',
            'description': 'Locus - a transient candidate - exhibits a microlensing-like variability',
        }
    ]


    def make_lc(self, locus):

        with warnings.catch_warnings():
            # The cast of locus.timeseries: astropy.table.Table to a pandas
            # dataframe results in the conversion of some integer-valued
            # columns to floating point represntation. This can result in a
            # number of noisy warning so we will catch & ignore them for the
            # next couple of lines.
            warnings.simplefilter("ignore", astropy.table.TableReplaceWarning)
            df = locus.timeseries.to_pandas()

        data = df[['ant_mjd', 'ztf_fid', 'ztf_magpsf', 'ztf_sigmapsf']]
        
        dn = data.dropna()
        times=dn['ant_mjd'][dn['ztf_fid']==1]
        mags = dn['ztf_magpsf'][dn['ztf_fid']==1]
        mags_err = dn['ztf_sigmapsf'][dn['ztf_fid']==1]
        flxs = mag_to_flux(mags)
        flx_errs = magerr_to_fluxerr(mags, mags_err)

    def is_known_other_phenomenon(self, locus, locus_params):
        """
        Method to check the locus' pre-existing parameters indicated that it has
        been identified or is likely to be a variable of a type other than microlensing

        :param locus:
        :param locus_params:
        :return: boolean
        """

        # Default result is not a known variable
        known_var = False
        
        # Tunable detection thresholds.
        # Ref: Sokolovsky et al. 2016: https://ui.adsabs.harvard.edu/abs/2017MNRAS.464..274S/abstract
        period_peak_sn_threshold = 20.0  # Based on tests with ZTF alerts
        stetson_k_threshold = 0.8  # The expected K-value for a constant lightcurve with Gaussian noise

        # Check for periodicity
        if locus_params['properties']['feature_period_s_to_n_0_magn_r'] >= period_peak_sn_threshold:
            known_var = True

        # Check Stetson-K index
        if locus_params['properties']['feature_stetson_k_magn_r'] <= stetson_k_threshold:
            known_var = True

        # Check whether this event is associated with a GW event
        if 'plausible_gw_events_assoc' in locus_params['properties'].keys():
            known_var = True

        # If the alert has parameters from JPL Horizons, then it is likely cause by
        # a Solar System object
        if 'horizons_targetname' in locus_params['properties'].keys():
            known_var = True

        # Check whether the ANTARES crossmatch against known galaxy catalogs threw up any matches
        # The locus.catalog_objects attribute is a dictionary of lists of known objects for each
        # catalogs.  If a match has been found, then the key for the corresponding catalog will be
        # in the list of keys.  So we can use that to check for matches with galaxy catalogs.
        # Of those available in the list the Gemini NIR survey of known quasars is the closest
        if 'gnirs_dqs' in locus.catalog_objects.keys():
            known_var = True
                
        return known_var

    def calculate_eta(mag):
        """ Via puzle https://github.com/jluastro/puzle/blob/main/puzle/stats.py"""
        delta = np.sum((np.diff(mag)*np.diff(mag)) / (len(mag)-1))
        variance = np.var(mag)
        eta = delta / variance
        return eta

    def return_eta_residual_slope_offset():
        """ 
        Via puzle https://github.com/jluastro/puzle/blob/main/puzle/cands.py
        TODO is 6 months and a year - calculate slope and intercept based on real Rubin data
        """
        slope = 3.8187919463087248
        offset = -0.07718120805369133
        return slope, offset

    def is_microlensing_candidate(self, locus, times, mags, errors):
        """
        Example of a set of Microlensing detection criteria
        """
        if len(times) < 10:  # Too few data points
            return False

        # Extract the full parameter set from the locus and the alert
        locus_params = locus.to_dict()

        # Use the pre-calculated properties of the locus to eliminate those
        # which show signs of variability, e.g. in their periodicity signature or
        # the Stetson-K index
        known_var = self.is_known_other_phenomenon(locus, locus_params)
        if known_var:
            return False

        # Sort data by time
        sorted_idx = np.argsort(times)
        times, mags, errors = times[sorted_idx], mags[sorted_idx], errors[sorted_idx]

        # 1. Check for smoothness (low skewness means symmetric light curve)
        # TODO: Check for threshold with parallax and maybe remove or lower threshold
        if abs(skew(mags)) > 1:
            return False

        # TODO is 6 months and a year - calculate this based on percentile of real data
        eta_thresh = 1.255 # Avg from ZTF level 2 (low eta)
        # Do check for existance since if there's only one band of data, only one will exist
        eta_r_exists = 'feature_eta_e_magn_r' in locus_params['properties'].keys()
        eta_g_exists = 'feature_eta_e_magn_g' in locus_params['properties'].keys()
        eta_r = locus_params['properties']['feature_eta_e_magn_r']
        eta_g = locus_params['properties']['feature_eta_e_magn_g']
        if eta_r_exists and eta_g_exists:
            if eta_r >= eta_thresh and eta_g >= eta_thresh:
                return False
        elif eta_r_exists:
            if eta_r >= eta_thresh:
                return False
        elif eta_g_exists:
            if eta_g >= eta_thresh:
                return False

        # 2. Check variability (microlensing should have a clear peak)
        # Decrease threshold with longer baseline
        # Q for broker - 365 days or full lightcurve?
        if np.ptp(mags) < 0.5:  # Peak-to-peak magnitude difference
            return False

        flxs = mag_to_flux(mags)
        flx_errs = magerr_to_fluxerr(mags, errors)

        # 3. Perform a lightweight template fit (Paczyński model)

        popt, chi2_paczynski = fit_paczynski(times, mags, flxs, flx_errs)
        resid = flxs - paczynski(times, *popt)
        chi2 = np.sum((resid / flx_errs) ** 2) / len(times)

        # 4. Apply a simple chi2 threshold
        if chi2 >= 2:  # Poor-fit light curves fails
            return False
        # except RuntimeError:
        #     return False  # Fit failed

        eta_resid = self.calculate_eta(resid)
        eta_slope, eta_offset = self.return_eta_residual_slope_offset()
        if eta_r_exists and eta_g_exists:
            if (eta_resid < eta_r*eta_slope + eta_offset) and (eta_resid < eta_g*eta_slope + eta_offset):
                return False
        elif eta_r_exists:
            if (eta_resid < eta_r*eta_slope + eta_offset):
                return False
        elif eta_g_exists:
            if (eta_resid < eta_g*eta_slope + eta_offset):
                return False

        # TODO - Rache potentially query full lightcurve if not already there and if possible

        # TODO - Natasha add parallax microlensing fit

        return True
    def run(self, locus):
        print('Processing Locus:', locus.locus_id)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", astropy.table.TableReplaceWarning)
            df = locus.timeseries.to_pandas()

        data = df[['ant_mjd', 'ztf_fid', 'ztf_magpsf', 'ztf_sigmapsf']].dropna()

        
        
        # Split into g-band and i-band
        for band in [1, 2]:  # 1 = g-band, 2 = i-band
            band_data = data[data['ztf_fid'] == band]
            times, mags, errors = band_data['ant_mjd'].values, band_data['ztf_magpsf'].values, band_data['ztf_sigmapsf'].values
            
            if self.is_microlensing_candidate(locus, times, mags, errors):
                print(f'Locus {locus.locus_id} is a microlensing candidate in band {band}')
                locus.tag('microlensing_candidate')
        
        
        return