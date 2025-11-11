# --- KMTNet Algorithm ---
# Authors: Atousa Kalantari, Somayeh Khakpash

"""
Implementation of the two parameter PSPL fitting procedure from Kim et al. (2018)
(https://iopscience.iop.org/article/10.3847/1538-3881/aaa47b).

This method models the microlensing single-lens light curve. The standard model involves three nonlinear parameters:
    - t0 (time of maximum magnification)
    - tE (Einstein ring crossing time)
    - u0 (minimum impact parameter)

To reduce computational cost, the algorithm operates in two regimes—high magnification and low magnification—
reducing the problem to just two nonlinear parameters: t0 and teff, where teff = tE * u0.

Model definitions:
    High Magnification regime:
        Ft_high(t, f_1, f_0, t0, t_eff) = |f_1| * Q^(-1/2) + |f_0|
        where Q = 1 + ((t - t0) / t_eff)^2

    Low Magnification regime:
        Ft_low(t, f_1, f_0, t0, t_eff) = |f_1| * [1 - (1 + Q/2)^(-2)]^(-1/2) + |f_0|
        where Q = 1 + ((t - t0) / t_eff)^2

Grid search setup:
    - t_eff ∈ [1, 100], stepped as:    t_eff_{k+1} = (1 + δ) * t_eff_k
    - For each t_eff, t0 grid:         t0_{k, l+1} = t0_{k, l} + δ * t_eff_k
    - Step size:                       δ = 1/3

For each (t0, t_eff) grid point, fitting is performed using data within the window: t0 ± Z * t_eff, with Z = 7.
Each window must contain at least 10 data points.

The best-fit parameters (t0, t_eff, f0, f1) are identified from this grid search.

Subsequently, a linear “flat” model is also fit to the data in the region t0_best ± 7 * t_eff_best to obtain chi2_flat.

We then compute the delta chi-squared metric as:
    delta_chi_squared_kmt = |chi_mlens - chi2_linearfit| / chi2_linearfit

Empirically, if delta_chi_squared_kmt > 0.9, we consider the event a microlensing candidate.


The chosen metric threshold (>0.9), window size (Z = 7) and at least 10 datapoints, step size (delta = 1/3), and teff range ([1, 100]) are based on Rubin cadence and observed days.
***These values must be checked and validated to ensure they are optimal.***
"""


# imports
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import warnings
# import weightedstats as ws

# --- Model Definitions ---
# Functions for high and low magnification regimes,
# params: (t, f_1, f_0, t0, t_eff)

def Ft_high(t, f_1, f_0, t0, t_eff):
    # High-mag analytic approximation 
    Q = 1 + ((t - t0) / t_eff)**2
    return np.abs(f_1) * (Q**(-1.0 / 2)) + np.abs(f_0)

def Ft_low(t, f_1, f_0, t0, t_eff):
     # Low-mag analytic approximation
    Q = 1 + ((t - t0) / t_eff)**2
    return np.abs(f_1) * (1 - (1 + Q / 2)**-2)**(-1.0 / 2) + np.abs(f_0)

 # --- Chi2 Functions for Minimization ---
def chi2_high(f_params, t, flux, flux_err, t0, teff):
    # Compute chi^2 for the high-mag model (minimize over f1, f0)
    f_1, f_0 = f_params
    model = Ft_high(t, f_1, f_0, t0, teff)
    inv_sigma2 = 1.0 / (flux_err**2)
    return np.sum((flux - model)**2 * inv_sigma2)

def chi2_low(f_params, t, flux, flux_err, t0, teff):
    # Compute chi^2 for the low-mag model
    f_1, f_0 = f_params
    model = Ft_low(t, f_1, f_0, t0, teff)
    inv_sigma2 = 1.0 / (flux_err**2)
    return np.sum((flux - model)**2 * inv_sigma2)

def run_kmtnet_fit(times, fluxes, flux_errors):

    """
    Fit a two-parameter PSPL (KMTNet-Algorithm) model to a light curve,
    using a grid search over t0 and t_eff,
    and return Delta_Chi2 along with best-fit parameters.
    """

    # Ensure inputs are numpy arrays
    times = np.asarray(times)
    fluxes = np.asarray(fluxes)
    flux_errors = np.asarray(flux_errors)

    # Filter out NaNs from input data before creating DataFrame
    valid_data_indices = ~np.isnan(times) & ~np.isnan(fluxes) & ~np.isnan(flux_errors) & (flux_errors > 0)
    times, fluxes, flux_errors = times[valid_data_indices], fluxes[valid_data_indices], flux_errors[valid_data_indices]

    if len(times) < 5: # KMTNet fit needs sufficient data points
        return None, None, None # Not enough data for KMTNet fit

    data_df = pd.DataFrame({
        'time': times,
        'flux': fluxes,
        'flux_err': flux_errors
    })



    # --- Grid Search: Build t0-teff grid for nonlinear fitting ---
    teff_min, teff_max = 1, 50
    teff_list, t0_tE_list = [], []
    current_teff = teff_min

    # Build teff grid (teff_{k+1} = (1 + delta) * teff_k)
    while current_teff <= teff_max:
        teff_list.append(current_teff)
        delta = 1/5 if current_teff < 1 else 1/3
        current_teff *= (1 + delta)

    # For each teff, build the grid of t0 values
    t0_min, t0_max = np.min(times), np.max(times)
    for teff in teff_list:
        t0_current = t0_min
        while t0_current <= t0_max:
            t0_tE_list.append([t0_current, teff])
            delta = 1/5 if teff < 1 else 1/3
            t0_current += delta * teff

    # If no grid was produced, exit
    if not t0_tE_list: return None, None, None


    # --- Main Grid Fit Loop ---
    param1, param2 = [], []  # Will store fit results for high and low mag regimes
    f_initial = [0.01, 0.99] # Initial guess for f_1, f_0

    for i, (t0_val, teff_val) in enumerate(t0_tE_list):
        # For each grid point, select data within the relevant window (t0 ± 7 teff)
        df_i = data_df[(data_df['time'] > (t0_val - 7 * teff_val)) & (data_df['time'] < (t0_val + 7 * teff_val))]

        if len(df_i) < 10:
            continue # Skip if not enough data in interval

        # Prepare arguments for minimize (t, flux, flux_err, t0_val, teff_val)
        args = (df_i['time'].values, df_i['flux'].values, df_i['flux_err'].values, t0_val, teff_val)

        try:
            # Fit the high-magnification model for current grid point
            result1 = minimize(chi2_high, f_initial, args=args, method='BFGS')
            # Compute chi2 for the entire dataset using best-fit parameters
            model_diff1 = data_df['flux'].values - Ft_high(data_df['time'].values, result1.x[0], result1.x[1], t0_val, teff_val)
            chi2_all1 = np.sum((model_diff1)**2 * (1.0 / (data_df['flux_err'].values**2)))


            # Fit the low-magnification model
            result2 = minimize(chi2_low, f_initial, args=args, method='BFGS')
            model_diff2 = data_df['flux'].values - Ft_low(data_df['time'].values, result2.x[0], result2.x[1], t0_val, teff_val)
            chi2_all2 = np.sum((model_diff2)**2 * (1.0 / (data_df['flux_err'].values**2)))

            # Store: [index, t0, teff, f1, f0, local_chi2, window_npts, global_chi2, npts]
            param1.append([i, t0_val, teff_val, result1.x[0], result1.x[1], result1.fun,
                            len(df_i), chi2_all1, len(data_df)])
            param2.append([i, t0_val, teff_val, result2.x[0], result2.x[1], result2.fun,
                            len(df_i), chi2_all2, len(data_df)])
        except Exception as e:
            # If optimization fails, print a warning but continue
            print(f"Warning: KMTNet minimize failed for iteration {i}: {e}")
            continue 

    # If no fits were successful, exit
    if not param1 and not param2: return None, None, None 

    # --- Select Best-fit Parameters ---
    # Find best fits: lowest chi2 on entire dataset for both regimes
    min_value1 = min(param1, key=lambda x: x[7])
    min_value2 = min(param2, key=lambda x: x[7])

    # Use the regime (high or low mag) with the best global chi2
    if min_value1 < min_value2:
              min_value = min_value1
              param = param1
              F_t = Ft_high
              which_regim = 'high'

    else:
              min_value = min_value2
              param = param2
              F_t = Ft_low
              which_regim = 'low'

    # Extract the precise parameters where chi2 is minimized
    for sublist in param:
        if sublist[7] == min_value[7]:
            parameter = sublist



    chi_mlens = parameter[7] # Minimum chi2 for microlensing fit
    t0 = parameter[1]
    t_eff = parameter[2]
    f1 = parameter[3]
    f0 = parameter[4]

    # --- Linear Fit for Flat Model ---
    # Fit a constant flux to the light curve in the same window
    data_df_interval = data_df[(data_df['time'] > (t0 - 7 * t_eff)) & (data_df['time'] < (t0 + 7 * t_eff))]

    if len(data_df_interval) == 0:
        # Fallback to global mean if window is empty
        mean_flux_interval = np.mean(data_df['flux'].values)
    else:
        # Ensure weights are valid (not zero or inf)
        weights = 1.0 / (data_df_interval['flux_err'].values**2)
        valid_weights_indices = ~np.isinf(weights) & ~np.isnan(weights) & (weights > 0)
        if np.sum(valid_weights_indices) > 0:
            # mean_flux_interval = ws.weighted_mean(data_df_interval['flux'].values[valid_weights_indices],
            #                                       weights[valid_weights_indices])
            mean_flux_interval = np.mean(data_df_interval['flux'].values[valid_weights_indices],
                                                  weights=weights[valid_weights_indices])
        else:
            mean_flux_interval = np.mean(data_df['flux'].values)

    # Compute chi2 for flat line fit
    chi2_linearfit = np.sum((data_df['flux'] - mean_flux_interval)**2 / (data_df['flux_err']) ** 2)


     # --- Compute Metric and Return ---
    if chi2_linearfit == 0:
        delta_chi_squared_kmt = 0
    else:
        delta_chi_squared_kmt = (abs(chi_mlens - chi2_linearfit) / chi2_linearfit)

    # Return: delta chi2, best-fit physical params
    # If delta_chi_squared_kmt > 0.9, the light curve would be a microlensing candidate.
    return delta_chi_squared_kmt, (t0, t_eff, f1, f0), which_regim
