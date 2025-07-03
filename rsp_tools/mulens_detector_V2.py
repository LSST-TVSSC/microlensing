import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.lines as mlines
import itertools
from astropy.coordinates import SkyCoord
import astropy.units as u
import argparse
from antares_filter.KMTNET_Algorithm import run_kmtnet_fit, Ft_high, Ft_low
from lsst.rsp import get_tap_service
from lsst.daf.butler import Butler
import lsst.afw.display as afw_display
import lsst.sphgeom as sphgeom
import lsst.geom as geom
from tqdm import tqdm
from lsst.utils.plotting import (get_multiband_plot_colors,
                                 get_multiband_plot_symbols,
                                 get_multiband_plot_linestyles)

def chunked(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def process_light_curve(dict_lc):
            ID = list(dict_lc.keys())[0]
            rband_ind = np.where(dict_lc[ID]['band']=='r')
            time = dict_lc[ID]['time'][rband_ind]
            flux = dict_lc[ID]['flux'][rband_ind]
            flux_err = dict_lc[ID]['flux_err'][rband_ind]
            delta_chi_squared_kmt, (t0, t_eff, f1, f0), which_regim = run_kmtnet_fit(time, flux, flux_err)
            return delta_chi_squared_kmt, (t0, t_eff, f1, f0), which_regim



def run(args):

    # Set up the data butler
    service = get_tap_service("tap")
    assert service is not None
    butler = Butler('dp1', collections="LSSTComCam/DP1")
    assert butler is not None


    # Query for data in region
    ra_cen = args.ra
    dec_cen = args.dec
    radius = args.radius
    region = sphgeom.Region.from_ivoa_pos(f"CIRCLE {ra_cen} {dec_cen} {radius}")

    colors = ['purple', 'r', 'g', 'orange', 'blue', 'k']
    all_bands = ['r', 'i', 'g', 'z', 'u', 'y']

    # Fetch a table of all Objects within this search radius
    query = "SELECT ra, dec, diaObjectId, "\
        "r_psfFluxMean, nDiaSources, r_scienceFluxMean, "\
        "r_psfFluxNdata, r_psfFluxMax,r_psfFluxMin, r_psfFluxLinearSlope "\
        "FROM dp1.DiaObject "\
        "WHERE CONTAINS (POINT('ICRS', ra, dec), "\
        "CIRCLE('ICRS'," + str(ra_cen) + ", "\
        + str(dec_cen) + ", 2)) = 1 "

    job = service.submit_job(query)
    job.run()
    job.wait(phases=['COMPLETED', 'ERROR'])
    print('Job phase is', job.phase)
    if job.phase == 'ERROR':
        job.raise_if_error()
    assert job.phase == 'COMPLETED'
    objtab = job.fetch_result().to_table()
    objtab = objtab[objtab['r_scienceFluxMean'] < 10**5.5]

    n_chuncks = int(len(np.array(objtab['diaObjectId'].data))/500)+1
    
    with open('./mulens_candidates.dat', 'w') as fout:
        fout.write('# diaObjectId   coord_ra  coord_dec  delta_chi_squared_kmt  t0  t_eff  f1  f0 which_regim\n')
    
        chunk_counter = 0
        IDs = np.array(objtab['diaObjectId'].data)
        for chunk in chunked(IDs, 500):
            print('Starting chunck %i/%i:'%(chunk_counter, n_chuncks))
            ids_str = ",".join(str(i) for i in chunk)  
            
            
            query = f"""
                SELECT fsodo.diaObjectId, fsodo.coord_ra, fsodo.coord_dec,
                       fsodo.visit, fsodo.band,
                       fsodo.psfDiffFlux, fsodo.psfDiffFluxErr,
                       fsodo.psfFlux AS psfFlux, fsodo.psfFluxErr,
                       vis.expMidptMJD
                FROM dp1.ForcedSourceOnDiaObject AS fsodo
                JOIN dp1.Visit AS vis ON vis.visit = fsodo.visit
                WHERE fsodo.diaObjectId IN ({ids_str})
            """
    
            job = service.submit_job(query)
            job.run()
            job.wait(phases=['COMPLETED', 'ERROR'])
    
            
            if job.phase == 'ERROR':
                job.raise_if_error()
            assert job.phase == 'COMPLETED'
    
            forced_sources = job.fetch_result().to_table()
            all_data = forced_sources #[forced_sources['band']=='r']
            all_data_list = []
            
            for ID in ids_str.split(','):
                target = objtab[objtab['diaObjectId']== int(ID)]
                if len(all_data[all_data['diaObjectId'] == int(ID)][all_data[all_data['diaObjectId'] == int(ID)]['band']=='r'])<10:
                    continue
                all_data_list.append(
                    {ID: {'ra_coor':target['ra'].data[0],
                          'dec_coor':target['dec'].data[0],
                          'band': np.array(all_data[all_data['diaObjectId'] == int(ID)]['band'].data),
                         'time': np.array(all_data[all_data['diaObjectId'] == int(ID)]['expMidptMJD'].data),
                         'flux': np.array(all_data[all_data['diaObjectId'] == int(ID)]['psfFlux'].data),
                         'flux_err': np.array(all_data[all_data['diaObjectId'] == int(ID)]['psfFluxErr'].data),
                          'fluxes_diff': np.array(all_data[all_data['diaObjectId'] == int(ID)]['psfDiffFlux'].data),
                         'flux_err_diff': np.array(all_data[all_data['diaObjectId'] == int(ID)]['psfDiffFlux'].data)
                          }}
                )
            print('Found %i out of 500 objects with at least 5 flux measurements in r band.'%(len(all_data_list)))
    
            
            counter = 0
            for data in tqdm(all_data_list):
                delta_chi_squared_kmt, (t0, t_eff, f1, f0), which_regim = process_light_curve(data)
                if delta_chi_squared_kmt > 0.9:
                    fout.write(f"{ID}  {dict_lc[ID]['ra_coor']}  {dict_lc[ID]['dec_coor']}  {delta_chi_squared_kmt}  {t0}  {t_eff}  {f1}  {f0} {which_regim}\n")


                    fig, axs = plt.subplots(1, 2, figsize=(15,5))
                    ID = list(data.keys())[0]
                    for b, band in enumerate(all_bands):
                        band_ind = np.where(data['band'] == band)
                        rband = np.where(data['band'] == 'r')
                        times = data['time'][band_ind]
                        fluxes = data['flux'][band_ind]
                        flux_err = data['flux_err'][band_ind]
                    
                        axs[1].errorbar(times, 
                                        fluxes, 
                                        yerr=flux_err, 
                                        fmt='o',
                                        color=colors[b], 
                                        label=all_bands[b]+' band')
                        #plt.errorbar(times, fluxes_diff,yerr= flux_err_diff, fmt='o', color=colors[b], label=all_bands[b]+' band')
                    
                    axs[0].errorbar(data['time'][rband], 
                                   data['flux'][rband], 
                                    yerr=data['flux_err'][rband],
                                    fmt='o',
                                   color='r', 
                                   label='r band')
                    times_new = np.linspace(min(times_new),
                                           max(times_new),
                                           100)
                    if which_regim=='high':
                        model = Ft_high(times_new, f1, f0, t0, t_eff)
                    if which_regim=='low':
                        model = Ft_low(times_new, f1, f0, t0, t_eff)
                    
                    axs[0].plot(times_new, 
                                model,
                                color='g',
                               label = 'Best-fit model')
                    plt.suptitle('ID=%i, metric=%.3f'%(target["diaObjectId"].data[0], delta_chi_squared_kmt))
                    axs[0].set_ylabel('Flux')
                    axs[0].set_xlabel('Time (MJD)')
                    axs[1].set_xlabel('Time (MJD)')
                    
                    axs[0].legend()
                    axs[1].legend()
                    plt.savefig('./../candidates/%i.png'%target["diaObjectId"].data[0])
                    plt.close()
                    

            print('Finished Chunck %i/%i.'%(chunk_counter, n_chuncks))
            chunk_counter +=1
            


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('ra', help='Field center RA in decimal degrees')
    parser.add_argument('dec', help='Field center Dec in decimal degrees')
    parser.add_argument('radius', help='Search radius in decimal degrees')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = get_args()
    run(args)
