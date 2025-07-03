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
    query = f"""
                SELECT
                  fsodo.diaObjectId,
                  fsodo.coord_ra,
                  fsodo.coord_dec,
                  COUNT(*) AS positive_flux_count
                FROM dp1.ForcedSourceOnDiaObject AS fsodo
                JOIN dp1.Visit AS vis ON fsodo.visit = vis.visit
                WHERE
                  fsodo.band = 'r'
                  AND fsodo.psfFlux > 0
                  AND CONTAINS(
                    POINT('ICRS', fsodo.coord_ra, fsodo.coord_dec),
                    CIRCLE('ICRS', {ra_cen}, {dec_cen}, {radius})
                  ) = 1
                GROUP BY
                  fsodo.diaObjectId,
                  fsodo.coord_ra,
                  fsodo.coord_dec
                HAVING
                  COUNT(*) > 10
            """
    
    job = service.submit_job(query)
    job.run()
    job.wait(phases=['COMPLETED', 'ERROR'])
    print('Job phase is', job.phase)
    if job.phase == 'ERROR':
        job.raise_if_error()
    assert job.phase == 'COMPLETED'
    objtab = job.fetch_result().to_table()
    
    
    
    # multiprocessing.set_start_method('spawn')
    num_cores = 4
    print('Number of cores are:', num_cores)
    n_chuncks = int(len(np.array(objtab['diaObjectId'].data))/500)+1
    
    with open('./mulens_candidates.dat', 'w') as fout:
        fout.write('# diaObjectId   coord_ra  coord_dec  delta_chi_squared_kmt  t0  t_eff  f1  f0\n')
    
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
            all_data = forced_sources[forced_sources['band']=='r']
            all_data_list = []
            for ID in ids_str.split(','):
                target = objtab[objtab['diaObjectId']== int(ID)]
                all_data_list.append(
                    {ID: {'ra_coor':target['coord_ra'],
                          'dec_coor':target['coord_dec'],
                         'time': np.array(all_data[all_data['diaObjectId'] == int(ID)]['expMidptMJD'].data),
                         'flux': np.array(all_data[all_data['diaObjectId'] == int(ID)]['psfFlux'].data),
                         'flux_err': np.array(all_data[all_data['diaObjectId'] == int(ID)]['psfFluxErr'].data)
                          }}
                )
    
            
            counter = 0
            while len(all_data_list) > 0:
                print('Step: ', counter)
                num_to_process = min(num_cores, len(all_data_list) )
                
                batch = all_data_list[:num_to_process]
                
                # Remove the processed objects from the list
                all_data_list = all_data_list[num_to_process:]
                
                pool = multiprocessing.Pool(processes=num_to_process)
                    
                results = pool.map(process_light_curve, batch)
                pool.close()
                pool.join()
    
                for res in results:
                    if res is not None:
                        fout.write(res)
    
                counter += 1
    
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
