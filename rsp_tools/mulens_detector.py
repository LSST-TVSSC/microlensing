import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.lines as mlines
import itertools
from astropy.coordinates import SkyCoord
import astropy.units as u
import argparse

from lsst.rsp import get_tap_service
from lsst.daf.butler import Butler
import lsst.afw.display as afw_display
import lsst.sphgeom as sphgeom
import lsst.geom as geom
from lsst.utils.plotting import (get_multiband_plot_colors,
                                 get_multiband_plot_symbols,
                                 get_multiband_plot_linestyles)

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

    # Fetch a table of all Objects within this search radius
    query = """
            SELECT objectId, coord_ra, coord_dec,
                   u_psfMag, u_cModelMag, g_psfMag, g_cModelMag,
                   r_psfMag, r_cModelMag, i_psfMag, i_cModelMag,
                   z_psfMag, z_cModelMag, y_psfMag, y_cModelMag,
                   refExtendedness
            FROM dp1.Object
            WHERE CONTAINS(POINT('ICRS', coord_ra, coord_dec),
                  CIRCLE('ICRS', {}, {}, {})) = 1
            ORDER BY objectId
            """.format(ra_cen, dec_cen, radius)
    job = service.submit_job(query)
    job.run()
    job.wait(phases=['COMPLETED', 'ERROR'])
    print('Job phase is', job.phase)
    if job.phase == 'ERROR':
        job.raise_if_error()
    assert job.phase == 'COMPLETED'
    objtab = job.fetch_result().to_table()

    # Select objects with valid iband magnitudes
    mask = objtab['i_psfMag'] != np.nan
    objtab[mask]

    # For each selected target, retrieve the DP1 lightcurve and run the microlensing fitter
    source_radius = 0.5 # Arcseconds
    for source in objtab[mask]:
        query = "SELECT fsodo.diaObjectId, fsodo.coord_ra, fsodo.coord_dec, " \
                "fsodo.visit, fsodo.detector, fsodo.band, " \
                "fsodo.psfDiffFlux, fsodo.psfDiffFluxErr, " \
                "fsodo.psfFlux as psfFlux, fsodo.psfFluxErr, " \
                "vis.expMidptMJD " \
                "FROM dp1.ForcedSourceOnDiaObject as fsodo " \
                "JOIN dp1.Visit as vis ON vis.visit = fsodo.visit " \
                "WHERE CONTAINS (POINT('ICRS', coord_ra, coord_dec), " \
                "CIRCLE('ICRS', " + str(source.coord_ra) + ", " \
                + str(source.coord_dec) + ", " + str(source_radius) + ")) = 1 "
        job = service.submit_job(query)
        job.run()
        job.wait(phases=['COMPLETED', 'ERROR'])
        print('Job phase is', job.phase)
        if job.phase == 'ERROR':
            job.raise_if_error()
        assert job.phase == 'COMPLETED'

        # This returns a table of the lightcurves of nearby objects
        forced_sources = job.fetch_result().to_table()

        # KMTNet algorithm goes here!


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
