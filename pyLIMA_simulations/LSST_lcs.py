import numpy as np
import matplotlib.pyplot as plt

import simulate_lc as smlc

time = np.arange(2461000,2462000,50)
ra = 270
dec = -28

#adjust tetelescope ZP, u,g,r,i,z,y
zps = [27.8,28,28.2,28.4,28,27.2]
source_mags = [24.8,23,22.2,23.2,23.8,24.2]
blend_mags =  [23.8,22.5,22,23.6,23,24.2]

u = smlc.simulate_one_LSST_band(time,band_name='u')
g = smlc.simulate_one_LSST_band(time,band_name='g')
r = smlc.simulate_one_LSST_band(time,band_name='r')
i = smlc.simulate_one_LSST_band(time,band_name='i')
z = smlc.simulate_one_LSST_band(time,band_name='z')
y = smlc.simulate_one_LSST_band(time,band_name='y')

ml_event =   smlc.simulate_a_microlensing_event(ra,dec,[u,g,r,i,z,y])
#lcs in fluxes!
lcs,ml_model,pyLIMA_parameters = smlc.simulate_lcs(ml_event,model_choice='USBL',source_mags=source_mags,blend_mags=blend_mags,zps=zps)

for ind,lc in enumerate(lcs):
    plt.scatter(time,zps[ind]-2.5*np.log10(lc))

plt.gca().invert_yaxis()
plt.show()

breakpoint()
