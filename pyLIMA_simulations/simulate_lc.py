import numpy as np
import pyLIMA.event
import pyLIMA.models.PSPL_model
import pyLIMA.models.USBL_model
import pyLIMA.simulations.simulator



def simulate_one_LSST_band(time,band_name='g'):

    telo = pyLIMA.simulations.simulator.simulate_a_telescope("Rubin_"+band_name, timestamps=time, location="Earth", astrometry=False,)
    
    return telo
    
def simulate_a_microlensing_event(ra,dec,telescopes_list):

    ml_event = pyLIMA.event.Event(ra=ra, dec=dec)
    ml_event.name = "Rubin_" + str(ra) + "_" + str(dec)
    
    for telo in telescopes_list:
        ml_event.telescopes.append(telo)

    return ml_event
    
def refine_microlensing_parameters(themodel, microlensing_parameters, time):
   

    # Ensure t0 is within the observations
    t0min, t0max = np.percentile(time, [20, 80])

    new_ml_parameters = np.copy(microlensing_parameters)
    new_ml_parameters[0] = np.random.uniform(t0min, t0max)       

    new_tE = 10 ** (np.random.normal(1.2, 0.3))
    new_ml_parameters[2] = new_tE

    if themodel.model_type()=='USBL':

        # Rho 
        new_ml_parameters[3] = 10 ** np.random.uniform(-4.5, -1.30)

        # Impact parameter u0 -- the following ensures some caustic signal
        new_ml_parameters[1] = np.random.uniform(-1 * new_ml_parameters[3], 1 * new_ml_parameters[3])

        # Separation 
        new_ml_parameters[4] = 10 ** np.random.uniform(-0.5, 0.5)

        # Mass ratio
        new_ml_parameters[5] = np.random.uniform(10**-5,1)


    return new_ml_parameters

def refine_flux_parameters(source_mag,blend_mag,zp=27.85):

    f_source = 10 ** ((zp - source_mag) / 2.5)

    fblend = 10 ** ((zp - blend_magg) / 2.5)

    return f_source,f_blend
    
def simulate_lcs(ml_event,model_choice='PSPL',source_mags=[],blend_mags = [],zps=[]):

    if model_choice == 'PSPL':

        themodel = pyLIMA.models.PSPL_model.PSPLmodel(ml_event,blend_flux_parameter='fblend',)

    if model_choice == 'USBL':

        choice = np.random.choice(["central_caustic", "second_caustic", "third_caustic"])
        themodel = pyLIMA.models.USBL_model.USBLmodel(ml_event, origin=[choice, [0, 0]],blend_flux_parameter='fblend')
    
    the_parameters = (pyLIMA.simulations.simulator.simulate_microlensing_model_parameters(themodel))
    all_time = []
    for tel in ml_event.telescopes:
       all_time.append(tel.lightcurve['time'].value)
       
       
    all_time = np.unique(all_time)    
    the_parameters = refine_microlensing_parameters(themodel, the_parameters, all_time)
    pyLIMA_parameters = themodel.compute_pyLIMA_parameters(the_parameters)

    lcs = []
    
    for ind,telo in enumerate(ml_event.telescopes):
    
        try:
            new_fs, new_fb = refine_flux_parameters(source_mags[ind],blend_mags[ind],zp=zps[ind])
            pyLIMA_parameters['fsource_'+str(telo.name)] = new_fs
            pyLIMA_parameters['fblend_'+str(telo.name)] = new_fb
        except:
            pass

        magnification = themodel.model_magnification(telo, pyLIMA_parameters)
        flux = pyLIMA_parameters['fsource_'+str(telo.name)] * magnification +  pyLIMA_parameters['fblend_'+str(telo.name)]
        
        lcs.append(flux)
        
        
    return lcs,themodel,pyLIMA_parameters

