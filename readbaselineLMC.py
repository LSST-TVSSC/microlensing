import pandas as pd 
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u


year=float(365.24)
fov =float(3.5/2.0)## in degree    
###############################################################################################################################
###############################################################################################################################
#0               1       2          3                               5                    7                   9
#observationId,fieldRA,fieldDec,observationStartMJD,flush_by_mjd,visitExposureTime,band,filter,rotSkyPos,rotSkyPos_desired,
#numExposures,airmass,seeingFwhm500,seeingFwhmEff,seeingFwhmGeom,skyBrightness,night,slewTime,visitTime,slewDistance,
 #10           11                      13                 14           15       16     17          18        19              
#fiveSigmaDepth,altitude,azimuth,paraAngle,pseudoParaAngle,cloud,moonAlt,sunAlt,scheduler_note,target_name,
#    20           21       22         23        24           25     26    27      28               29        
#target_id,observationStartLST,rotTelPos,rotTelPos_backup,moonAz,sunAz,sunRA,sunDec,moonRA,moonDec,moonDistance,
#solarElong,moonPhase,cummTelAz,observation_reason,science_program,cloud_extinction
###############################################################################################################################
###############################################################################################################################
fil1=open("./counterOpsim.txt","w")
fil1.close()
nam1=['ID','RA', 'DEC', 't', 'filter', 'sig5']#6 
#      0     1     2     3      6        20
df= pd.read_csv('./output.csv',delimiter=',',header=None,skiprows=[0,],usecols=[0,1,2,3,6,20],names=nam1)
nm=int(len(df['RA']))

print("number of data points:  ",nm)
ra =np.zeros((nm))
dec=np.zeros((nm))
ra =df['RA']
dec=df['DEC']
tmin=float(np.min(df['t']))
print (df['RA'],    np.min(df['RA']),    np.max(df['RA']) )
print (df['DEC'],   np.min(df['DEC']),   np.max(df['DEC']))
print (df['t'],     np.min(df['t']),     np.max(df['t'])  )
print (df['filter'],np.min(df['filter']),np.max(df['filter']))
print (df['sig5'],  np.min(df['sig5']),  np.max(df['sig5']))
print ("**********************************************************************************")

###############################################################################################################################


'''
for i in range(801): 
    for j in range(81): 
        lon= float(-100.0+ 0.25*i)
        lat=  float(-10.0+ 0.25*j)
        fil1=open("./OPSIM/opm{0:.2f}_{1:.2f}.txt".format(lon, lat),"w")
        fil1.close()
'''        
       

###############################################################################################################################

    
nr=np.zeros((801,81))
nfil=int(0)    
for k in range(nm): 
    if(float(df['t'][k]-tmin)>0.0): 
        c=SkyCoord(ra=float(ra[k]), dec=float(dec[k]), unit=(u.degree ,  u.degree)  )
        l=float(c.galactic.l.degree)##[0, 360.0]
        b=float(c.galactic.b.degree)##[-90, 90]
        nfil=-1
        if(  df['filter'][k]=='u'):  nfil=0
        elif(df['filter'][k]=='g'):  nfil=1
        elif(df['filter'][k]=='r'):  nfil=2
        elif(df['filter'][k]=='i'):  nfil=3
        elif(df['filter'][k]=='z'):  nfil=4
        elif(df['filter'][k]=='y'):  nfil=5
        else: 
            print ("Big error , filter",  df['filter'][k], nfil)
            input("Enter a number  ")
        nok=0;  
        for i in range(801): 
            for j in range(81): 
                lon=float(-100.0+ 0.25*i)
                lat= float(-10.0+ 0.25*j)          
                if(lon<0.0): lon2=float(360.0+lon)## that is between [0, 360] degree
                else:        lon2=float(lon)
                dist=np.sqrt((lon2-l)**2.0+(lat-b)**2.0)##degree
                if(dist<fov or dist==fov):##degree
                    nok+=1; 
                    fil1=open("./OPSIM/opm{0:.2f}_{1:.2f}.txt".format(lon, lat),"a+")  
                    tst=np.array([ ra[k], dec[k], float(df['t'][k]-tmin) , int(nfil), float(df['sig5'][k])  ])##5
                    np.savetxt(fil1, tst.reshape(-1,5), fmt= "%.3f  %.3f  %.9f  %d   %.6f")
                    fil1.close()
                    nr[i,j]+=1
        print(ra[k], dec[k], float(df['t'][k]-tmin)/year ,   nfil, l, b, df['sig5'][k]   , nok ,    k) 
        print("*******************************************************************************")
                
print(nr) 






for i in range(801): 
    for j in range(81): 
        tst=np.array([i, j, nr[i,j] ])
        fil2=open("./counterOpsim.txt","a+")
        np.savetxt(fil2, tst.reshape(-1,3),fmt="%d  %d  %d")
        fil2.close()

###############################################################################################################################
###############################################################################################################################















