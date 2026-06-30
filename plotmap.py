import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import math
import pylab
from matplotlib import cm
from matplotlib import colors
cmap=plt.get_cmap('viridis')
mpl.rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rc('text', usetex=True)
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams["font.size"] = 11.5
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Computer Modern Sans"]
mpl.rcParams["text.usetex"] = True
mpl.rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
mpl.rcParams['text.usetex'] = True
from astropy.coordinates import SkyCoord
import astropy.units as u



def tickfun( start, x, dd0):
    return( start+x*dd0)
    

######################################################################################################

tt=9
v =np.zeros((tt))
nx=801
ny=81
ddeg=0.25

par=np.zeros((nx*ny,3))
par=np.loadtxt("./counterOpsim.txt")

Mapp=np.zeros((ny, nx))
k=0; 
for i in range(801): 
    for j in range(81): 
        Mapp[j,i]= np.log10(par[k,2]+1.0)
        k+=1
print(Mapp, np.min(Mapp) , np.max(Mapp) )

##################################################################################################

plt.cla()
plt.clf()
fig=plt.figure(figsize=(10,3.5),facecolor='w')
ax=fig.add_subplot(111)
plt.imshow(Mapp,cmap=cmap,extent=(0.0,801.0,0.0,81.0),interpolation='nearest',aspect='auto', origin='lower')
plt.clim()
minn=np.min(Mapp)
maxx=np.max(Mapp)
step=float((maxx-minn)/(tt-1.0));
for m in range(tt):
    v[m]=round(float(minn+m*step),1)
cbar=plt.colorbar(orientation='horizontal',shrink=0.8,pad=0.1,ticks=v)
cbar.ax.tick_params(labelsize=16)
plt.clim(v[0]-0.05*step,v[tt-1]+0.05*step)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.title(r"$\log_{10}[\rm{No.}~\rm{of}~\rm{visits}]$", fontsize=16)
plt.xlim(nx*0.0 , nx*1.0)
plt.ylim(ny*0.0 , ny*1.0)
ticc=np.array([nx*0.1, nx*0.3, nx*0.5, nx*0.7, nx*0.9 ])
ax.set_xticks(ticc,labels=[int(j) for j in tickfun(-100.0, ticc,ddeg) ])
ticc=np.array([ny*0.1, ny*0.3, ny*0.5, ny*0.7, ny*0.9 ])
ax.set_yticks(ticc,labels=[int(j) for j in tickfun(-10.0,  ticc,ddeg) ])
#ax.invert_xaxis()
#ax.set_aspect('equal', adjustable='box')
plt.xlabel(r"$\rm{l}(\rm{deg})$",fontsize=22,labelpad=0.05)
plt.ylabel(r"$\rm{b}(\rm{deg})$",fontsize=22,labelpad=0.05)
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./visit.jpg", dpi=200)
    
    
    
####################################################################################################################3        
    
coor=np.zeros((nx*ny,2))    
k=0    
for i in range(nx): 
    for j in range(ny): 
        lon=float(-100.0+ 0.25*i)
        lat=float( -10.0+ 0.25*j) 
        if(lon<0.0): lon2=float(lon+360.0)   
        else:        lon2=lon; 
        c=SkyCoord(l=float(lon2), b=float(lat), unit=(u.degree ,  u.degree), frame='galactic')
        ceq=c.transform_to('icrs')
        RA= float(ceq.ra.degree)
        DEC=float(ceq.dec.degree)
        coor[k,0], coor[k,1]= RA, DEC
        k+=1
       
       
####################################################################################################################3           
        
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8,6))
ax1=fig.add_subplot(111)
plt.scatter(coor[:,0], coor[:,1], marker="o", s=7.2, color='g', alpha=0.7)##"o", markersize=3, color=col[tst[i,4]])
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
ax1.set_aspect('equal', adjustable='box')
#plt.xlim(RA0 ,  RA1)
#plt.ylim(DEC0 , DEC1)
plt.xlabel(r"$\rm{Right}~\rm{Ascention}(\rm{deg})$",fontsize=20,labelpad=0.05)
plt.ylabel(r"$\rm{Declination}(\rm{deg})$",fontsize=20,labelpad=0.05)
fig=plt.gcf()
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./RADECGB.jpg" , dpi=200)    
    
####################################################################################################################3    
       
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8,6))
ax1=fig.add_subplot(111)
scatt=plt.scatter(coor[:,0], coor[:,1], marker="o", s=7.2, c=np.log10(par[:,2]+1.0),cmap=cmap, alpha=0.7)
plt.colorbar(scatt, label=r"$\log_{10}[\rm{No.}~\rm{of}~\rm{visits}]$")
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
#ax1.set_aspect('equal', adjustable='box')
plt.xlabel(r"$\rm{Right}~\rm{Ascention}(\rm{deg})$",fontsize=20,labelpad=0.05)
plt.ylabel(r"$\rm{Declination}(\rm{deg})$",fontsize=20,labelpad=0.05)
fig=plt.gcf()
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./RADECGB2.jpg" , dpi=200)    
    
    
####################################################################################################################3        
        
    
    
    
