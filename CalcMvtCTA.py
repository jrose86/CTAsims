### Calculate Optical Path Delays, Projected Baselines, and Noise Models for all pairs of telescopes
### Made by Sahar Nikkhah
### Modified by Mike Lisa and Mackenzie Scott
### Modified by Josie Rose to add uv

### josie adapted 3 nov 2023 to run with input parameters from analysis macro
### 18 nov 2023 josie - evaluating baseline, opd at each frame
### 21 nov 2023 josie - print to just ONE file per pair, add full uv calcs, clean up unnecessary parts
### 22 nov 2023 josie - calculations with ONLY uvw coords geometry for everything
### 23 dec 2023 josie - adapted for CTA coordinates to simulate baselines and uv coverage

### Import Packages
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astropy.coordinates import AltAz
from astropy.coordinates import SkyCoord
from astropy.coordinates import ITRS
import numpy as np
#import matplotlib.pyplot as plt
import math
from datetime import datetime
import subprocess
import sys # allows to read in args when executing script

# read in parameters from command line - for run length, source, and local start time
longestRun = sys.argv[1]
nframes = sys.argv[2]
frameSize = sys.argv[3]
source1 = sys.argv[4] # needs 2 inputs since the source name has a space
source2 = sys.argv[5]
locStartDate = sys.argv[6] # again 2 inputs bc of space
locStartTime = sys.argv[7]

# just a check that everything is read in as it should be
print("INPUTS longest run:", longestRun, " #frames:", nframes, " frame size:", frameSize, " source:", source1, source2, " start date, time:", locStartDate, locStartTime)

# Configure source and time from input parameters
TargetName = source1 + ' ' + source2
start_observing_time = Time(locStartDate + ' ' + locStartTime)
start_observing_time += TimeDelta(25200,format='sec') # Convert to UTC time
#dt = TimeDelta(float(Window)*int(FrameLong), format='sec')
dt = TimeDelta(float(longestRun), format='sec')
end_observing_time = start_observing_time + dt
Target = SkyCoord.from_name(TargetName)

print()
print(TargetName, "target coords", Target)

### Location of VERITAS, telescopes, and Length of cables
Veritas      = EarthLocation(lat='31.675', lon='-110.952', height='1270')  # location of veritas
#from cta website, north array coords: Longitude: 17º 53′ 31.218″ West   Latitude: 28º 45′ 43.7904″ North, from google search of MAGIC, alt = 2187m above sea level
# south array coords: Latitude: 24º 41′ 0.34″ South   Longitude: 70º 18′ 58.84″ West, alt = 2640m above sea level (from whatismyelevation.com)
CTAnorth     = EarthLocation(lat='28.762', lon='17.892', height='2187')
CTAsouth     = EarthLocation(lat='-24.683', lon='70.316', height='2640')
# telLocs      = np.array([[135.48, -8.61, 12.23], [44.1, -47.7, 4.4], [29.4, 60.1, 9.8], [-35.9, 11.3, 7.0]]) # old McGill/ASIIP coordinates
# telLocs      = np.array([[135.49, -8.72, 7.23], [45.145, -49.115, -0.94], [28.86, 60.73, 4.51], [-36.112, 11.696, 1.63]]) # these are the 2011 database coordinates
telLocs        = np.array([[135.48, -8.61, 12.23], [44.836, -49.601, 5.102], [29.335, 60.022, 10.636], [-35.885, 11.742, 6.417]]) # measured by Dave 2023, to cm precision

northCoords = np.array([[-205,-30,0], [-175,160,0], [-75,-180,0], [-50,-75,0], [-70,35,0], [0,0,0], [65,-35,0], [50,70,0], [-5,215,0], [170,220,0], [195,75,0], [310,-20,0], [195,-145,0]])


CableDelays  = [676.8e-9, 585.0e-9, 955.0e-9, 1063.7e-9] #Cable lengths for all 4 telescopes
bucketLength = 4e-9
lat = 0.553 # latitude of telescopes in rad for uv

# There are six PAIRS of telescopes, and I index them as follows:
# index 0=(1,2); 1=(1,3); 2=(1,4); 3=(2,3); 4=(2,4); 5=(3,4);

# --> Important: the first telescope of a pair is considered to be the origin, so these are the angles relative to that
Dvec    = np.empty([6,3])  # The Distance between the telescopes as a 3D vector
Tpair   = np.empty([6,2])  # identifies which is the first and the second telescope

##Determine the distance and angle between all 6 pairs of telescopes
iwhich = 0
for t1 in range(0,4):
    for t2 in range(t1+1,4):
        Tpair[iwhich,0] = t1
        Tpair[iwhich,1] = t2
        for icomp in range(0,3):
            Dvec[iwhich,icomp] = telLocs[t2,icomp] - telLocs[t1,icomp]
        iwhich += 1

print("Dvec: " , Dvec)
print()

# placeholder arrays for everything we're going to calculate
times = np.empty([0]) # time of day (in seconds)
frames = np.empty([0]) # each frame we calculate baseline, opd, etc for
u = [np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0])] # u coord, E/W proj plane of baseline
v = [np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0])] # v coord, N/S proj plane of baseline
w = [np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0])] # w coord, separation btwn telescopes in direction of light (same as opd)

# =========================== Iterate calculations over the whole run ==============================================================================

observing_time = start_observing_time
frame = 1
while observing_time < end_observing_time:
    #Configure altitude and azimuth of star at given time
    ITRScoords = Target.transform_to(ITRS(obstime=observing_time))
    hrang = Veritas.lon.rad - ITRScoords.spherical.lon.rad
    if abs(hrang) > 2*math.pi:
        if hrang > 2*math.pi:
            hrang = hrang - (2*math.pi)
        else:
            hrang = hrang + (2*math.pi)
    decl = ITRScoords.spherical.lat.rad

    for iwhich in range(0,6):
        # calculate uv coords
        u[iwhich] = np.append(u[iwhich], (Dvec[iwhich,0]*math.cos(hrang)) + (Dvec[iwhich,1]*-1*math.sin(lat)*math.sin(hrang)) + (Dvec[iwhich,2]*math.cos(lat)*math.sin(hrang)))
        v[iwhich] = np.append(v[iwhich], (Dvec[iwhich,0]*math.sin(hrang)*math.sin(decl)) + (Dvec[iwhich,1]*math.sin(lat)*math.cos(hrang)*math.sin(decl) + Dvec[iwhich,1]*math.cos(lat)*math.cos(decl)) + (Dvec[iwhich,2]*-1*math.cos(lat)*math.cos(hrang)*math.sin(decl) + Dvec[iwhich,2]*math.sin(lat)*math.cos(decl)))
        w[iwhich] = np.append(w[iwhich], ((Dvec[iwhich,0]*-1*math.sin(hrang)*math.cos(decl)) + (Dvec[iwhich,1]*-1*math.sin(lat)*math.cos(hrang)*math.cos(decl) + Dvec[iwhich,1]*math.cos(lat)*math.sin(decl)) + (Dvec[iwhich,2]*math.cos(lat)*math.cos(hrang)*math.cos(decl) + Dvec[iwhich,2]*math.sin(lat)*math.sin(decl)))/3e8)

    sectime = observing_time
    sectime.format = 'cxcsec'
    sectime.out_subfmt = 'decimal'
    times = np.append(times,sectime.value - start_observing_time.value)
    observing_time+=TimeDelta(frameSize,format='sec') # evaluate for each frame!
    frames = np.append(frames, frame)
    frame+=1

#print("check frames: ", frames)
### The pairs and time will be used in plot titles
pairLabel = ['T1T2','T1T3','T1T4','T2T3','T2T4','T3T4']
ttt = Time(start_observing_time, format='iso', out_subfmt='date_hms')
timeString = ttt.value

# ===================== save text files of data =====================================================================

# single new output file replaces all the old ones, but we still have one file per telescope pair
for iwhich in range(0,6):
    print("saving file for", pairLabel[iwhich])
    outfile = open('delays/pyinfo'+pairLabel[iwhich]+'.txt','w');
    #outfile.write('frame#   time in run (s)   u coord   v coord   baseline(m)   opd(ns)\n') # printing this line messes up reading in the TNtuple
    outfile.write('0   0   0   0   0   0\n'); # this seems silly but you MUST do it for the TNtuple not to skip the first row (has to be floats and not strings though?? -- needs more work!)
    for ipt in range(0, times.size):
        outfile.write('%d   %f   %f   %f   %f   %f\n' % (frames[ipt], times[ipt], u[iwhich][ipt], v[iwhich][ipt], math.sqrt((u[iwhich][ipt]*u[iwhich][ipt])+(v[iwhich][ipt]*v[iwhich][ipt])), -1.0e9*w[iwhich][ipt]))
    outfile.close()

print("python all done!")
