#import sys
#sys.path.append('/home/python_libs/')
#sys.path.append('/home/python_libs/lib64/python2.6/site-packages/')
#sys.path.append('/home/python/katcp-0.3.5/')
#sys.path.append('/home/python/construct-2.5.1/')
#sys.path.append('/home/python/construct-2.5.1/build')
#sys.path.append('/home/python/spead-0.5.0/build/')
#sys.path.append('/home/python/six-1.6.1/')
#sys.path.append('/home/python/corr-0.7.3/')
#sys.path.append('/home/python_libs/lib/python/')



#import numpy, scipy, katcp, corr,struct,time
#import matplotlib.pyplot as plt
#from operator import xor
#from optparse import OptionParser

import numpy 
import matplotlib.pyplot as plt #http://matplotlib.org/users/pyplot_tutorial.html
from matplotlib.pyplot import *
import os
#%%
#setup
path = os.getcwd()
f = 'files/test20s_0.dat'  #153 rows, 256 colums
bin_file = os.path.join(path, f)
tmp_ascii = f.split('.')[0] + '.txt'
print "Bin File Path: %s" % bin_file
print "Tmp Ascii File: %s" % tmp_ascii
navgd = 10;
print "navgd %d" %navgd
#%%
#spec = numpy.zeros((256)) #makes an zero array, but appends add to the zero array, so [0,0,0..][1,2,4..]
spec=[]
tmp = numpy.zeros(256)
#%%
#average the first 10 scans and plot
i=0
specave=numpy.zeros(256)
tmp = numpy.zeros(256)
infile = open(bin_file, "r")
while i < navgd: #only reads for navgs rows
	tmp = numpy.fromfile(infile, dtype='>I', count=256)
	print 'scan', i	
	specave = numpy.add(specave, tmp)  
	i=i+1

infile.close()
specave=specave/navgd  #averaged spectra

plt.figure()
plt.plot(specave[0:256]) #[0:65]... channel 64 and 191, 256 are bad
axis([0, 256, 1E6, 2E7])
xlabel('Index')
ylabel('Channel Count')
title('Channel Count by Channel')
#plt.annotate('local max', xytext=(64, 1.5E7)) #should work
plt.show()
#%%
#generate time base
t=130.072 #ms
i=0
spec=[]; time=[]
f = open(bin_file, "r")
#while (len(numpy.fromfile(f, dtype='>I', count=256)) != '\n' ): #does not seem to work, till scan77. must be a [] or EOL
for i in range (100000): #scans, each 0.130s. the longest is say 2h/0.130s = 55,384 
   try:
     tmp = numpy.fromfile(f, dtype='>I', count=256)
     #print tmp
     if len(tmp) == 0:
       break
     else:
       spec = numpy.append(spec, tmp)
       time.append(t*i/1000.0)
       #print 'scan#', i, 't=', time[i] 	
       i += 1
   except ValueError:
     pass
f.close()
#print '# of scans: ', len(spec)/256  #number of scans
#print len(time), min(time), max(time)
#%%
#generate rows. Each scan or row has 256 channels (3.90625 MHz each)
row=[]
scans = len(spec)/256
for i in range (int(scans)):
   row.append(spec[256*i:255+256*i])
   #print 'row#', i, 256*i, 255+256*i
#print len(row), len(row[0]), len(time), scans   
#%%
#put all the data into channels or frequency bins (3.90625 Mhz each)
ch=[]
for i in range (len(row[0])):
  tmp=[]
  for j in range (len(row)):
      tmp.append(row[j][i])
  ch.append(tmp)
#print len(ch), len(ch[0])
#%%
#write each channel to a tmp file so it can be read by allan.c program
for i in range (len(ch)): #256 channels
   f = open(tmp_ascii, 'w')  #keep overriding
   print 'ch: ', i
   for j in range (len(row)): 
       f.write('%8.5f %11.3f \n' %(time[j], ch[i][j] )) 
   f.close() 
   #read ch.txt and pipe to allan.c and append three output columns
#%%
#colorwheel
from itertools import permutations
from random import sample
Nlines = 200
color_lvl = 8
rgb = np.array(list(permutations(range(0,256,color_lvl),3)))/255.0
colors = sample(rgb,Nlines)
#%%
#plot all 256 channels vs time
print len(time), len(ch[0])
plt.figure()
#plt.plot(spec);
for i in range (0, len(time), 1):
  if ((i != 64) or (i !=191)):   
    plt.plot(time, ch[i], label='Ch', color=colors[i])
#legend = plt.legend(loc='upper center', shadow=True, fontsize='x-large')
axis([0, 20, 1E6, 2E7])
xlabel('time (s)')
ylabel('Channel Count')
title('Individual Channel Time series')
plt.show()
#plt.close() #close plot
#plt.semilogy(time[0:len(time)-1], ch[0], label='Model length'); #blue is default, there is one to many time points, so remove last time
#ax0.set_title('Set default color cycle to rgby')
#plt.semilogy(time[0:len(time)-1], ch[0], 'ro'); #there is one to many time points, so remove last time
