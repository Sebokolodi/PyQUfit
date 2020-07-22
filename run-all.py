import os
import glob
import numpy



data = glob.glob('CYGDATA/DATA/data-*.txt')
noise = 'CYGDATA/NOISE/noise-data.txt'
freq =  'CYGDATA/FREQ/f.txt' 


list.sort(data, key=lambda x: int(x.split('-')[-1].split('.txt')[0]))
data = ','.join(data[0:2])


os.system('./qufit.py -json config.json -f %r -noise %r -nparam 5 -plot True -out TEST -pref FIT-QU -indata %r'
         %(freq, noise, data))








