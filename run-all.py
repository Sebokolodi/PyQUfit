import os
import glob
import numpy



Q = glob.glob('CYGDATA/Q/q*.txt')
U = glob.glob('CYGDATA/U/u*.txt')
I = glob.glob('CYGDATA/I/i*.txt')

list.sort(Q, key=lambda x: int(x.split('-')[-1].split('.txt')[0]))
list.sort(U, key=lambda x: int(x.split('-')[-1].split('.txt')[0]))
list.sort(I, key=lambda x: int(x.split('-')[-1].split('.txt')[0]))

#f = 'CYGDATA/FREQ/f.txt'
#nq = 'CYGDATA/NOISE/noise-q.txt'
#nu = 'CYGDATA/NOISE/noise-u.txt'
#ni = 'CYGDATA/NOISE/noise-i.txt'

Q = ','.join(Q[0:5])
U = ','.join(U[0:5])
I = ','.join(I[0:5])
#print('%r'%Q)


os.system('./qufit.py -json config.json -f CYGDATA/FREQ/f.txt -nq CYGDATA/NOISE/noise-q.txt -nu CYGDATA/NOISE/noise-u.txt -ni \
          CYGDATA/NOISE/noise-i.txt -nparam 5 -plot True -out TEST -pref FIT-QU -q %r -u %r -i %r'%(Q, U, I))





