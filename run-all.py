import os
import glob
import numpy


#source_name = '0004-5122'
#source_name = '2058-4905'
source_name = '2103-5225'
model = 'm11'
nparam =  6

input_dir = '../data/%s'%source_name + '/'
output_dir = '../test/%s'%source_name + '/'


for i in range(1, 3):
    data = input_dir + '%s.lobe%d.allflux.RMtools.txt'%(source_name,i)

    outname = source_name + 'lobe%d'%i + '%s'%model
    os.system('./qufit.py -json config.json -nparam %d -plot True -phi 2000 -dphi 5 -out %s -pref %s -indata %r' %(nparam, output_dir, outname, data))




