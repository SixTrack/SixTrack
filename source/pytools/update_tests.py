import os
from os.path import isfile, join
from shutil import copyfile


#testToUpdateArry = ['beam-HO_6D-newstyle','beam-HO_6D-oldstyle', 'beam-HO_6D-simple-oldstyle', 'beam-HO_6D-simple-newstyle', 'beam-HO_6D-simple-newstyle-orbit',
#'beam-HO_LR-newstyle','beam-HO_LR-oldstyle','beam-HO_LR-ibbc-oldstyle', 'beam-HO_LR-ibbc-newstyle', 'scatter_collimation','beam-HO_6D-simple-newstyle'    ]
testToUpdateArry = ['beam-HO_6D-simple-newstyle-coupling']
singeltr = '../../build/SixTrack_cmakesix_BUILD_TESTING_defaultcompiler_defaultbuildtype/SixTrack_50002_crlibm_rn_Linux_gfortran_static_x86_64_64bit_double'
no_stf = '../../build/SixTrack_cmakesix_BUILD_TESTING_NOSTF_defaultcompiler_defaultbuildtype/SixTrack_50002_crlibm_rn_nostf_Linux_gfortran_static_x86_64_64bit_double'

for testToUpdate in testToUpdateArry:
	os.chdir('../../test/'+testToUpdate)
	os.system(singeltr)
	os.system(no_stf)
	onlyfiles = [f for f in os.listdir('.') if isfile(join('.', f))]
	print(onlyfiles)
	for f in onlyfiles:
		if(f.endswith('.canonical')):
			fileName = f[:-10]
			copyfile(fileName, f)
			