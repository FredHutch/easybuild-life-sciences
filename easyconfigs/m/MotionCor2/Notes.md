This file had to be hacked for the install part.

CUDA version is searched from the dependencies list, (this is bad)
The file name format changed.

# old
self.motioncor2_bin = 'MotionCor2_%s-Cuda%s' % (self.motioncor2_verstring, cuda_short_ver)

# new
self.motioncor2_bin = 'MotionCor2_%s_Cuda%s-%s' % (self.motioncor2_verstring, cuda_short_ver, '02-15-2020')

/app/software/EasyBuild/4.3.4/lib/python2.7/site-packages/easybuild/easyblocks/m/motioncor2.py
