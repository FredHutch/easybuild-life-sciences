#!/usr/bin/env python3

# very crud methon to update EasyConfig dependencie list. 
# Create Module list like this:
# find /app/software/EasyBuild/4.9.4/easybuild/easyconfigs/ -name \*\.eb -print | sed 's;^/.*/;;' >Mod_List
# cat Mod_list | sed 's;^/.*/;;' | grep -e 2023a -e 12\.3\.0 >2023a.list 

import os
import sys
from os import listdir
from os.path import isfile, join

from toolchain import toolchain

tc_toolchain = '2023a'

mod_list = [ 'matplotlib', 'Biopython', 'Pysam', 'pybedtools', 'PyYAML', 'PyTables',
'h5py', 'numba', 'awscliv2', 'Pillow', 'scikit-learn', 'scikit-build', 'scikit-bio',
'scikit-image', 'scikit-optimize', 'Porechop', 'python-igraph', 'statsmodels', 'Blosc',
'dask', 'PostgreSQL', 'typing-extensions', 'Arrow', 'ICU', 'Pandoc', 'FreeTDS',
'OpenJPEG', 'OpenBLAS', 'Tk', 'libxml2', 'libxslt', 'libffi', 'Qt5', 'cURL',
'libGLU', 'Mesa', 'netCDF', 'igraph', 'snappy', 'freetype',
]
module_list = {}

def create_module_list(tc):
    """ traverse the EasyBuild easyconfigs
    """
    eb_install_path = '/app/software/EasyBuild'
    eb_ver = '4.9.4'
    eb_easyconfig_path = os.path.join(eb_install_path, eb_ver, 'easybuild/easyconfigs')
    if not os.path.isdir(eb_easyconfig_path) :
       print(f'could not find easyconfigs: {eb_easyconfig_path}')
       sys.exit(1)
    for top_dir in list(map(chr, range(ord('a'), ord('z')+1))):
        for mod_dir in os.listdir(os.path.join(eb_easyconfig_path, top_dir)):
            for eb_config in os.listdir(os.path.join(eb_easyconfig_path, top_dir, mod_dir)):
                if tc.tc_filter(eb_config):
                    mod_suffix = eb_config.removeprefix(mod_dir)[1:][:-3]
                    mod_version = tc.tc_trim(mod_suffix)
                    if mod_dir in module_list:
                        module_list[mod_dir].append(mod_version)
                    else:
                        module_list[mod_dir] = [mod_version]

if __name__ == '__main__':
   tc = toolchain('2023a')
   create_module_list(tc)
   for mod in mod_list:
      if mod in module_list:
          for mod_ver in module_list[mod]:
              print(f"    ('{mod}', '{mod_ver}'),")
      else:
         print(f"    ('{mod}', 'Not Found'),")

