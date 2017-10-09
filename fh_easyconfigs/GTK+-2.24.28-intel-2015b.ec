easyblock = 'ConfigureMake'

name = 'gtk+'
version = '2.24.28'

homepage = 'https://developer.gnome.org/gtk+/stable/'
description = """
 The GTK+ 2 package contains libraries used for creating graphical user interfaces for applications. 
"""

toolchain = {'name': 'intel', 'version': '2015b'}

source_urls = ['http://ftp.gnome.org/pub/GNOME/sources/%(namelower)s/%(version_major_minor)s']
sources = [SOURCELOWER_TAR_XZ]

# 2016-01-26 CJB Bumped Pango from 1.37.1 to 1.38.1 due to dependencies conflicts.
dependencies = [
    ('ATK', '2.16.0'),
    ('gdk-pixbuf', '2.31.4'),
    ('Pango', '1.38.1'),
]

moduleclass = 'vis'
