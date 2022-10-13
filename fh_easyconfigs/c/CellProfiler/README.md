## CellProfiler Build Notes 4.1.1

Add appMenu to OS during build. Installed on Rhino's
apt-get install appmenu-gtk3-module
sudo apt-get install appmenu-gtk2-module appmenu-gtk3-module

# this worked once, rebuild wx with gtk-3 installed then it failed.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/gtk-3.0/modules
