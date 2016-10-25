#! /bin/bash
# command for deploying proxmox container
prox --mem 8G --disk 8G --cores 8 --no-bootstrap --runlist easybuild.runlist new myhost

