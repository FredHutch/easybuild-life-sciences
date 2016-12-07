#! /bin/bash
# command for deploying proxmox container
if [[ -z "$1" ]]; then
  echo "usage: prox.sh <newhostname>"
  exit 1
fi
prox new --mem 8G --disk 8G --cores 8 --no-bootstrap --runlist easybuild.runlist $1

