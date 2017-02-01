#! /bin/bash
# command for deploying proxmox container
if [[ -d sandbox ]]; then
  if [[ -z "$1" ]]; then
    echo "usage: prox.sh <newhostname>"
    exit 1
  fi
  prox new --mem 8G --disk 8G --cores 8 --no-bootstrap --runlist sandbox/easybuild.runlist $1
else
  echo "folder sandbox does not exist. Please switch to the root of the easybuild-life-sciences repos."
fi
