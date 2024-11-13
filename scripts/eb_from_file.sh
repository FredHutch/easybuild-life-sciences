#!/bin/bash
# Build software with EasyBuild
# <build.txt> containes a list of easyconfigs one per line. minus the ".eb" extension

MOD_DIR=/app/modules/all

function convertsecs() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 printf "%02d:%02d:%02d\n" $h $m $s
}

echo '==' `date +"%Y-%m-%d %H:%M"` START >built.log

for x in `cat build.txt`; do
    start_epoch=`date +%s`
    start_time=`date +"%Y-%m-%d %H:%M"`
    modpath=${x//://}.lua
    easyconfig=${x//:/-}.eb

    if [[ -f ${MOD_DIR}/${modpath} ]]; then
       echo Found: ${modpath}
    else
       echo "== ${start_time} Installing: ${easyconfig}"
       status=`eb easyconfig --robot easyconfigs 2>&1`
       elapsed=$(( `date +%s` - start_epoch ))
       if [[ $? == 0 ]]; then
           echo "==" `date +"%Y-%m-%d %H:%M"` Built: ${x} elapsed: $(convertsecs $elapsed)  >> built.log
       else
           echo "== Error" $easyconfig >> build.err 
           echo $status >> build.err 
       fi
    fi
done
echo "==" `date +"%Y-%m-%d %H:%M"` End >>build.err
echo "==" `date +"%Y-%m-%d %H:%M"` End >>built.log

