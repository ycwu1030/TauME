#!/bin/bash

datafolder="/Users/ycwu/Workingspace/MC-Generators/MG5_aMC_v3_5_4/bin/data_pion_tauEDM"
delphesdir="/Users/ycwu/Workingspace/Misc/delphes/build/readers"
delphescard="/Users/ycwu/Library/CloudStorage/OneDrive-Personal/Projects/000.Codings/TauOO/scripts/gen_card.tcl"
rootdir="/Users/ycwu/Library/CloudStorage/OneDrive-Personal/Projects/000.Codings/TauOO/test"
rootnameprefix="delphes_events_pion"
n=6
dflist=(0.0 100.0 -100.0 200.0 -200.0 500.0 -500.0)
dfnames=(0x0 100x0 m100x0 200x0 m200x0 500x0 m500x0)
# n=8
# dflist=(0.0 0.2 -0.2 0.4 -0.4 0.8 -0.8 1.0 -1.0 1.5 -1.5 2.0 -2.0 5.0 -5.0 10.0 -10.0 20.0 -20.0 50.0 -50.0 100.0 -100.0 200.0 -200.0 500.0 -500.0)
# dfnames=(0x0 0x2 m0x2 0x4 m0x4 0x8 m0x8 1x0 m1x0 1x5 m1x5 2x0 m2x0 5x0 m5x0 10x0 m10x0 20x0 m20x0 50x0 m50x0 100x0 m100x0 200x0 m200x0 500x0 m500x0)

for id in $(seq 0 1 $n)
  do
    df=${dflist[$id]}
    dfn=${dfnames[$id]}
    datadir=${datafolder}/${dfn}
    files=$(ls ${datadir})
    fileid=0
    rootid=0
    filelist=""
    for file in $files
      do
        gunzip -k ${datadir}/${file}
        filenoext=${file%.*}
        filelist="${filelist} ${datadir}/${filenoext}"
        fileid=$[$fileid+1]
        if [ $fileid -eq 10 ]; then
          ${delphesdir}/DelphesLHEF ${delphescard} ${rootdir}/${rootnameprefix}_${dfn}_${rootid}.root $filelist
          rm ${datadir}/*.lhe
          echo $filelist
          rootid=$[$rootid+1]
          fileid=0
          filelist=""
        fi
      done
      if [ $fileid -gt 0 ]; then
        ${delphesdir}/DelphesLHEF ${delphescard} ${rootdir}/${rootnameprefix}_${dfn}_${rootid}.root $filelist
          rm ${datadir}/*.lhe
          echo $filelist
          rootid=$[$rootid+1]
          fileid=0
          filelist=""
      fi
done
