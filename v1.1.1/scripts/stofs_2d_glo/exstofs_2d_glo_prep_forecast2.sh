#!/bin/bash
###############################################################################
#                                                                             #
# This script is the preprocessor for the STOFS that runs under the ADCIRC    #
# model. This script is the tide-only forecast for the STOFS.                 #
# It sets some shell script variables for export to child scripts             #
# and copies some generally used files to the work directory.                 #
# After this the actual preprocessing is performed by the following scripts:  #
#                                                                             #
# Remarks :                                                                   #
#                                                                             #
#                                                                 Dec, 2011   #
#                                                                 May, 2016   #
#                                                                 Jun, 2020   #
#                                                                 Apr, 2021   #
#                                                                 Jun, 2022   #
#                                                                             #
###############################################################################
# Start of stofs_2d_glo_prep_forecast2.sh.ecf script ------------------------ #
# 0.  Preparation
# 0.a Basic modes of operation

  seton='-xa'
  setoff='+xa'
  set $seton

  msg="Starting stofs_2d_glo_prep_forecast2 scripts"
  echo "$msg"
  postmsg "$jlogfile" "$msg"

  spinh=162
  wndh=3
  nowh=6
  lsth=180
  export ncpu=$NCPU

# --------------------------------------------------------------------------- #
# 1.  Set times
# 1.a Set all necessary times
#     YMDH     :  current time cycle in yyyymmddhh format
#     time_beg :  begin time of run (normally, -6 hour nowcast or the most recent hotstart time)
#     time_now :  current time
#     time_end :  ending time of run ($lsth hour forecast)

  export date=$PDY
  export YMDH=${PDY}${cyc}
  export nback=20

  time_now=$YMDH
  time_end=`$NDATE $lsth $YMDH`

# --------------------------------------------------------------------------- #
# 2.  Get hotstart and prep files 

  if [[ ! -f $COMOUTrerun/${RUN}_prep.68.nc ]]; then
     cpreq $COMOUT/${RUN}.${cycle}.hotstart ${time_now}.hotstart
  else
     cpreq $COMOUTrerun/${RUN}_prep.68.nc ${time_now}.hotstart
  fi 
  cpreq $COMOUTrerun/${RUN}_prep.61.nc fort.61.nc
  cpreq $COMOUTrerun/${RUN}_prep.63.nc fort.63.nc

# --------------------------------------------------------------------------- #
# 3.  Copy stofs_2d_glo partmesh, grid, and template files

  ln -s $FIXstofs/${RUN}_attr fort.13
  ln -s $FIXstofs/${RUN}_grid fort.14
  ln -s $FIXstofs/${RUN}_body fort.24
  ln -s $FIXstofs/${RUN}_rotm fort.rotm
  ln -s $FIXstofs/${RUN}_elev_stat elev_stat.151
  ln -s $FIXstofs/${RUN}_elev_stat vel_stat.151
  cpreq $FIXstofs/${RUN}_prep.15 ${RUN}_fort.15

# 3.a Copy stofs_2d_glo_nod_equi file and tar file

  if [[ -f $COMGES/${RUN}_nod_equi ]]; then
     cpreq $COMGES/${RUN}_nod_equi .
  fi
  if [[ -f $COMGES/${RUN}_${ncpu}.tar.gz ]]; then
     cpreq $COMGES/${RUN}_${ncpu}.tar.gz .
     tar xvzf ${RUN}_${ncpu}.tar.gz
     export err=$?; err_chk
  fi

# --------------------------------------------------------------------------- #
# 4.  Set time and file for tide-only forecast
# 4.a Set parameters for tide-only forecast input file

  ncdump -v time ${time_now}.hotstart > hotstart.out
  export err=$?; err_chk
  time_hotstart=$(grep 'time = [0-9]' hotstart.out | awk '{print $3}')
  rm hotstart.out

  if [ `expr $(echo "scale=0; $time_hotstart/($wndh*3600)" | bc) % 2` = 0 ]; then
     ihot=568
     cpreq ${time_now}.hotstart fort.68.nc
  else
     ihot=567
     cpreq ${time_now}.hotstart fort.67.nc
  fi

  fcstd=$(echo "scale=5; ($time_hotstart)/86400" | bc)
  rnday=$(echo "scale=5; $fcstd+$lsth/48" | bc)
  nout=5
  touts=$(echo "scale=5; $rnday-($nowh+$lsth)/24" | bc)
  toutf=$rnday
  nhstar=5
  nhsinc=1800

  exec 5<&0 < ${RUN}_nod_equi
      read hh dd mm yyyy
      mm=$(printf "%02d" $mm)
      dd=$(printf "%02d" $dd)
      hh=$(printf "%02d" $hh)
      read con1 fft1 facet1
      read con2 fft2 facet2
      read con3 fft3 facet3
      read con4 fft4 facet4
      read con5 fft5 facet5
      read con6 fft6 facet6
      read con7 fft7 facet7
      read con8 fft8 facet8

  sed -e "s/cycle/$time_now/g" \
      -e "s/ihot/$ihot/g" \
      -e "s/rnday/$rnday/g" \
      -e "s/fft1/$fft1/g" -e "s/facet1/$facet1/g" \
      -e "s/fft2/$fft2/g" -e "s/facet2/$facet2/g" \
      -e "s/fft3/$fft3/g" -e "s/facet3/$facet3/g" \
      -e "s/fft4/$fft4/g" -e "s/facet4/$facet4/g" \
      -e "s/fft5/$fft5/g" -e "s/facet5/$facet5/g" \
      -e "s/fft6/$fft6/g" -e "s/facet6/$facet6/g" \
      -e "s/fft7/$fft7/g" -e "s/facet7/$facet7/g" \
      -e "s/fft8/$fft8/g" -e "s/facet8/$facet8/g" \
      -e "s/nout/$nout/g" \
      -e "s/touts/$touts/g" -e "s/toutf/$toutf/g" \
      -e "s/nhstar/$nhstar/g" -e "s/nhsinc/$nhsinc/g" \
      -e "s/hh/$hh/g" -e "s/dd/$dd/g" \
      -e "s/mm/$mm/g" -e "s/yyyy/$yyyy/g" \
                ${RUN}_fort.15 | \
  sed -n "/DUMMY/!p" > fort.15
  rm ${RUN}_fort.15

  if [ ! -f fort.15 ]; then
     echo "FATAL ERROR: Tide-only forecast input file did not exist"
     err_exit
  else
     echo "Tide-only forecast input file existed"
  fi
 
# --------------------------------------------------------------------------- #
# 5.  Execute stofs_2d_glo_adcprep for tide-only forecast

  export pgm="stofs_2d_glo_adcprep for tide-only forecast"
  . prep_step
  startmsg
  if [[ ! -s ${RUN}_${ncpu}.tar.gz ]]; then
     mpiexec -n 1 -ppn 1 $EXECstofs/${RUN}_adcprep --np $ncpu --partmesh >> $pgmout 2>errfile
     export err=$?; err_chk
     mpiexec -n 1 -ppn 1 $EXECstofs/${RUN}_adcprep --np $ncpu --prepall >> $pgmout 2>errfile
     export err=$?; err_chk
     filelist="partmesh.txt PE*/fort.14 PE*/fort.18 PE*/fort.13 PE*/fort.24 PE*/elev_stat.151 PE*/vel_stat.151"
     tar cvzf ${RUN}_${ncpu}.tar.gz $filelist
     export err=$?; err_chk
     cpfs ${RUN}_${ncpu}.tar.gz $COMGES/.
  else
     mpiexec -n 1 -ppn 1 $EXECstofs/${RUN}_adcprep --np $ncpu --prep15 >> $pgmout 2>errfile
     export err=$?; err_chk
  fi
  
# --------------------------------------------------------------------------- #
# 6.  Execute stofs_2d_glo_padcirc for tide-only forecast

  export pgm="stofs_2d_glo_padcirc for tide-only forecast"
  . prep_step
  startmsg
  mpiexec -n $ncpu -ppn 127 --cpu-bind core $EXECstofs/${RUN}_padcirc >> $pgmout 2>adcirc.err
  export err=$?; err_chk
  if [[ -n $(grep 'ADCIRC stopping' adcirc.err) || -n $(grep 'ADCIRC Terminating' adcirc.err) ]]; then
     echo "FATAL ERROR: ADCIRC_RUN crashed and terminated"
     err_exit
  else
     echo "ADCIRC_RUN completed normally"
  fi

# --------------------------------------------------------------------------- #
# 7.  Send files to $COM

  if [ $SENDCOM = YES ]; then
     echo "Copying fort.61.nc to $COMOUT/${RUN}.${cycle}.points.htp.nc"
     cpfs fort.61.nc             $COMOUT/${RUN}.${cycle}.points.htp.nc
     echo "Copying fort.63.nc to $COMOUT/${RUN}.${cycle}.fields.htp.nc"
     cpfs fort.63.nc             $COMOUT/${RUN}.${cycle}.fields.htp.nc
  fi

  if [ $SENDDBN = YES ]; then
     $DBNROOT/bin/dbn_alert MODEL ${DBN_ALERT_TYPE} $job $COMOUT/${RUN}.${cycle}.points.htp.nc
     $DBNROOT/bin/dbn_alert MODEL ${DBN_ALERT_TYPE} $job $COMOUT/${RUN}.${cycle}.fields.htp.nc
  fi

  msg="Completing stofs_2d_glo_prep_forecast2 script"
  echo "$msg"
  postmsg "$jlogfile" "$msg"

# End of stofs_2d_glo_prep_forecast2.sh.ecf script -------------------------- #
