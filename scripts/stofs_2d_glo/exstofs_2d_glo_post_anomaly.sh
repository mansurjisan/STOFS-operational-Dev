#!/bin/bash
###############################################################################
#                                                                             #
# This script is the postprocessor for the STOFS that runs under the ADCIR    #
# model. It sets some shell script variables for export to child scripts      #
# and copies some generally used files to the work directory.                 #
# After this the actual preprocessing is performed by the following scripts:  #
#                                                                             #
# Remarks :                                                                   #
#                                                                             #
#                                                                 Sep, 2023   #
#                                                                             #
###############################################################################
# Start of exstofs_2d_glo_post_anomaly.sh.ecf script ------------------------ #
# 0.  Preparation
# 0.a Basic modes of operation

  seton='-xa'
  setoff='+xa'
  set $seton

  msg="Starting stofs_2d_glo_post_anomaly script"
  echo "$msg"
  postmsg "$jlogfile" "$msg"

  wndh=3
  nowh=6
  lsth=180

# --------------------------------------------------------------------------- #
# 1.  Set times
# 1.a Set all necessary times
#     YMDH     :  current time cycle in yyyymmddhh format
#     time_beg :  begin time of run (normally, -6 hour hindcast or the most recent restart time)
#     time_now :  current time
#     time_end :  ending time of run ($lsth hour forecast)

  export date=$PDY
  export YMDH=${PDY}${cyc}
  export nback=20

  time_beg=`$NDATE -6 $YMDH`
  time_now=$YMDH
  time_end=`$NDATE $lsth $YMDH`

  ymdh=$time_beg
  while [ $ymdh -lt $time_end  ]; do
	  for mins in $(seq -f "%02g" 6 6 54); do
		  echo ${ymdh}${mins} >> ymdhm.txt
	  done
	  ymdh=`$NDATE 1 $ymdh`
	  echo ${ymdh}00 >> ymdhm.txt
  done

# --------------------------------------------------------------------------- #
# 2.  Copy output files from $COMIN

  if [ -f $COMIN/${RUN}.${cycle}.points.cwl.nc ]; then
     cpreq $COMIN/${RUN}.${cycle}.points.cwl.nc cwl.fort.61.nc
     cpreq $COMIN/${RUN}.${cycle}.points.htp.nc htp.fort.61.nc
  fi

  ncdump -f f cwl.fort.61.nc > tmp.cdl
  export err=$?; err_chk
  time_number=$(grep -P "time = UNLIMITED" tmp.cdl | grep -o '[0-9]\+')
  station_number=$(grep -P "station =" tmp.cdl | grep -o '[0-9]\+')
  sed -i '31d' tmp.cdl; sed -i '1,/zeta =/!d' tmp.cdl

# --------------------------------------------------------------------------- #
# 3.  Create water level files for each station

  export pgm="stofs_2d_glo_anomaly"
  . prep_step
  startmsg

  extract_year=$(echo $time_now | cut -c1-4)
  extract_year2=$(echo $time_now | cut -c3-4)
  extract_month=$(echo $time_now | cut -c5-6)
  extract_day=$(echo $time_now | cut -c7-8)
  
  mpiexec -n 1 -ppn 1 $EXECstofs/${RUN}_anomaly "cwl.fort.61.nc" "htp.fort.61.nc" "${extract_day}" "${cyc}" "$FIXstofs/${RUN}_station.ctl" "ymdhm.txt" "mdl.t${cyc}z.txt" >> $pgmout 2>errfile
  export err=$?; err_chk

  if [ ! -d $DATA/model ]; then
	mkdir $DATA/model
  fi
  mv *.cwl $DATA/model/.
  cpreq mdl.${cycle}.txt $DATA/model/${RUN}_1hcwl.txt

# --------------------------------------------------------------------------- #
# 4.  Compute anomaly for each station

  if [ ! -d $DATA/data ]; then
	mkdir $DATA/data $DATA/database $DATA/log $DATA/msl $DATA/msl/plots $DATA/msl/maps
  fi
  cpreq $FIXstofs/${RUN}_station.ctl $DATA/data/.
  ln -s $FIXstofs/${RUN}_cron.bnt $DATA/data/cron.bnt
  ln -s $FIXstofs/${RUN}_ft03.dta $DATA/data/ft03.dta
  ln -s $FIXstofs/${RUN}_ft07.dta $DATA/data/ft07.dta
  ln -s $FIXstofs/*.gif $DATA/data/.

  if [ -f $COMIN/database.tar.gz ]; then
	cpreq $COMIN/database.tar.gz $DATA/.
      	tar xvzf database.tar.gz
        export err=$?; err_chk
  else
	mkdir $DATA/database
  fi
 
  sed -i '1d' $DATA/data/${RUN}_station.ctl
  while read obs; do
	  nosid=$(echo $obs | awk '{print $2}')
	  msl2mllw=$(echo $obs | awk '{print $3}') 
	  echo $nosid $msl2mllw
	  if [ -f $DCOMIN/${nosid}.xml ]; then
		  awk -v m=$msl2mllw '{print $4,$5",",($6+m)*3.2808","}' $DCOMIN/${nosid}.xml > $DATA/database/${nosid}.csv
          else
		  echo "NULL" > $DATA/database/${nosid}.csv
	  fi
  done < $DATA/data/${RUN}_station.ctl
  
  mpiexec -n 1 -ppn 1 $USHstofs/etweb_database.tcl >> $pgmout 2>errfile
  export err=$?; err_chk

  extract_now=$(echo ${extract_month}/${extract_day}/${extract_year} ${cyc}:00:00)
  mpiexec -n 1 -ppn 1 $USHstofs/etweb_extract.tcl --all --date ${extract_now} >> $pgmout 2>errfile
  export err=$?; err_chk

  rm $DATA/database/*.csv
  if [ -d $DATA/database ]; then
       tar cvzf database.tar.gz database/
       export err=$?; err_chk
       cpfs database.tar.gz $COMOUT/.
  fi
 
# --------------------------------------------------------------------------- #
# 5.  Combine water level and anomaly for each station

  if [ ! -d $DATA/tmp ]; then
	mkdir $DATA/tmp 
  fi

  extract_date=$(echo ${extract_year2}${extract_month}${extract_day}${cyc})
  while read line; do
        abbr=$(echo $line | awk '{print $1}')
        if [ -f $DATA/msl/plots/${abbr}_${extract_date}.anom ]; then
        	anom=$(tail -n 1 $DATA/msl/plots/${abbr}_${extract_date}.anom | awk '{print $2}')
        	if [ ${anom} != 0.00 ]; then
        		while read line1; do
        			awk -f $USHstofs/inter.awk $DATA/msl/plots/${abbr}_${extract_date}.anom > tmp.txt
        			export err=$?; err_chk
        		done < $DATA/msl/plots/${abbr}_${extract_date}.anom
                        
        		grep -e ${time_beg} -A 4000 tmp.txt > $DATA/tmp/anom.tmp
        		sed -i '1d' $DATA/tmp/anom.tmp
        		cons_anom=$(tail -n 1 tmp.txt | awk '{print $2}')
        		cut -f2 -d' ' $DATA/tmp/anom.tmp > $DATA/tmp/${abbr}.6manom
        		line_number=$(wc -l < $DATA/tmp/${abbr}.6manom)
        		line_filt=$((${time_number}-${line_number}))
              		for i in $(seq 1 ${line_filt}); do
        			echo $cons_anom >> $DATA/tmp/${abbr}.6manom
        		done
 
        		paste $DATA/model/${abbr}.cwl $DATA/tmp/${abbr}.6manom | awk '{print $2 + $3}' > $DATA/tmp/${abbr}.acwl
        		export err=$?; err_chk
        		awk -f $USHstofs/transpose.awk $DATA/tmp/${abbr}.acwl >> ${station_number}_stations.csv 
        		export err=$?; err_chk
        		rm tmp.txt $DATA/tmp/anom.tmp
        	else
        		awk '{print $2}' $DATA/model/${abbr}.cwl > $DATA/tmp/${abbr}.acwl
        		export err=$?; err_chk
        		awk -f $USHstofs/transpose.awk $DATA/tmp/${abbr}.acwl >> ${station_number}_stations.csv 
        		export err=$?; err_chk
        	fi
        else
        	awk '{print $2}' $DATA/model/${abbr}.cwl > $DATA/tmp/${abbr}.acwl
        	export err=$?; err_chk
        	awk -f $USHstofs/transpose.awk $DATA/tmp/${abbr}.acwl >> ${station_number}_stations.csv 
        	export err=$?; err_chk
        fi
  done < $DATA/data/${RUN}_station.ctl
  
# --------------------------------------------------------------------------- #
# 6.  Generate NetCDF file for combined water level with anomaly

  for i in $(seq 1 ${time_number}) ; do
        awk -v i=${i} '{print $i","}' ${station_number}_stations.csv >> tmp2.txt
        export err=$?; err_chk
  done
  sed -i '$ s/,/;/g' tmp2.txt; echo '}' >> tmp2.txt
  cat tmp.cdl tmp2.txt > cwl.fort.61.cdl
  ncgen -o acwl.fort.61.nc cwl.fort.61.cdl
  export err=$?; err_chk
  rm *.cdl *.csv tmp.txt tmp2.txt

# --------------------------------------------------------------------------- #
# 7.  Send files to $COM

  if [ $SENDCOM = YES ]; then
     echo "Copying fort.61.nc to $COMOUT/${RUN}.${cycle}.points.cwl.nc"
     cpfs acwl.fort.61.nc        $COMOUT/${RUN}.${cycle}.points.cwl.nc
     echo "Copying fort.63.nc to $COMOUT/${RUN}.${cycle}.fields.cwl.nc"
     cpfs cwl.fort.63.nc         $COMOUT/${RUN}.${cycle}.fields.cwl.nc
     if [ $SENDDBN = YES ]; then
        $DBNROOT/bin/dbn_alert MODEL ${DBN_ALERT_TYPE} $job $COMOUT/${RUN}.${cycle}.points.cwl.nc
        $DBNROOT/bin/dbn_alert MODEL ${DBN_ALERT_TYPE} $job $COMOUT/${RUN}.${cycle}.fields.cwl.nc
     fi
  fi

  msg="Completing stofs_2d_glo_post_anomaly script"
  echo "$msg"
  postmsg "$jlogfile" "$msg"

# End of exstofs_2d_glo_post_anomly.sh.ecf script -------------------------- #
