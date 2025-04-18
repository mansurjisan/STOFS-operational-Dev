#!/bin/bash

############################################################################
#  Name: stofs_3d_atl_create_surface_forcing_gfs.sh                        #
#  This script read the NCEP/GFS data to create the GFS based surface      #
#  forcing files, stofs_3d_atl.t12z.gfs.{air,prc,rad}.nc for the nowcast   #
#  and forecast simuations.                                                #
#                                                                          #
#  Remarks:                                                                #
#                                                        September, 2022   #
############################################################################


# ---------------------------> Begin ...
 set -x

echo 'stofs_3d_atl_create_surface_forcing_gfs.sh started at UTC' `date -u +%Y%m%d%H`


# ---------------------------> directory/file names


  dir_wk=${DATA_prep_gfs}/

  mkdir -p $dir_wk
  cd $dir_wk
  rm -fr $dir_wk/*

  mkdir -p ${COMOUTrerun}


  pgmout=pgmout_gfs.$$
  

# ---------------------------> Global Variables
  fn_nco_update_time_varName=${FIXstofs3d}/stofs_3d_atl_gfs_input_nco_update_var.nco
  fn_txt_sflux_inputs=${FIXstofs3d}/stofs_3d_atl_sflux_inputs.txt

  fn_gfs_rad_schism=sflux_rad_1.0001.nc
  #fn_gfs_rad_date_tag=${RUN}.gfs.rad.nfcast.${PDYHH_FCAST_BEGIN:0:8}.${cycle}.nc
  fn_gfs_rad_std=${RUN}.${cycle}.gfs.rad.nc

  fn_gfs_prc_schism=sflux_prc_1.0001.nc
  #fn_gfs_prc_date_tag=${RUN}.gfs.prc.nfcast.${PDYHH_FCAST_BEGIN:0:8}.${cycle}.nc
  fn_gfs_prc_std=${RUN}.${cycle}.gfs.prc.nc

  fn_gfs_air_schism=sflux_air_1.0001.nc
  #fn_gfs_air_date_tag=${RUN}.gfs.air.nfcast.${PDYHH_FCAST_BEGIN:0:8}.${cycle}.nc
  fn_gfs_air_std=${RUN}.${cycle}.gfs.air.nc
   

# ---------------> Region of interest
    LONMIN=-98.5
    LONMAX=-59.5
    LATMIN=8.0
    #    LATMAX=46.3 (updated to following: for sake of HRRR intpl)
    #   LATMAX=50.52
    LATMAX=50.751


    set +x
# ---------------> Dates
   yyyymmdd_today=${PDYHH_FCAST_BEGIN:0:8}
   yyyymmdd_prev=${PDYHH_NCAST_BEGIN:0:8}

   echo yyyymmdd_prev=$yyyymmdd_prev

    iyr=`echo ${yyyymmdd_prev} | cut -c1-4`
    imon=`echo ${yyyymmdd_prev} | cut  -c5-6`
    iday=`echo ${yyyymmdd_prev} | cut -c7-8`
    #ihr=12  # 00
    #hr_1st_file=0
 
# --------------------------> default: create list 
    list_fn_yest_t06z=`ls ${COMINgfs}/gfs.${yyyymmdd_prev}/06/atmos/gfs.t06z.pgrb2.0p25.f006`                >> $pgmout 2> errfile
    list_fn_yest_t12z=`ls ${COMINgfs}/gfs.${yyyymmdd_prev}/12/atmos/gfs.t12z.pgrb2.0p25.f00{1,2,3,4,5,6}`    >> $pgmout 2> errfile
    list_fn_yest_t18z=`ls ${COMINgfs}/gfs.${yyyymmdd_prev}/18/atmos/gfs.t18z.pgrb2.0p25.f00{1,2,3,4,5,6}`    >> $pgmout 2> errfile
    list_fn_today_t00z=`ls ${COMINgfs}/gfs.${yyyymmdd_today}/00/atmos/gfs.t00z.pgrb2.0p25.f00{1,2,3,4,5,6}`  >> $pgmout 2> errfile
    list_fn_today_t06z=`ls ${COMINgfs}/gfs.${yyyymmdd_today}/06/atmos/gfs.t06z.pgrb2.0p25.f00{1,2,3,4,5,6}`  >> $pgmout 2> errfile

     # form file list of today:t12z
     str_base_t12z=${COMINgfs}/gfs.${yyyymmdd_today}/12/atmos/gfs.t12z.pgrb2.0p25.f
     list_hr_oi_001_050=`seq -f "%03g" 1 1 50`

     list_fn_today_t12z=''
     for str_hhh in $list_hr_oi_001_050
     do
        fn_k=${str_base_t12z}${str_hhh}
        list_fn_today_t12z=$list_fn_today_t12z' '${fn_k}
     done

     # concatenate dir/file names
     LIST_fn_all_1="${list_fn_yest_t06z} "
     LIST_fn_all_1+="${list_fn_yest_t12z[@]} "
     LIST_fn_all_1+="${list_fn_yest_t18z[@]} "
     LIST_fn_all_1+="${list_fn_today_t00z[@]} "
     LIST_fn_all_1+="${list_fn_today_t06z[@]} "
     LIST_fn_all_1+="${list_fn_today_t12z[@]}"

     #for a in ${LIST_fn_all};  do echo $a; done


# ---------------------------> backup list
  list_fn_bk_1=`ls ${COMINgfs}/gfs.${yyyymmdd_prev}/06/atmos/gfs.t06z.pgrb2.0p25.f006`
  list_fn_bk_2=`ls ${COMINgfs}/gfs.${yyyymmdd_prev}/12/atmos/gfs.t12z.pgrb2.0p25.f00{1,2,3,4,5,6,7,8,9}`
  #list_fn_bk_3=`ls ${COMINgfs}/gfs.${yyyymmdd_prev}/12/atmos/gfs.t12z.pgrb2.0p25.f0{1,2,3,4,5,6,7,8}?`
  list_fn_bk_3=`ls ${COMINgfs}/gfs.${yyyymmdd_prev}/12/atmos/gfs.t12z.pgrb2.0p25.f0{1,2,3,4,5,6}?`
  list_fn_bk_4=`ls ${COMINgfs}/gfs.${yyyymmdd_prev}/12/atmos/gfs.t12z.pgrb2.0p25.f07{0,1,2,3,4}`

  LIST_fn_all_2="${list_fn_bk_1[@]} "
  LIST_fn_all_2+="${list_fn_bk_2[@]} "
  LIST_fn_all_2+="${list_fn_bk_3[@]} "
  LIST_fn_all_2+="${list_fn_bk_4[@]}"

  echo; echo "list_1"
  A=$LIST_fn_all_1; for a in ${A[@]}; do echo $a; done

  echo; echo "list_2"
  A=$LIST_fn_all_2; for a in ${A[@]}; do echo $a; done


# ------------------------> check size
list_route_no=(1 2)
for flag_route_no in ${list_route_no[@]}; do

 echo $flag_route_no
 if [[ $flag_route_no == 1 ]]; then
    list_wk=$LIST_fn_all_1
 else
    list_wk=$LIST_fn_all_2
 fi

 #echo ${list_wk[@]}


 FILESIZE=500000000
 LIST_fn_final=''
 for fn_gfs_k_sz in ${list_wk[@]}
 do
   echo "Processing:: " $fn_gfs_k_sz

   if [ -s $fn_gfs_k_sz ]; then
      filesize=`wc -c $fn_gfs_k_sz | awk '{print $1}' `

      if [ $filesize -ge $FILESIZE ];
      then
         LIST_fn_final+="${fn_gfs_k_sz} "
         echo "File size OK: $fn_gfs_k_sz : filesize $filesize GE $FILESIZE"
      else
         echo "WARNING: " $fn_gfs_k_sz ": filesize $filesize less than $FILESIZE"
         echo "WARNING: " $fn_gfs_k_sz ": filesize $filesize less than $FILESIZE"  >> $pgmout  # >> $jlogfile
      fi

   else
      echo "WARNING: "  $fn_gfs_k_sz " does not exist"
      echo "WARNING: "  $fn_gfs_k_sz " does not exist"  >> $jlogfile
   fi
 done


  if [[ $flag_route_no == 1 ]]; then
    LIST_fn_final_qa_sz_1=$LIST_fn_final
  else
    LIST_fn_final_qa_sz_2=$LIST_fn_final
  fi

done # for flag_route_no i


# ----------> combine if needed
# merge if missing
 N_list_target=74   # minCr=74: 73,74,75: 3.5, 3.541667, 3.583333

 A1=($LIST_fn_final_qa_sz_1)
 B2=($LIST_fn_final_qa_sz_2)

  N_list_1=${#A1[@]}; echo $N_list_1
  N_list_2=${#B2[@]}; echo $N_list_2;


if [[ ${N_list_1} -gt 1 ]]; then
 
  LIST_fn_final_qa_sz=(${A1[@]})

  if [[ ${N_list_1} -lt ${N_list_target} ]] && [[ ${N_list_2} -gt ${N_list_1} ]]; then
    echo "N_list_1 = $N_list_1"; echo "N_list_2 = $N_list_2"

    n_diff_1_2=$((${N_list_2}-${N_list_1}))

    # error   LIST_fn_final_qa_sz=${A1[@]} ${B2[@]:$N_list_1:$n_diff_1_2}
    LIST_fn_final_qa_sz=(${A1[@]} ${B2[@]:$N_list_1:$n_diff_1_2})

    echo "combined: LIST_fn_1 & 2: "
    for a in ${LIST_fn_final_qa_sz[@]}; do echo $a; done

  else
	  # LIST_fn_final_qa_sz=(${A1[@]})
    echo "List from LIST_fn_final_qa_sz_1"
 
  fi

elif  [[ ${N_list_2} -gt 1 ]]; then
  echo "List from LIST_fn_final_qa_sz_2"	
  LIST_fn_final_qa_sz=(${B2[@]})

else
	LIST_fn_final_qa_sz=()

fi

  set -x

# ---------------------> process data


 list_var_oi='TMP:2 m above|RH:2 m above|SPFH:2 m above|PRMSL|PRATE|UGRD:10 m above|VGRD:10 m above|ALBDO:surface|DSWRF:surface|USWRF:surface|DLWRF:surface|ULWRF:surface'

 rm -f *_voi*.
 rm -f *_sflux.nc

 ihr=12  # 00
 hr_1st_file=0

 # FYI
 let cnt="hr_1st_file-1"
 for fn_gfs_k in ${LIST_fn_final_qa_sz[@]}
 do
   let cnt=$cnt+1
   str_xxx_cnt=`seq -f "%03g" $cnt 1 $cnt`
   ln -sf $fn_gfs_k sorce_gfs_no_${str_xxx_cnt}
 done

N_dim_cr_min_cntList=48 
N_LIST_fn_final_qa_sz=${#LIST_fn_final_qa_sz[@]}

fn_merged_sflux=gfs_merge_v1.nc

echo; echo "N_LIST_fn_final_qa_sz = ${N_LIST_fn_final_qa_sz}"; echo 

if [[ ${N_LIST_fn_final_qa_sz} -gt ${N_dim_cr_min_cntList} ]]; then

  echo "N_LIST_fn_final_qa_sz = ${N_LIST_fn_final_qa_sz}"
  echo   

 let cnt="hr_1st_file-1"
 for fn_gfs_k in ${LIST_fn_final_qa_sz[@]}
 do

   let cnt=$cnt+1

   str_xxx_cnt=`seq -f "%03g" $cnt 1 $cnt`
   echo "Processing($str_xxx_cnt): " $fn_gfs_k

   # FYI:
   # ln -sf $fn_gfs_k sorce_gfs_no_${str_xxx_cnt} 


   fn_varOI=GFS_voi_${str_xxx_cnt}.grb2
      $WGRIB2  -s  $fn_gfs_k  | egrep "$list_var_oi" | $WGRIB2  -i  $fn_gfs_k  -grib  $fn_varOI  >> $pgmout 2> errfile
      export err=$?; #err_chk

   fn_roi=iGFS_voi_rio_${str_xxx_cnt}.grb2
      $WGRIB2  $fn_varOI  -small_grib ${LONMIN}:${LONMAX} ${LATMIN}:${LATMAX} $fn_roi   >> $pgmout 2> errfile
      export err=$?; #err_chk

   fn_0_rnVar=GFS_voi_rio_0rename_${str_xxx_cnt}.nc
      $WGRIB2  $fn_roi -netcdf $fn_0_rnVar  >> $pgmout 2> errfile
      export err=$?; #err_chk

   fn_out=GFS_sflux_no_${str_xxx_cnt}.nc


   str_time=`echo '"'days since $iyr-$imon-$iday 00:00:00'"'`
   let hr_cnt_since_hr00=${ihr}+${cnt}

     ncap2 -Oh -s "tin=${hr_cnt_since_hr00}"  -s "time@units=$str_time"  -s "time@base_date ={ $iyr, $imon, $iday, 0}" -S $fn_nco_update_time_varName -v ${fn_0_rnVar}  $fn_out   >> $pgmout 2> errfile
     export err=$?; #err_chk

 done

# merge GFS_sflux_no_xxx.nc
# fn_merged_sflux=gfs_merge_v1.nc
 rm -f ${fn_merged_sflux};

 echo fn_merged_sflux= $fn_merged_sflux

    rm -f $fn_merged_sflux
    find . -size 0  -exec rm -f {} \;

    list_GFS_sflux_no=`ls GFS_sflux_no_*.nc`
    #if [ ! -z "$list_GFS_sflux_no???" ]; then
    if [ ! -z "$list_GFS_sflux_no" ]; then
      ncrcat -O  GFS_sflux_no_*.nc  $fn_merged_sflux
    fi

fi   # if [[ ${N_LIST_fn_final_qa_sz} -gt ${N_dim_cr_min_cntList} ]]; then


# ---------------------------------> QC & archive
#list_var_ori=(elev2D.th TEM_3D.th SAL_3D.th uv3D.th TEM_nu SAL_nu)
#list_var_std=(elev2dth tem3dth sal3dth uv3dth temnu salnu) 

N_dim_cr_min=48
N_dim_cr_max=74
list_fn_sz_cr=(2000000)
list_end_time_step=(4.0)
list_offset_time=(1.0)

fn_ori=${fn_merged_sflux}
fn_std_1=${fn_gfs_rad_std}
fn_std_2=${fn_gfs_prc_std}
fn_std_3=${fn_gfs_air_std}

fn_std=${fn_std_1}
rm -f ${fn_std}

k=0
   if [[ -s ${fn_ori} ]]; then 
       sz_k=$((`wc -c ${fn_ori} | awk '{print $1}'`)) 
       
       if [[ ${sz_k} -gt ${list_fn_sz_cr[$k]} ]]; then 
            dim_k=`ncdump -h  ${fn_ori}  | grep "time = UNLIMITED" | awk -F'(' '{print $2}' | awk -F' ' '{print $1}'`; 
       else
            sz_k=$((0))
            dim_k=$((0))
       fi
   else
            sz_k=$((0))
            dim_k=$((0))
   fi
   echo "dim=${dim_k}, sz_k-bytes=${sz_k}, sz_cr=${list_fn_sz_cr[$k]}"  

 
   time_end_step=${list_end_time_step[$k]}
   time_offset=${list_offset_time[$k]}

   if [[ ${dim_k} -ge ${N_dim_cr_max} ]]; then
      ln -sf ${fn_ori} ${fn_std} 
      cpreq -pf ${fn_ori} ${COMOUTrerun}/${fn_std_1}
      cpreq -pf ${fn_ori} ${COMOUTrerun}/${fn_std_2}
      cpreq -pf ${fn_ori} ${COMOUTrerun}/${fn_std_3}
      echo "done: method - non-backup"

   elif [[ ${dim_k} -ge ${N_dim_cr_min} ]]; then
      ncap2 -s "time(-1)=${time_end_step}" ${fn_ori} -O ${fn_std}
      cpreq -pf ${fn_std} ${COMOUTrerun}/${fn_std_1}
      cpreq -pf ${fn_std} ${COMOUTrerun}/${fn_std_2}
      cpreq -pf ${fn_std} ${COMOUTrerun}/${fn_std_3}
      echo "done: method - backup 1"    

   else 
      if [[ -f  ${COMOUT_PREV}/rerun/${fn_std} ]]; then 

      echo "WARNING: GFS data not available; using backup data $COMOUT_PREV/rerun/${fn_std}"

      myr=${yyyymmdd_prev:0:4}
      mmon=${yyyymmdd_prev:4:2}
      mday=${yyyymmdd_prev:6:2}

      rm -f tmp*.nc

      fn_prev=prev_${fn_std}
      cpreq -pf ${COMOUT_PREV}/rerun/${fn_std} ${fn_prev}
      
      N_fn_prev=`ncdump -h ${fn_prev} | grep "time = UNLIMITED" | awk -F'(' '{print $2}' | awk -F' ' '{print $1}'`;

      if [[ ${N_fn_prev} -ge ${N_dim_cr_max} ]]; then
          fn_tmp1=tmp1_${fn_std}
          ncks -d time,12,,1 ${fn_prev} -O ${fn_tmp1}
          ncap2 -s "time=time-${time_offset}" ${fn_tmp1} tmp2_${fn_std}
          echo "time_offset= ${time_offset}, time_end_step= ${time_end_step}"

      fi    

      ncap2 -s "time(-1)=${time_end_step}" tmp2_${fn_std}  -O ${fn_std}
      ncatted -O -h -a units,time,o,c,"days since ${myr}-${mmon}-${mday} 00:00:00" ${fn_std}
      ncatted -O -h -a base_date,time,o,I,"${myr}, ${mmon}, ${mday}, 0" ${fn_std}
           
         cpreq -pf ${fn_std} ${COMOUTrerun}/${fn_std_1}
         cpreq -pf ${fn_std} ${COMOUTrerun}/${fn_std_2}
         cpreq -pf ${fn_std} ${COMOUTrerun}/${fn_std_3}
         echo "done: method - backup 2"   
      
      else
         msg="FATAL ERRORS: Missing critical GFS data in both current and previous day's ${COMINgfs}/gfs.${yyyymmdd_prev}/${cyc}/atmos/gfs.t${cyc}z.pgrb2.0p25.f000-074 as well as backup data ${COMOUT_PREV}/rerun/stofs_3d_atl.t12z.gfs.{rad,prc,air}"
         echo -e ${msg}
         exit 9

      fi

   fi # if [[ ${dim_k} -ge ${N_dim_cr_max} ]]



echo
echo "The script stofs_3d_atl_create_surface_forcing_gfs.sh completed at date/time: " `date`
echo 


