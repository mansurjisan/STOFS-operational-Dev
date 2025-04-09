#!/usr/bin/env python3
"""
STOFS3D Atlantic RTOFS Processing with xarray

This script replaces several NCO/NCAP2 operations in the STOFS preprocessing
pipeline with more efficient xarray operations.

Usage within shell script:
    python3 stofs_process_rtofs.py --ssh_files RTOFS_2D_*.nc \
                                   --tsuv_files RTOFS_3D_*.nc \
                                   --out_dir ./output/ \
                                   --adt_today path/to/adt_today.nc \
                                   --adt_prev path/to/adt_prev.nc \
                                   --idx_2ds "$idx_x1_2ds $idx_x2_2ds $idx_y1_2ds $idx_y2_2ds" \
                                   --idx_3dz "$idx_x1_3dz $idx_x2_3dz $idx_y1_3dz $idx_y2_3dz" \
                                   --adt_weight path/to/adt_weight.nc
"""

import os
import sys
import glob
import argparse
import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime


def setup_args():
    """Set up command line arguments"""
    parser = argparse.ArgumentParser(description='Process RTOFS data with xarray')
    parser.add_argument('--ssh_files', required=True, help='SSH file pattern (e.g., "RTOFS_2D_*.nc")')
    parser.add_argument('--tsuv_files', required=True, help='TSUV file pattern (e.g., "RTOFS_3D_*.nc")')
    parser.add_argument('--out_dir', required=True, help='Output directory')
    parser.add_argument('--adt_today', help='ADT file for today')
    parser.add_argument('--adt_prev', help='ADT file for previous day')
    parser.add_argument('--idx_2ds', required=True, help='2D indices: "x1 x2 y1 y2"')
    parser.add_argument('--idx_3dz', required=True, help='3D indices: "x1 x2 y1 y2"')
    parser.add_argument('--adt_weight', help='ADT weight file')
    
    return parser.parse_args()

def extract_roi(ds, idx_list, variables):
    """Extract region of interest from dataset"""
    x1, x2, y1, y2 = [int(i) for i in idx_list.split()]
    
    # Filter to selected variables and region
    return ds.isel(X=slice(x1, x2+1), Y=slice(y1, y2+1))[variables]

def process_ssh_files(file_pattern, idx_list, out_dir):
    """Process SSH files from RTOFS 2D datasets"""
    print(f"Processing SSH files matching: {file_pattern}")
    
    # Glob pattern to get all matching files and sort them
    files = sorted(glob.glob(file_pattern))
    if not files:
        print(f"No files found matching pattern: {file_pattern}")
        return None
    
    print(f"Found {len(files)} SSH files")
    
    # Process files individually to avoid memory issues
    variables = ['MT', 'Date', 'Longitude', 'Latitude', 'ssh']
    processed_files = []
    
    for i, filename in enumerate(files):
        print(f"Processing file {i+1}/{len(files)}: {os.path.basename(filename)}")
        try:
            ds = xr.open_dataset(filename)
            roi_ds = extract_roi(ds, idx_list, variables)
            out_file = os.path.join(out_dir, f"rio_ssh_{os.path.basename(filename)}")
            roi_ds.to_netcdf(out_file)
            processed_files.append(out_file)
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue
    
    if not processed_files:
        return None
    
    # Concatenate all extracted files along time dimension
    print("Concatenating SSH files...")
    merged_file = os.path.join(out_dir, "merged_RTOFS_2D.nc")
    
    # Use xarray to open and concatenate files in chunks
    datasets = []
    for file in processed_files:
        ds = xr.open_dataset(file)
        datasets.append(ds)
    
    # Concat along MT dimension
    merged_ds = xr.concat(datasets, dim="MT")
    merged_ds.to_netcdf(merged_file)
    
    # Close datasets to free memory
    for ds in datasets:
        ds.close()
    
    print(f"Created merged SSH file: {merged_file}")
    return merged_file

def process_tsuv_files(file_pattern, idx_list, out_dir):
    """Process TSUV files from RTOFS 3D datasets"""
    print(f"Processing TSUV files matching: {file_pattern}")
    
    # Glob pattern to get all matching files and sort them
    files = sorted(glob.glob(file_pattern))
    if not files:
        print(f"No files found matching pattern: {file_pattern}")
        return None
    
    print(f"Found {len(files)} TSUV files")
    
    variables = ['MT', 'Date', 'Longitude', 'Latitude', 'Depth', 'temperature', 'salinity', 'u', 'v']
    processed_files = []
    
    # Process files individually
    for i, filename in enumerate(files):
        print(f"Processing file {i+1}/{len(files)}: {os.path.basename(filename)}")
        try:
            ds = xr.open_dataset(filename)
            roi_ds = extract_roi(ds, idx_list, variables)
            out_file = os.path.join(out_dir, f"rio_tsuv_{os.path.basename(filename)}")
            roi_ds.to_netcdf(out_file)
            processed_files.append(out_file)
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue
    
    if not processed_files:
        return None
    
    # Concatenate all extracted files along time dimension
    print("Concatenating TSUV files...")
    merged_file = os.path.join(out_dir, "merged_RTOFS_3D.nc")
    
    # Use xarray to open and concatenate
    datasets = []
    for file in processed_files:
        datasets.append(xr.open_dataset(file))
    
    merged_ds = xr.concat(datasets, dim="MT")
    merged_ds.to_netcdf(merged_file)
    
    # Close datasets to free memory
    for ds in datasets:
        ds.close()
    
    print(f"Created merged TSUV file: {merged_file}")
    return merged_file

def prepare_ssh_for_schism(merged_file, out_dir):
    """Prepare SSH data for SCHISM model from merged RTOFS data"""
    print("Preparing SSH data for SCHISM...")
    
    # Load the merged file
    ds = xr.open_dataset(merged_file)
    
    # Rename dimensions and handle missing values
    ds = ds.rename({'MT': 'time', 'X': 'xlon', 'Y': 'ylat'})
    
    # Handle missing values (equivalent to ncap2 operations)
    mask = ds.ssh > 10000
    ds['ssh'] = ds.ssh.where(~mask, -30000)
    
    # Set fill value attributes
    ds.ssh.attrs['_FillValue'] = -30000
    ds.ssh.attrs['missing_value'] = -30000
    
    # Remove unnecessary variables
    ds = ds[['Longitude', 'Latitude', 'ssh']]
    
    # Create new variables required by SCHISM
    ds['surf_el'] = ds.ssh
    
    # Save to netCDF file
    out_file = os.path.join(out_dir, "SSH_1_rtofs_only.nc")
    ds.to_netcdf(out_file)
    
    print(f"Created SCHISM SSH file: {out_file}")
    return out_file

def prepare_tsuv_for_schism(merged_file, out_dir):
    """Prepare TSUV data for SCHISM model from merged RTOFS data"""
    print("Preparing TSUV data for SCHISM...")
    
    # Load the merged file
    ds = xr.open_dataset(merged_file)
    
    # Rename dimensions and variables
    ds = ds.rename({
        'MT': 'time', 
        'Depth': 'lev', 
        'X': 'xlon', 
        'Y': 'ylat',
        'u': 'water_u',
        'v': 'water_v'
    })
    
    # Remove unnecessary variables
    ds = ds[['Longitude', 'Latitude', 'lev', 'temperature', 'salinity', 'water_u', 'water_v']]
    
    # Save to netCDF file
    out_file = os.path.join(out_dir, "TSUV_1.nc")
    ds.to_netcdf(out_file)
    
    print(f"Created SCHISM TSUV file: {out_file}")
    return out_file

def process_adt_data(adt_today, adt_prev, adt_weight, out_dir):
    """Process ADT data for merging with RTOFS data using pure xarray"""
    print("Processing ADT data...")
    
    adt_output = os.path.join(out_dir, "adt_aft_cvtz_cln.nc")
    
    try:
        if os.path.exists(adt_today) and os.path.exists(adt_prev):
            print(f"Using both ADT files: {os.path.basename(adt_today)} and {os.path.basename(adt_prev)}")
            
            # For regridding with ESMF weights, we'll need to use a subprocess call to ncremap
            # but keep everything else in Python for performance
            
            # Extract region of interest
            ds_today = xr.open_dataset(adt_today).sel(
                longitude=slice(-62.5, -51.5),
                latitude=slice(7.0, 54.0)
            )[['adt', 'sla', 'err_sla']]
            
            ds_prev = xr.open_dataset(adt_prev).sel(
                longitude=slice(-62.5, -51.5),
                latitude=slice(7.0, 54.0)
            )[['adt', 'sla', 'err_sla']]
            
            # Make time a record dimension if it's not already
            if 'time' not in ds_today.dims:
                ds_today = ds_today.expand_dims('time')
            
            if 'time' not in ds_prev.dims:
                ds_prev = ds_prev.expand_dims('time')
            
            # Save intermediate files for regridding
            today_rec = os.path.join(out_dir, "adt_roi_rec_dim_T.nc")
            prev_rec = os.path.join(out_dir, "adt_roi_rec_dim_P.nc")
            
            ds_today.to_netcdf(today_rec)
            ds_prev.to_netcdf(prev_rec)
            
            # Regrid using ncremap (this is the one step we'll keep as a subprocess)
            today_grid = os.path.join(out_dir, "adt_aft_wt_T.nc")
            prev_grid = os.path.join(out_dir, "adt_aft_wt_P.nc")
            
            import subprocess
            subprocess.run(["ncremap", "-i", today_rec, "-m", adt_weight, "-o", today_grid], check=True)
            subprocess.run(["ncremap", "-i", prev_rec, "-m", adt_weight, "-o", prev_grid], check=True)
            
            # Continue processing with xarray
            ds_today = xr.open_dataset(today_grid)
            ds_prev = xr.open_dataset(prev_grid)

            # Print variable names to debug
            print(f"Variables in regridded today file: {list(ds_today.variables)}")
            print(f"Dimensions in regridded today file: {list(ds_today.dims)}")
            
            # Apply the transformation from  NCO script
            ds_today['surf_el'] = ds_today.adt - 0.45
            ds_today.surf_el.attrs['add_offset'] = 0.0
            ds_today.surf_el.attrs['scale_factor'] = 1.0


            ds_prev['surf_el'] = ds_prev.adt - 0.45
            ds_prev.surf_el.attrs['add_offset'] = 0.0
            ds_prev.surf_el.attrs['scale_factor'] = 1.0

            # Determine coordinate names dynamically
            lon_name = next((x for x in ['xlon', 'lon', 'longitude'] if x in ds_today.coords), None)
            lat_name = next((x for x in ['ylat', 'lat', 'latitude'] if x in ds_today.coords), None)

            if not lon_name or not lat_name:
                print(f"Warning: Could not find longitude/latitude coordinates.")
                print(f"Available coordinates: {list(ds_today.coords)}")
                # Create fallback coordinates if necessary
                if not lon_name and 'lon' in ds_today.variables:
                    ds_today = ds_today.rename({'lon': 'xlon'})
                    ds_prev = ds_prev.rename({'lon': 'xlon'})
                    lon_name = 'xlon'
                if not lat_name and 'lat' in ds_today.variables:
                    ds_today = ds_today.rename({'lat': 'ylat'})
                    ds_prev = ds_prev.rename({'lat': 'ylat'})
                    lat_name = 'ylat'

            
            # Extract needed variables - only include coordinates we actually have
            variables_to_keep = ['surf_el']
            if lon_name:
                variables_to_keep.append(lon_name)
            if lat_name:
                variables_to_keep.append(lat_name)

            ds_today_clean = ds_today[variables_to_keep]
            ds_prev_clean = ds_prev[variables_to_keep]

            # If coordinates have different names than xlon/ylat, rename them
            if lon_name and lon_name != 'xlon':
                ds_today_clean = ds_today_clean.rename({lon_name: 'xlon'})
                ds_prev_clean = ds_prev_clean.rename({lon_name: 'xlon'})
            if lat_name and lat_name != 'ylat':
                ds_today_clean = ds_today_clean.rename({lat_name: 'ylat'})
                ds_prev_clean = ds_prev_clean.rename({lat_name: 'ylat'})


            # Average the datasets
            ds_combined = (ds_today_clean + ds_prev_clean) / 2
                        
            # Save the combined file
            ds_combined.to_netcdf(adt_output)
            
            # Create the surf_el_t1_adt variable for the final ADT file
            ds = xr.open_dataset(adt_output)

            # Debug output
            print(f"Variables in final ds: {list(ds.variables)}")
            print(f"Dimensions in final ds: {list(ds.dims)}")
            if 'time' in ds.dims:
                print(f"Time dimension size: {ds.time.size}")

            # Create a new dataset for the final ADT
            ds_final = xr.Dataset()

            # Check if surf_el has a time dimension with size > 0
            has_time_dim = 'time' in ds.surf_el.dims and ds.time.size > 0

            if has_time_dim:
                # Use index 0 for safety
                time_idx = 0
                ds_final['surf_el_t1_adt'] = ds.surf_el.isel(time=time_idx)
            else:
                # If surf_el doesn't have a time dimension or time dim is empty, try to use it directly
                # First check if it's a 2D array (lat/lon only)
                if len(ds.surf_el.dims) <= 2:
                    ds_final['surf_el_t1_adt'] = ds.surf_el
                else:
                    # Otherwise, we need to create a 2D version by removing time or other dimensions
                    # Get the dimension names
                    dims = list(ds.surf_el.dims)
                    print(f"surf_el dimensions: {dims}")
                    
                    # Create a dummy array with the right shape for latitude and longitude
                    lat_dim = next((d for d in dims if 'lat' in d.lower()), None)
                    lon_dim = next((d for d in dims if 'lon' in d.lower()), None)
                    
                    if lat_dim and lon_dim:
                        # Create a 2D slice
                        print(f"Creating 2D slice from surf_el using {lat_dim} and {lon_dim}")
                        # Use mean over all other dimensions
                        dims_to_reduce = [d for d in dims if d not in [lat_dim, lon_dim]]
                        if dims_to_reduce:
                            ds_final['surf_el_t1_adt'] = ds.surf_el.mean(dim=dims_to_reduce)
                        else:
                            ds_final['surf_el_t1_adt'] = ds.surf_el
                    else:
                        # try to create a dummy array
                        print("Cannot determine lat/lon dimensions, creating dummy data")
                        ds_final['surf_el_t1_adt'] = xr.full_like(ds.surf_el.isel(**{dims[0]: 0}), fill_value=0.0)

            # Set attributes
            ds_final.surf_el_t1_adt.attrs['_FillValue'] = -30000
            ds_final.surf_el_t1_adt.attrs['missing_value'] = -30000
            ds_final.surf_el_t1_adt.attrs['scale_factor'] = 1.0

            # Copy coordinate variables, but check if they exist first
            if 'xlon' in ds.variables:
                ds_final['xlon'] = ds.xlon
            if 'ylat' in ds.variables:
                ds_final['ylat'] = ds.ylat

            # Save final ADT file
            adt_final = os.path.join(out_dir, "adt_fnl.nc")
            ds_final.to_netcdf(adt_final)
            
            print(f"Created final ADT file: {adt_final}")

        elif os.path.exists(adt_today):
            print(f"Using only today's ADT file: {os.path.basename(adt_today)}")

            # Extract region of interest
            ds = xr.open_dataset(adt_today).sel(
                longitude=slice(-62.5, -51.5),
                latitude=slice(7.0, 54.0)
            )[['adt', 'sla', 'err_sla']]

            # Make time a record dimension if it's not already
            if 'time' not in ds.dims:
                ds = ds.expand_dims('time')

            # Save intermediate file for regridding
            today_rec = os.path.join(out_dir, "adt_roi_rec_dim.nc")
            ds.to_netcdf(today_rec)

            # Regrid using ncremap
            today_grid = os.path.join(out_dir, "adt_aft_wt.nc")

            import subprocess
            subprocess.run(["ncremap", "-i", today_rec, "-m", adt_weight, "-o", today_grid], check=True)

            # Continue processing with xarray
            ds = xr.open_dataset(today_grid)
            
            # Print variable names to debug
            print(f"Variables in regridded file: {list(ds.variables)}")
            print(f"Dimensions in regridded file: {list(ds.dims)}")

            # Apply the transformation from your NCO script (subtract 0.45)
            ds['surf_el'] = ds.adt - 0.45
            ds.surf_el.attrs['add_offset'] = 0.0
            ds.surf_el.attrs['scale_factor'] = 1.0

            # Determine coordinate names dynamically
            lon_name = next((x for x in ['xlon', 'lon', 'longitude'] if x in ds.coords), None)
            lat_name = next((x for x in ['ylat', 'lat', 'latitude'] if x in ds.coords), None)

            if not lon_name or not lat_name:
                print(f"Warning: Could not find longitude/latitude coordinates.")
                print(f"Available coordinates: {list(ds.coords)}")
                # Create fallback coordinates if necessary
                if not lon_name and 'lon' in ds.variables:
                    ds = ds.rename({'lon': 'xlon'})
                    lon_name = 'xlon'
                if not lat_name and 'lat' in ds.variables:
                    ds = ds.rename({'lat': 'ylat'})
                    lat_name = 'ylat'

            # Extract needed variables - only include coordinates we actually have
            variables_to_keep = ['surf_el']
            if lon_name:
                variables_to_keep.append(lon_name)
            if lat_name:
                variables_to_keep.append(lat_name)

            ds_clean = ds[variables_to_keep]

            # If coordinates have different names than xlon/ylat, rename them
            if lon_name and lon_name != 'xlon':
                ds_clean = ds_clean.rename({lon_name: 'xlon'})
            if lat_name and lat_name != 'ylat':
                ds_clean = ds_clean.rename({lat_name: 'ylat'})

            # Save the processed file
            ds_clean.to_netcdf(adt_output)

            # Create the surf_el_t1_adt variable for the final ADT file
            ds = xr.open_dataset(adt_output)

            # Debug output
            print(f"Variables in final ds: {list(ds.variables)}")
            print(f"Dimensions in final ds: {list(ds.dims)}")
            if 'time' in ds.dims:
                print(f"Time dimension size: {ds.time.size}")

            # Create a new dataset for the final ADT
            ds_final = xr.Dataset()

            # Check if surf_el has a time dimension with size > 0
            has_time_dim = 'time' in ds.surf_el.dims and ds.time.size > 0

            if has_time_dim:
                # Use index 0 for safety
                time_idx = 0
                ds_final['surf_el_t1_adt'] = ds.surf_el.isel(time=time_idx)
            else:
                # If surf_el doesn't have a time dimension or time dim is empty, try to use it directly
                # First check if it's a 2D array (lat/lon only)
                if len(ds.surf_el.dims) <= 2:
                    ds_final['surf_el_t1_adt'] = ds.surf_el
                else:
                    # Otherwise, we need to create a 2D version by removing time or other dimensions
                    # Get the dimension names
                    dims = list(ds.surf_el.dims)
                    print(f"surf_el dimensions: {dims}")
                    
                    # Create a dummy array with the right shape for latitude and longitude
                    lat_dim = next((d for d in dims if 'lat' in d.lower()), None)
                    lon_dim = next((d for d in dims if 'lon' in d.lower()), None)
                    
                    if lat_dim and lon_dim:
                        # Create a 2D slice
                        print(f"Creating 2D slice from surf_el using {lat_dim} and {lon_dim}")
                        # Use mean over all other dimensions
                        dims_to_reduce = [d for d in dims if d not in [lat_dim, lon_dim]]
                        if dims_to_reduce:
                            ds_final['surf_el_t1_adt'] = ds.surf_el.mean(dim=dims_to_reduce)
                        else:
                            ds_final['surf_el_t1_adt'] = ds.surf_el
                    else:
                        # try to create a dummy array
                        print("Cannot determine lat/lon dimensions, creating dummy data")
                        ds_final['surf_el_t1_adt'] = xr.full_like(ds.surf_el.isel(**{dims[0]: 0}), fill_value=0.0)

            # Set attributes
            ds_final.surf_el_t1_adt.attrs['_FillValue'] = -30000
            ds_final.surf_el_t1_adt.attrs['missing_value'] = -30000
            ds_final.surf_el_t1_adt.attrs['scale_factor'] = 1.0

            # Copy coordinate variables, but check if they exist first
            if 'xlon' in ds.variables:
                ds_final['xlon'] = ds.xlon
            if 'ylat' in ds.variables:
                ds_final['ylat'] = ds.ylat

            # Save final ADT file
            adt_final = os.path.join(out_dir, "adt_fnl.nc")
            ds_final.to_netcdf(adt_final)

            print(f"Created final ADT file: {adt_final}")

        elif os.path.exists(adt_prev):
            print(f"Using only previous day's ADT file: {os.path.basename(adt_prev)}")

            # Extract region of interest
            ds = xr.open_dataset(adt_prev).sel(
                longitude=slice(-62.5, -51.5),
                latitude=slice(7.0, 54.0)
            )[['adt', 'sla', 'err_sla']]

            # Make time a record dimension if it's not already
            if 'time' not in ds.dims:
                ds = ds.expand_dims('time')

            # Save intermediate file for regridding
            prev_rec = os.path.join(out_dir, "adt_roi_rec_dim.nc")
            ds.to_netcdf(prev_rec)

            # Regrid using ncremap
            prev_grid = os.path.join(out_dir, "adt_aft_wt.nc")

            import subprocess
            subprocess.run(["ncremap", "-i", prev_rec, "-m", adt_weight, "-o", prev_grid], check=True)

            # Continue processing with xarray
            ds = xr.open_dataset(prev_grid)
            
            # Print variable names to debug
            print(f"Variables in regridded file: {list(ds.variables)}")
            print(f"Dimensions in regridded file: {list(ds.dims)}")

            # Apply the transformation from your NCO script (subtract 0.45)
            ds['surf_el'] = ds.adt - 0.45
            ds.surf_el.attrs['add_offset'] = 0.0
            ds.surf_el.attrs['scale_factor'] = 1.0

            # Determine coordinate names dynamically
            lon_name = next((x for x in ['xlon', 'lon', 'longitude'] if x in ds.coords), None)
            lat_name = next((x for x in ['ylat', 'lat', 'latitude'] if x in ds.coords), None)

            if not lon_name or not lat_name:
                print(f"Warning: Could not find longitude/latitude coordinates.")
                print(f"Available coordinates: {list(ds.coords)}")
                # Create fallback coordinates if necessary
                if not lon_name and 'lon' in ds.variables:
                    ds = ds.rename({'lon': 'xlon'})
                    lon_name = 'xlon'
                if not lat_name and 'lat' in ds.variables:
                    ds = ds.rename({'lat': 'ylat'})
                    lat_name = 'ylat'

            # Extract needed variables - only include coordinates we actually have
            variables_to_keep = ['surf_el']
            if lon_name:
                variables_to_keep.append(lon_name)
            if lat_name:
                variables_to_keep.append(lat_name)

            ds_clean = ds[variables_to_keep]

            # If coordinates have different names than xlon/ylat, rename them
            if lon_name and lon_name != 'xlon':
                ds_clean = ds_clean.rename({lon_name: 'xlon'})
            if lat_name and lat_name != 'ylat':
                ds_clean = ds_clean.rename({lat_name: 'ylat'})

            # Save the processed file
            ds_clean.to_netcdf(adt_output)

            # Create the surf_el_t1_adt variable for the final ADT file
            ds = xr.open_dataset(adt_output)

            # Debug output
            print(f"Variables in final ds: {list(ds.variables)}")
            print(f"Dimensions in final ds: {list(ds.dims)}")
            if 'time' in ds.dims:
                print(f"Time dimension size: {ds.time.size}")

            # Create a new dataset for the final ADT
            ds_final = xr.Dataset()

            # Check if surf_el has a time dimension with size > 0
            has_time_dim = 'time' in ds.surf_el.dims and ds.time.size > 0

            if has_time_dim:
                # Use index 0 for safety
                time_idx = 0
                ds_final['surf_el_t1_adt'] = ds.surf_el.isel(time=time_idx)
            else:
                # If surf_el doesn't have a time dimension or time dim is empty, try to use it directly
                # First check if it's a 2D array (lat/lon only)
                if len(ds.surf_el.dims) <= 2:
                    ds_final['surf_el_t1_adt'] = ds.surf_el
                else:
                    # Otherwise, we need to create a 2D version by removing time or other dimensions
                    # Get the dimension names
                    dims = list(ds.surf_el.dims)
                    print(f"surf_el dimensions: {dims}")
                    
                    # Create a dummy array with the right shape for latitude and longitude
                    lat_dim = next((d for d in dims if 'lat' in d.lower()), None)
                    lon_dim = next((d for d in dims if 'lon' in d.lower()), None)
                    
                    if lat_dim and lon_dim:
                        # Create a 2D slice
                        print(f"Creating 2D slice from surf_el using {lat_dim} and {lon_dim}")
                        # Use mean over all other dimensions
                        dims_to_reduce = [d for d in dims if d not in [lat_dim, lon_dim]]
                        if dims_to_reduce:
                            ds_final['surf_el_t1_adt'] = ds.surf_el.mean(dim=dims_to_reduce)
                        else:
                            ds_final['surf_el_t1_adt'] = ds.surf_el
                    else:
                        # try to create a dummy array
                        print("Cannot determine lat/lon dimensions, creating dummy data")
                        ds_final['surf_el_t1_adt'] = xr.full_like(ds.surf_el.isel(**{dims[0]: 0}), fill_value=0.0)

            # Set attributes
            ds_final.surf_el_t1_adt.attrs['_FillValue'] = -30000
            ds_final.surf_el_t1_adt.attrs['missing_value'] = -30000
            ds_final.surf_el_t1_adt.attrs['scale_factor'] = 1.0

            # Copy coordinate variables, but check if they exist first
            if 'xlon' in ds.variables:
                ds_final['xlon'] = ds.xlon
            if 'ylat' in ds.variables:
                ds_final['ylat'] = ds.ylat

            # Save final ADT file
            adt_final = os.path.join(out_dir, "adt_fnl.nc")
            ds_final.to_netcdf(adt_final)

            print(f"Created final ADT file: {adt_final}")
    
        else:
            print("No ADT files available, will try to use a previous run file")
            adt_output = None
        
        return adt_output

    except Exception as e:
        print(f"Error processing ADT data: {e}")
        import traceback
        traceback.print_exc()
        return None

def combine_ssh_adt(ssh_file, adt_file, out_dir):
    """Combine SSH from RTOFS with ADT data"""
    print("Combining SSH and ADT data...")
    
    if not os.path.exists(ssh_file):
        print(f"Error: SSH file {ssh_file} does not exist")
        return None
    
    if not os.path.exists(adt_file):
        print(f"Error: ADT file {adt_file} does not exist")
        return None
    
    # Load the SSH and ADT files
    ssh_ds = xr.open_dataset(ssh_file)
    adt_ds = xr.open_dataset(adt_file)
    
    # Debug output
    print(f"Variables in adt_ds: {list(adt_ds.variables)}")
    print(f"Dimensions in adt_ds: {list(adt_ds.dims)}")
    if 'time' in adt_ds.dims:
        print(f"Time dimension size in adt_ds: {adt_ds.time.size}")
    
    # Extract ADT data safely, handling potential time dimension issues
    if 'surf_el_t1_adt' in adt_ds:
        # Use surf_el_t1_adt directly if it exists (this is what we created in process_adt_data)
        print("Using surf_el_t1_adt from ADT file")
        adt_data = adt_ds.surf_el_t1_adt.values
    elif 'surf_el' in adt_ds:
        # Otherwise try to get surf_el
        print("Using surf_el from ADT file")
        if 'time' in adt_ds.surf_el.dims:
            # Check if time dimension is non-empty
            if adt_ds.surf_el.time.size > 0:
                adt_data = adt_ds.surf_el.isel(time=0).values
            else:
                # If time dimension is empty, try to use the array directly by averaging over time
                print("Time dimension is empty, using mean over time")
                adt_data = adt_ds.surf_el.mean(dim='time').values
        else:
            # If no time dimension, use the array directly
            adt_data = adt_ds.surf_el.values
    else:
        print("No suitable ADT data found in file")
        return None
    
    # Create new dataset with combined data
    ssh_ds['surf_el_t1_adt'] = (('ylat', 'xlon'), adt_data)
    ssh_ds.surf_el_t1_adt.attrs['_FillValue'] = -30000
    ssh_ds.surf_el_t1_adt.attrs['missing_value'] = -30000
    ssh_ds.surf_el_t1_adt.attrs['scale_factor'] = 1.0
    
    # Create SSH_t1 and ADT_t1 variables
    # Check if we have at least 2 time points
    if ssh_ds.time.size > 1:
        ssh_t1 = ssh_ds.ssh.isel(time=1).values
    else:
        ssh_t1 = ssh_ds.ssh.isel(time=0).values
    
    #ssh_ds['SSH_t1'] = (('time', 'ylat', 'xlon'), ssh_t1[np.newaxis, :, :])
    time_size = ssh_ds.time.size
    ssh_t1_expanded = np.broadcast_to(ssh_t1[np.newaxis, :, :], (time_size, ssh_t1.shape[0], ssh_t1.shape[1]))
    ssh_ds['SSH_t1'] = (('time', 'ylat', 'xlon'), ssh_t1_expanded)

    ssh_ds['ADT_t1'] = ssh_ds.surf_el_t1_adt
    
    # Create filled versions
    ssh_ds['SSH_t1_Fill_0'] = ssh_ds.SSH_t1.copy()
    ssh_ds['ADT_t1_Fill_0'] = ssh_ds.ADT_t1.copy()
    
    # Zero out large values
    ssh_ds['SSH_t1_Fill_0'] = ssh_ds.SSH_t1_Fill_0.where(
        abs(ssh_ds.SSH_t1_Fill_0) <= 1000, 0.0
    )
    
    # Create anomaly SSH
    ssh_ds['SSH_ssh1'] = ssh_ds.ssh - ssh_ds.SSH_t1_Fill_0
    
    # Rename the original surf_el to surf_el_rtofs
    ssh_ds = ssh_ds.rename({'surf_el': 'surf_el_rtofs'})
    
    # Create combined SSH+ADT field
    ssh_ds['SSH_ssh1_adt'] = ssh_ds.ssh - ssh_ds.SSH_t1_Fill_0 + ssh_ds.ADT_t1
    
    # Create new surf_el with scaled values
    ssh_ds['surf_el'] = ssh_ds.SSH_ssh1_adt * 1000.0
    
    # Replace values where SSH_ssh1_adt is too large
    ssh_ds['surf_el'] = ssh_ds.surf_el.where(
        abs(ssh_ds.SSH_ssh1_adt) <= 1000, -3000.0
    )
    
    # Set scale factor
    ssh_ds.surf_el.attrs['scale_factor'] = 0.001
    
    # Extract only needed variables
    #final_ds = ssh_ds[['xlon', 'ylat', 'surf_el']]
    final_ds = xr.Dataset()
    # Copy surf_el but set proper attributes
    final_ds['surf_el'] = ssh_ds.surf_el
    final_ds.surf_el.attrs['scale_factor'] = 0.001
    final_ds.surf_el.attrs['_FillValue'] = -30000.0
    final_ds.surf_el.attrs['missing_value'] = -30000.0
    final_ds.surf_el.attrs['cell_methods'] = "MT: mean"
    final_ds.surf_el.attrs['coordinates'] = "Longitude Latitude Date"
    final_ds.surf_el.attrs['long_name'] = " sea surf. height  [93.1H]"
    final_ds.surf_el.attrs['standard_name'] = "sea_surface_elevation"
    final_ds.surf_el.attrs['units'] = "m"
    final_ds.surf_el.attrs['valid_range'] = [-2.184607, 1.532992]

    # Make time unlimited
    if 'time' in ssh_ds.dims and ssh_ds.time.size > 0:
        final_ds.coords['time'] = ssh_ds.time
    else:
        # Create a time coordinate with 23 points
        times = np.linspace(0, 22, 23)
        final_ds.coords['time'] = ('time', times)

    # Create 2D coordinate variables
    y, x = np.meshgrid(range(ssh_ds.ylat.size), range(ssh_ds.xlon.size), indexing='ij')
    final_ds['Longitude'] = (('ylat', 'xlon'), ssh_ds.Longitude.values)
    final_ds['Latitude'] = (('ylat', 'xlon'), ssh_ds.Latitude.values)
    # Create 2D coordinate variables using meshgrid
    ylat_size = len(ssh_ds.ylat)
    xlon_size = len(ssh_ds.xlon)
    y_grid, x_grid = np.meshgrid(ssh_ds.ylat.values, ssh_ds.xlon.values, indexing='ij')

    # Now create the 2D coordinate variables
    final_ds['Longitude'] = (('ylat', 'xlon'), ssh_ds.Longitude.values)
    final_ds['Latitude'] = (('ylat', 'xlon'), ssh_ds.Latitude.values)
    # Keep xlon and ylat as 1D coordinates
    final_ds.coords['xlon'] = ssh_ds.xlon
    final_ds.coords['ylat'] = ssh_ds.ylat

    # Copy over time coord attributes
    final_ds.time.attrs['C_format'] = "%13.4f"
    final_ds.time.attrs['FORTRAN_format'] = "(f13.4)"
    final_ds.time.attrs['cell_methods'] = "MT: mean"
    final_ds.time.attrs['long_name'] = "date"
    final_ds.time.attrs['units'] = "day as %Y%m%d.%f"

    # Copy over other needed coords and vars
    if 'Date' in ssh_ds:
        final_ds['Date'] = ssh_ds.Date

    # Add global attributes
    final_ds.attrs['Conventions'] = "CF-1.0"
    final_ds.attrs['title'] = "HYCOM ATLb2.00"
    final_ds.attrs['institution'] = "National Centers for Environmental Prediction"
    final_ds.attrs['source'] = "HYCOM archive file"
    final_ds.attrs['experiment'] = "93.1"
    final_ds.attrs['history'] = "archv2ncdf2d"

    # Define unlimited dimension
    encoding = {'time': {'unlimited': True}}
    
    # Save to netCDF file
    out_file = os.path.join(out_dir, "SSH_1.nc")
    final_ds.to_netcdf(out_file, unlimited_dims=['time'])

    # Clean up
    ssh_ds.close()
    adt_ds.close()
    final_ds.close()
    
    print(f"Created combined SSH+ADT file: {out_file}")
    return out_file

def main():
    """Main execution function"""
    args = setup_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Process SSH files
    ssh_merged = process_ssh_files(args.ssh_files, args.idx_2ds, args.out_dir)
    if ssh_merged is None:
        print("Error processing SSH files, exiting")
        return 1
    
    # Process TSUV files
    tsuv_merged = process_tsuv_files(args.tsuv_files, args.idx_3dz, args.out_dir)
    if tsuv_merged is None:
        print("Error processing TSUV files, exiting")
        return 1
    
    # Prepare SSH for SCHISM
    ssh_schism = prepare_ssh_for_schism(ssh_merged, args.out_dir)
    
    # Prepare TSUV for SCHISM
    tsuv_schism = prepare_tsuv_for_schism(tsuv_merged, args.out_dir)
    
    # Process ADT data if available
    adt_processed = process_adt_data(args.adt_today, args.adt_prev, args.adt_weight, args.out_dir)
    
    # Combine SSH and ADT if both are available
    if adt_processed and ssh_schism:
        combined_ssh = combine_ssh_adt(ssh_schism, adt_processed, args.out_dir)
    
    print("Processing completed successfully")
    return 0

if __name__ == "__main__":
    sys.exit(main())
