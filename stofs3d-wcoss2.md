# Running STOFS-3D Atlantic on WCOSS2: A Comprehensive Guide (In progress)

## Prerequisites
- Active WCOSS2 account
- Basic familiarity with Linux commands and job scheduling systems

## Model Package Location
The STOFS3D Atlantic package on WCOSS2 development machine is located at:
```bash
/lfs/h1/nos/estofs/noscrub/mansur.jisan/packages/IT-stofs.v2.1.0
```
Location of Operational STOFS on GitHub: 
```bash
https://github.com/noaa-ocs-modeling/STOFS-operational
```
## Configuration Steps

### 1. Script Modifications

#### A. Main Launch Script
Modify `launch_stofs.sh`:
- Replace all instances of `mansur.jisan` with your WCOSS2 username
- Located in the `ecf` directory

#### B. ECFlow Configuration Files
Update the following files in the ECF directory:
- `jstofs_3d_atl_prep_v2.ecf`
- `jstofs_3d_atl_now_forecast_v2.ecf`

In both files:
- Replace all instances of `mansur.jisan` with your WCOSS2 username

#### C. Job Scripts Configuration
Navigate to the `jobs` directory and modify:

1. `JSTOFS_3D_ATL_PREP` script:
   - Replace `mansur.jisan` with your WCOSS2 username
   - Update PDY (Production Day) on line 67:
   ```bash
   export PDY=20250112    # Use a date close to current operational date
   ```

2. `JSTOFS_3D_ATL_NOW_FORECAST` script:
   - Replace `mansur.jisan` with your WCOSS2 username
   - Update PDY on line 81:
   ```bash
   export PDY=20250112    # Use a date close to current operational date
   ```

## Running the Model
create directory and copy rerun file
```bash
cd /lfs/h1/nos/ptmp/mansur.jisan/com/stofs/v2.1/stofs/v2.1/stofs_3d_atl.20250113/rerun
cp /lfs/h1/ops/prod/com/stofs/v2.1/stofs_3d_atl.20250113/rerun/adt_aft_cvtz_cln.nc .
```
### 1. Launch the Process
```bash
cd ecf
./launch_stofs.sh
```

### 2. Monitoring Progress

#### Job Status
Monitor running jobs:
```bash
qstat -u $USER
```

#### Log Files
Track detailed progress in the `log` directory:
- Preprocessing logs: `prep_?????.out`
- Error logs: `prep_?????.err`
- Real-time status updates for individual steps

## Troubleshooting Tips
- Verify all username replacements before running
- Ensure PDY dates are within valid range
- Check log files immediately if jobs fail
- Monitor system resources during run

## Additional Notes
- Keep track of model completion times for future reference
- Regularly check for system maintenance schedules
- Back up important configuration changes

