The nested-grid module part of TIEGCM (TIEGCM itself is not included in this repository)

The code design largely follows the original TIEGCM 2.0 with reduced function calls for simplicity.
All netcdf procedures use the Fortran 90 interface and mpi procedures use the f08 interface.
The essential differences regarding the design of nested grid are commented in each file.
Please see the original TIEGCM codes without _ng extension for the explanation of physical processes.

The additional input parameters for TIEGCM-NG are explained as follows:
1. NESTING: logical, specifying whether nesting functionality is turned on; default is false (off)
2. NSTEP_NG: integer array, specifying how many sub-cyclings are used in each level
3. DLEV_NG: real array, specifying the spacing of vertical grids in each level
4. GLAT0_NG,GLAT1_NG,GLON0_NG,GLON1_NG: real array, specifying the top,bottom,left,right edges in each level; require top<bottom,left<right
5. DLAT_NG,DLON_NG: real array, specifying the horizontal grid size in each level
6. FILEOUT_NG: string array, specifying the output filename in each level; defaults are empty (no output)
7. VAROUT_NG: string array, specifying the output fields in each level; defaults are empty (no output)
8. NUDGE_NCPRE,NUDGE_NCPOST: string, specifying the common part of NUDGE_NCFILE filenames; default is empty
9. NUDGE_NCFILE: string array, specifying the varying part of NUDGE_NCFILE filenames; default is empty (no nudging), NUDGE_NCPRE+NUDGE_NCFILE+NUDGE_NCPOST together specify NUDGE_NCFILE filenames
10. NUDGE_FLDS: string array, specifying the fields to be nudged; defaults are 'TN','UN','VN','Z'
11. NUDGE_LBC,NUDGE_F4D: logical, specifying whether LBC or internal fields will be replaced with NUDGE_NCFILE fields; defaults are false
12. NUDGE_PERT: logical, specifying whether NUDGE_NCFILE fields are deviation fields; default is false
13. NUDGE_LEVEL: logical array, specifying which level (or global) to be nudged; the length is the number of nesting levels + 1; defaults are false
14. NUDGE_SPONGE: two real numbers, specifying the horizontal (in degree) and vertical (in scale height) relaxation span, this should not exceed NUDGE_NCFILE lon/lat/lev/ilev range; defaults are 0.0
15. NUDGE_DELTA: two real numbers, specifying the horizontal (in degree) and vertical (in scale height) exponential decaying factor; defaults are 1.0
16. NUDGE_POWER: two real numbers, specifying the horizontal (in degree) and vertical (in scale height) power in the exponent; defaults are 0.0; relaxation function takes the form of exp(-(NUDGE_SPONGE/NUDGE_DELTA)**NUDGE_POWER)
17. NUDGE_ALPHA: real, specifying the temporal relaxation factor; default is 0.0
18. NUDGE_TSHIFT: integer, specifying the time shift (second) of NUDGE_NCFILE (rarely used); default is 0

In case of erroneous settings, the program will check whether conflict input specifications exist and issue error message on exit.

The running procedure is the same as the original TIEGCM.