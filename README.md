The nested\-grid module part of TIEGCM \(TIEGCM itself is not included in this repository\)


The code design largely follows the original TIEGCM 2\.0 with reduced function calls for simplicity
All netcdf procedures use the Fortran 90 interface and mpi procedures use the f08 interface
The essential differences regarding the design of nested grid are commented in each file
Please see the original TIEGCM codes without \_ng extension for the explanation of physical processes

The additional input parameters for TIEGCM\-NG are explained as follows:
1. NESTING: logical, specifying whether nesting functionality is turned on; default is false \(off\)
2. NSTEP\_NG: integer array, specifying how many sub\-cyclings are used in each level
3. DLEV\_NG: real array, specifying the spacing of vertical grids in each level; it is suggested to have the same vertical grid across all nesting levels
4. GLAT0\_NG,GLAT1\_NG,GLON0\_NG,GLON1\_NG: real array, specifying the top,bottom,left,right edges in each level; require top&lt;bottom,left&lt;right
5. DLAT\_NG,DLON\_NG: real array, specifying the horizontal grid size in each level
6. FILEOUT\_NG: string array, specifying the output filename in each level
7. VAROUT\_NG: string array, specifying the output fields in each level; defaults are TN,O2,O1,Z,ZG \(mandatory output fields\)
8. NUDGE\_NCFILE: string, specifying the filename of external fields; default is empty \(no nudging\)
9. NUDGE\_FLDS: string array, specifying the fields to be nudged; defaults are empty \(no nudging\)
10. NUDGE\_LBC: logical array, specifying whether this field is to be used as the model lower boundary; defaults are false
11. NUDGE\_F4D: logical array, specifying whether this field is to be used to replace the internal model field; defaults are false
12. NUDGE\_LEVEL: logical array, specifying which level \(or global\) to be nudged; the length is the number of nesting levels \+ 1; defaults are false
13. NUDGE\_USE\_REFDATE: logical specifying whether NUDGE\_REFDATE will be used; default is false
	- if NUDGE\_USE\_REFDATE is false, then two integer fields should be included in NUDGE\_NCFILE: date and datesec
	- date is a 8\-digit integer specifying the year,month,day \(YYYYMMDD\), datesec measures the seconds counting from the current day \(00:00\)
	- if NUDGE\_USE\_REFDATE is true, then a integer field should be included in NUDGE\_NCFILE: time, which counts the seconds since NUDGE\_REFDATE
14. NUDGE\_REFDATE: two integers specifying the time of the zero point in NUDGE\_NCFILE; defaults are START\_YEAR,START\_DAY
15. NUDGE\_SPONGE: two real numbers, specifying the horizontal \(in degree\) and vertical \(in scale height\) relaxation span, this should not exceed NUDGE\_NCFILE lon/lat/lev/ilev range; defaults are 0\.0
16. NUDGE\_DELTA: two real numbers, specifying the horizontal \(in degree\) and vertical \(in scale height\) exponential decaying factor; defaults are 1\.0
17. NUDGE\_POWER: two real numbers, specifying the horizontal \(in degree\) and vertical \(in scale height\) power in the exponent; defaults are 1\.0; relaxation function takes the form of exp\(\-\(NUDGE\_SPONGE/NUDGE\_DELTA\)\*\*NUDGE\_POWER\)
18. NUDGE\_ALPHA: real, specifying the temporal relaxation factor; default is 0\.0

In case of erroneous settings, the program will check whether conflict input specifications exist and issue error message on exit

The running procedure is the same as the original TIEGCM
