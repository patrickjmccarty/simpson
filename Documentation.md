# SIMPSON 4 new features #

## Section `spinsys` ##
`dipole_ave N1 N2 bIS eta alpha beta gamma`   - for motionally averaged dipole-dipole interactions this allows to define asymmetry parameter `eta`

## Section `par` ##
`method value1 value2 value3` ... -  defines calculation methods for `fsimpson`.

|_value_|_feature_|
|:------|:--------|
|`direct`, `idirect`| both values equivalent, uses direct evaluation of `pulseq`, possibly with COMPUTE method if `acq_block` is used|
| `gcompute`, `igcompute` | both values equivalent, uses gamma-COMPUTE method, either in original `pulseq` syntax or in new syntax with `acq_block` – allows for manipulation before acquisition)|
|`time` |forces calculations to be done in time domain|
|`frequency`|forces calculations to be done in frequency domain – can be used only in combination with `acq_block` (both `direct` or `gcompute`)|
|`block_diag`|tries to use block-diagonal Hamiltonians whenever it is possible; considerable speedup for free evolution periods and when rf not applied on all channels|
|`FWTinterpolation`|triggers FWT interpolation method to be used for both amplitudes and frequencies – works with `gcompute` only|
|`FWT2interpolation`|triggers FWT interpolation method to be used just for frequencies, amplitudes are taken the same as in the source orientational set (nearest) – works with `gcompute` only|
|`ASGingterpolation`|triggers ASG interpolation method used after spectrum generation - works with gamma-COMPUTE in frequency domain only|
|`FWTASGinterpolation`, `FWT2ASGinterpolation`|combinations of the interpolation methods|
|`sparse`|uses sparse matix algorithms, needs simpson linked with Intel MKL libraries|
|`diag`, `pade`, `chebyshev`, `taylor`, `lanczos`|experimenting with different methods to calculate propagators (matrix exponentials)|

`num_cores  N`  – defines number of threads to be used for parallel calculations of parameter averages (mainly crystallite orientations). If not defined then all processor cores are used. Number of threads can also be limited by environment variable `SIMPSON_NUM_CORES`

`crystal_file filename N1 N2` – optional arguments N1 and N2 define index of  the first and the last orientation to be read from the filename

`crystal_file filename1 filename2` – used for FWT interpolation, defines source set of orientations (filename1) and larger target set of orientations (filename2). At the moment, only hemisphere Lebedev sets are allowed.

`rfprof_file filename`

`rfprof_file filename N1 N2` – optional arguments N1 and N2 define index of  the first and the last parameter to be read from the filename containing rf inhomogeneity profiles

`averaging_file filename`

`averaging_file filename N1 N2` – optional arguments N1 and N2 define index of  the first and the last parameter to be read from the filename containing set of parameter values over which the calculation should average. The filename is formated as follows: first line defines parameters in the same way as accepted for the `fsimpson` command, other lines define their values. Example:
|shift\_1\_iso|	dipole\_1\_2\_aniso|	weight|
|:------------|:-------------------|:------|
|-10p         |	-10000             |	0.2   |
|+10p         |	-10000             |	0.2   |
|0	           |-10000              |	0.2   |
|0	           |-8000               |	0.2   |
|0	           |-12000              |	0.2   |

`point_per_cycle N` – used for direct method in frequency domain, defines number of points to which the one cycle of acquisition Hamiltonian should be split

`ED_symmetry  val` – controls how to deal with possible excitation-detection symmetry: if val is -1 then it is not considered, if val is 0 then it is internally checked and used, if val is 1 then it is not checked and it is directly considered (user responsibility). Default is -1.

`oc_grad_level N` – defines accuracy level of optimal control gradient, defailt is 1 but L-BFGS methods requires at least 2 or higher

`oc_method` – defines optimal control method, can be CG for conjugated gradients (default), or L-BFGS

## Section `pulseq` ##
`acq_block {…}` – allows easy definition of  pulse manipulations during acquisition. The pulse sequence defined in parenthesis should contain one period of the total Hamiltonian. For gamma-COMPUTE the total length has to be the rotor period divided by the number of  gamma-angles. For COMPUTE the pulse sequence is repeated such that it fills up multiple of rotor periods.