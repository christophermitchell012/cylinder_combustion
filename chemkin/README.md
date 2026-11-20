# C5 equilibrium with CHEMKIN-II

The Fortran driver calls CKGML for standard molar Gibbs energies. The input deck uses the two NASA7 pentane records (X10/85) in Cantera 3.2.0's nasa_gas.yaml. It contains zero reactions because this comparison only needs thermodynamic properties.

Build with GNU Fortran from this directory. Set CK_SRC to the CHEMKIN source directory. The tested [source snapshot](https://github.com/siddarthc/CHOMBO-EBAMRReactive/tree/0885f23e3dccf8c2fb87c84637727d5b6d61363f/CHEMKIN) includes the April 1994 double-precision interpreter 3.6 and library 4.9.

```sh
gfortran -std=legacy "$CK_SRC/ckinterp.f" -o ckinterp
gfortran -std=legacy -I"$CK_SRC" -c "$CK_SRC/cklib.f" -o cklib.o
gfortran piona.f90 cklib.o -o piona
./ckinterp
echo 453.15 | ./piona
python compare.py
```

The interpreter creates chem.bin and chem.out. Its warning about missing reaction input is expected. The Fortran executable needs only chem.bin at runtime. The optional Python comparison requires NumPy and Cantera 3.2.0. Feed multiple temperatures on separate lines to obtain CSV. The 80% approach and 10% blend loading match the Python defaults. The wide sweep checks polynomial evaluation only.
