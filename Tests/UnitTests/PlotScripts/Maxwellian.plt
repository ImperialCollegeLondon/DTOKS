if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF VARIATION OF MAXWELLIAN WITH ENERGY AND TEMPERATURE

set output sprintf("Tests/UnitTests/Plots/Maxwellian3D.eps",prefix)
set logscale x
set logscale y
set logscale z
set xlabel('Energy (eV)')
set ylabel('Temperature (eV)')
set zlabel('Maxwellian (Arb)')
splot filename using 1:2:3

## 2 DIMENSIONAL PLOTS OF VARIATION OF MAXWELLIAN WITH ENERGY AND TEMPERATURE

set output sprintf("Tests/UnitTests/Plots/Maxwellian2D.eps",prefix)
set xlabel('Energy (eV)')
set zlabel('Maxwellian (Arb)')
plot filename using 1:3
