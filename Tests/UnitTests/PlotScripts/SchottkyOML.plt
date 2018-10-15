if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF VARIATION OF POTENTIAL WITH DUST TEMPERATURE AND ION TEMPERATURE FOR SchottkyOML MODEL

set output sprintf("Tests/UnitTests/Plots/SchottkyOML_Ti.eps",prefix)
set logscale y
# Plot positions
set xlabel('Dust Temperature (K)')
set ylabel('Ion Temperature (K)')
set zlabel('Normalised Potential (Arb)')
splot filename using 1:2:4

## 3 DIMENSIONAL PLOTS OF VARIATION OF POTENTIAL WITH DUST TEMPERATURE AND ELECTRON TEMPERATURE FOR SchottkyOML MODEL

set output sprintf("Tests/UnitTests/Plots/SchottkyOML_Te.eps",prefix)
set logscale y
# Plot positions
set xlabel('Dust Temperature (K)')
set ylabel('Electron Temperature (K)')
set zlabel('Normalised Potential (Arb)')
splot filename using 1:3:4
