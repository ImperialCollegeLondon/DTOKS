if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 2 DIMENSIONAL PLOTS OF VARIATION OF DELTATherm WITH Te

set output sprintf("Tests/UnitTests/Plots/DeltaTherm_Te.eps",prefix)
set logscale y
# Plot positions
set xlabel('Electron Temperature (eV)')
set ylabel('DeltaTherm (arb)')
plot filename using 1:2


