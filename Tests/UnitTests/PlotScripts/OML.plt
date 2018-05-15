if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF VARIATION OF POTENTIAL WITH TI/TE AND ION CHARGE

set output sprintf("Tests/UnitTests/Plots/OML.eps",prefix)
set logscale x
# Plot positions
set xlabel('Temperature Ratio Ti/Te')
set ylabel('Ion Charge (1/e)')
set zlabel('Normalised Potential (Arb)')
splot filename using 1:2:4
