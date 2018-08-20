if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF VARIATION OF POTENTIAL WITH GAMMA AND ION CHARGE FOR MOML MODEL

set output sprintf("Tests/UnitTests/Plots/MOML.eps",prefix)
set logscale x
set logscale y
# Plot positions
set xlabel('Temperature Ratio Ti/Te')
set ylabel('Gamma')
set zlabel('Normalised Potential (Arb)')
splot filename using 1:2:4
