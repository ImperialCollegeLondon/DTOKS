if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF VARIATION OF POTENTIAL WITH TI/TE

set output sprintf("Tests/UnitTests/Plots/DTOKScharging_Td_TiTe.eps",prefix)
set logscale y
# Plot positions
set xlabel('Dust Temperature (K)')
set ylabel('Temperature Ratio Ti/Te')
set zlabel('Normalised Potential (Arb)')
splot filename using 1:2:4

set output sprintf("Tests/UnitTests/Plots/DTOKScharging_dt_TiTe.eps",prefix)
set logscale x
set xlabel('DeltaTot (1/I_{th})')
set ylabel('Temperature Ratio Ti/Te')
set zlabel('Potential (Arb)')
splot filename using 3:2:4
