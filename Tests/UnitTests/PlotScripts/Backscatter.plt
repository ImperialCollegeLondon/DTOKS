if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF VARIATION OF COEFFICIENTS WITH POTENTIAL AND Ti

set output sprintf("Tests/UnitTests/Plots/BackscatterRE_TiPhi.eps",prefix)
set logscale x
# Plot positions
set xlabel('Ion Temperature (K)')
set ylabel('Potential')
set zlabel('RE (arb)')
splot filename using 3:4:5

set output sprintf("Tests/UnitTests/Plots/BackscatterRN_TiPhi.eps",prefix)
# Plot positions
set zlabel('RN (arb)')
splot filename using 3:4:6

## 2 DIMENSIONAL PLOTS OF VARIATION OF COEFFICIENTS WITH Ti AND Te

set output sprintf("Tests/UnitTests/Plots/BackscatterRE_Te.eps",prefix)
# Plot positions
set xlabel('Electron Temperature (K)')
set ylabel('RE (arb)')
plot filename using 2:5

set output sprintf("Tests/UnitTests/Plots/BackscatterRN_Te.eps",prefix)
# Plot Velocities
set ylabel('RN (arb)')
plot filename using 2:6

set output sprintf("Tests/UnitTests/Plots/BackscatterRE_Ti.eps",prefix)
# Plot positions
set xlabel('Ion Temperature (K)')
set ylabel('RE (arb)')
plot filename using 3:5

set output sprintf("Tests/UnitTests/Plots/BackscatterRN_Ti.eps",prefix)
# Plot Velocities
set ylabel('RN (arb)')
plot filename using 3:6

unset logscale x

## 2 DIMENSIONAL PLOTS OF VARIATION OF COEFFICIENTS WITH POTENTIAL

set output sprintf("Tests/UnitTests/Plots/BackscatterRE_Phi.eps",prefix)
# Plot positions
set xlabel('Potential')
set ylabel('RE (arb)')
plot filename using 4:5

set output sprintf("Tests/UnitTests/Plots/BackscatterRN_Phi.eps",prefix)
# Plot Velocities
set ylabel('RN (arb)')
plot filename using 4:6


