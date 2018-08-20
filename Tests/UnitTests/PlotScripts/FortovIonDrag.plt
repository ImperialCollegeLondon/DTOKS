if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF VARIATION OF DRAG FORCE WITH FLOW VELOCITY AND TE/TI

set output sprintf("Tests/UnitTests/Plots/FortovIonDrag_UTeTi.eps",prefix)
set logscale x
set logscale y
set logscale z
set format z "%s*10^{%S}"
set xlabel('Normalised ion flow velocity, u_i (1/v_T)')
set ylabel('Temperature Ratio Te/Ti')
set zlabel('Ion Drag Force (N)')
splot filename using 1:4:6

## 3 DIMENSIONAL PLOTS OF VARIATION OF DRAG FORCE WITH TE/TI AND POTENTIAL

#set output sprintf("Tests/UnitTests/Plots/FortovIonDrag_TeTiPhi.eps",prefix)
#set logscale x
#unset logscale y
#set xlabel('Temperature Ratio Te/Ti')
#set ylabel('Normalised potential (arb)')
#set zlabel('Ion Drag Force (N)')
#splot filename using 4:5:6
