from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm 
import pylab
import numpy as np
from math import exp

def A_func(w):
	return 0.678*w+1.543*w*w-1.212*w*w*w

def i_star_func(z):
	return 1.0-0.0946*z-0.305*z*z+0.950*z*z*z-2.2*z*z*z*z+1.150*z*z*z*z*z


def i_calc(Phi,Beta):
	z = Beta/(1.0+Beta)
	return i_star_func(z)

def A_calc(Phi,Beta,Lambda):
	if Lambda == 0:
		eta = -(Phi/Beta)
	else:
		eta = -(Phi/Beta)*(1+(Beta/4.0)*(1-exp(-4.0*Lambda/(Beta))))
	w = eta/(1.0+eta)

	return A_func(w)

def Arg(Phi,Beta,Lambda):
	return A_calc(Phi,Beta,Lambda)+(1-A_calc(Phi,Beta,Lambda))*i_calc(Phi,Beta)
	
fig = plt.figure()
ax = fig.gca(projection='3d')

x = np.linspace(0,5,100)
y = np.linspace(0,30,100)
X, Y = np.meshgrid( x, y)
i_star = i_star_func(x)
A= A_func(x)
Arg = Arg(-X,Y,0.0)
#pylab.plot(x,i_star)
#pylab.plot(x,A)
surf = ax.plot_surface(X,Y,Arg)
pylab.show()
