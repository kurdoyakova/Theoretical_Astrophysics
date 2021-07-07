import numpy as np
import math
import matplotlib.pyplot as plt
from math import sqrt, sin, pi, cos, log10, exp, acos, log
from scipy.integrate import dblquad, quad, romberg

h = 6.626e-34  		#Дж*сек
c = 3.0e8      		#м/сек
k = 1.38e-23   		#Дж*К
T_eff = 5800.0      #K

def solution(nu, mu):
	T_tau1 = lambda tau: T_eff*np.power(tau + 0.5, 0.25)
	T_tau2 = lambda tau: T_eff*np.power(3.0*tau/4.0 + 0.5, 0.25)
	T_tau3 = lambda tau: T_eff*np.power(3.0*tau/4.0 + sqrt(3)/4, 0.25)

	const = 2*h*(nu**3)/(c**2)
	const2 = h*nu/k

	Plank = const/(np.exp(const2/T_eff) - 1)
	S_nu = lambda tau, T, nu: const/(np.exp(const2/T(tau)) - 1)
	I_nu = lambda tau, T, nu: S_nu(tau,T,nu)*np.exp(-tau/mu)/mu

	I1 = quad(I_nu, 0, np.inf, args=(T_tau1, nu), epsabs = 1.5e-32)
	I2 = quad(I_nu, 0, np.inf, args=(T_tau2, nu), epsabs = 1.5e-32)
	I3 = quad(I_nu, 0, np.inf, args=(T_tau3, nu), epsabs = 1.5e-32)
	Eps = [I1[1], I2[1], I3[1]] 
	# print(Eps)
	return I1[0], I2[0], I3[0], Plank, Eps

def first(nu):
	n = 300
	start = 0.00006
	step = (1-start)/n
	mu = [start + i*step for i in range(n+1)]

	I1_0 = np.zeros(n+1); I2_0 = np.zeros(n+1); I3_0 = np.zeros(n+1)
	I1_1 = np.zeros(n+1); I2_1 = np.zeros(n+1); I3_1 = np.zeros(n+1)
	for i in range(n+1):
		I1_0[i], I2_0[i], I3_0[i], non1, eps1 = solution(nu, mu[i])
		I1_1[i], I2_1[i], I3_1[i], non2, eps2 = solution(nu, 1)
	print('nu = ', nu)
	print('eps for ShSh = ', eps1[0])
	print('eps for Edd = ', eps1[1])
	print('eps for Shand = ',eps1[2])
	print('-------------------------')
	fig, ax = plt.subplots()
	plt.plot(mu, I1_0/I1_1, label = 'Шварцшильда-Шустера')
	plt.plot(mu, I2_0/I2_1, label = 'Эддингтона')
	plt.plot(mu, I3_0/I3_1, label = 'Чандрасекара')

	plt.xlabel(r'$\mu$',fontsize=10)
	plt.ylabel(r'$\frac{I_{\nu} (0,\mu)}{I_{\nu} (0,1)}$', fontsize=15)
	ax.set_title(r'$\nu = $' + str(nu/1e14) + r'$^{14}$' + ' Гц')
	ax.legend()
	ax.grid()
	plt.show()

def second():
	n = 120
	start = 1.0e13
	step = (2.0e15-start)/n
	nu = [start + i*step for i in range(n+1)]

	I1 = np.zeros(n+1); I2 = np.zeros(n+1); I3 = np.zeros(n+1); Plank = np.zeros(n+1); eps = np.zeros((n+1,3))
	for i in range(n+1):
		I1[i], I2[i], I3[i], Plank[i], eps[i] = solution(nu[i], 1)
	for i in range(n+1):
		nu[i] = log10(nu[i])

	fig, ax = plt.subplots()
	plt.plot(nu, I1, label = 'Шварцшильда-Шустера')
	plt.plot(nu, I2, label = 'Эддингтона')
	plt.plot(nu, I3, label = 'Чандрасекара')
	plt.plot(nu, Plank, label = 'функция Планка')

	plt.xlabel(r'$\log (\nu) $', fontsize=12 )
	plt.ylabel(r'$I_{\nu}$', fontsize=12)
	ax.legend()
	ax.grid()
	plt.show()

def main():
	first(3.0e14)
	first(5.5e14)
	first(1.0e15)
	second()

main()