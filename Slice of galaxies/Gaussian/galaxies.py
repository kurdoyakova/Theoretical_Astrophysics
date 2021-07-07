#!/usr/bin/env python3

from scipy.optimize import curve_fit
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import os
import shutil


	# Функция Гаусса с помощью которой будем аппроксимировать срезы
def gauss(x, a, b, c):                  
	y = a*np.exp(-(x-b)**2/(2*c*c))
	return(y)

print('step')    
step=int(input())

print('filename')
filename=str(input())

data = fits.getdata(filename) # Открываем фитсы
x = np.linspace(0, 100, data.shape[0])

shutil.rmtree('plots')
os.mkdir('plots')    # Создаём дирректорию в которой у нас будут все картиночки
os.chdir('plots')

# Здесь будут храниться параметы функции Гаусса для каждого среза
a_list = []    
b_list = []
c_list = []
x2 = []

# curve_fit вычисляет наилучшие параметры для Гауссиана
#	...которые после добавляются в лист параметров
for i in range(1, data.shape[1]+1, step):
	params = curve_fit(gauss, x, data[:,i], maxfev=1000000,
	                   p0=(data[:,i].max(),50.0,10.0)) 
	x2.append(i)
	a_list.append(params[0][0])
	b_list.append(params[0][1])
	c_list.append((params[0][2])**2)
		
	ourslice = plt.scatter(x, data[:,i]) 
	gs = plt.plot(x,gauss(x,params[0][0],params[0][1],params[0][2])) # Наш Гауссиан
	s=str(i)
	plt.savefig(s)  # Сохраняем картинку
	plt.close()

for i in range(len(a_list)):  # Если параметры бесконечны, исправляем
	if abs(a_list[i]) > 100*abs(np.median(np.array(a_list))):
		a_list[i] = 0
	if abs(b_list[i]) > 100*abs(np.median(np.array(b_list))):
		b_list[i] = 0
	if abs(c_list[i]) > 100*abs(np.median(np.array(c_list))):
		c_list[i] = 0



plt.plot(a_list, color = 'r', label = "a, амплитуда")
plt.plot(b_list, color = 'g', label = "b, сдвиг от 0")
plt.plot(c_list, color = 'b', label = "c, полуширина")
plt.legend()
plt.savefig('Итоговая')
plt.close()

plt.plot(x2, a_list)    #graphics of parameters
plt.savefig('amplitude')
plt.close()

plt.plot(x2, b_list)
plt.savefig('shift from 0')
plt.close()

plt.plot(x2, c_list)
plt.savefig('half-width')
plt.close()
