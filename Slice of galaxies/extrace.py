from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from math import sqrt, log10, log, isnan
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS
import numpy as np

def func(x, c1, c2):
	# Линейная аппроксимация
    return c1*x + c2 

def Approximation(x,y, filename):
	x, y = zip(*sorted(zip(x,y)))
	popt, pcov = curve_fit(func, x, y)
	# Обязательно находим ошибки аппроксимации
	perr = np.sqrt(np.diag(pcov))   
	plt.plot(x, y, '-', label='HDF-N  ' + filename[5:-10], color='purple')
	plt.plot(x, func(np.asarray(x), *popt), 'g--',
         label='fit: c1=%5.5f, c2=%5.5f' % tuple(popt))
	plt.ylim(max(y)+0.5, min(y)-0.5)
	plt.xlabel(r'$r_e$ ', fontsize=12)
	plt.ylabel(r'$\mu_e$',fontsize=12)
	plt.legend()
	plt.show()
	return popt, perr

def main(file_name, z, scale, n, vvmin, vvmax, center):
	image_data = fits.getdata(file_name)
	image_data[image_data == 'nan'] = vvmin
	image_data[image_data <= 0 ] = vvmin
	plt.imshow(image_data, cmap='gray', norm=LogNorm(vmin = vvmin, vmax = vvmax))
	plt.title('HDF-N  ' + file_name[5:-10] + '  with z=' + str(z))
	# Рисуем линию разреза
	plt.hlines(center, 0, len(image_data[n])-1, color = 'r') 
	plt.colorbar()
	plt.show()

	maxx_value, maxx = 0, 0
	axe = np.zeros(len(image_data[n]))
	rad = np.zeros(len(image_data[n]))
	rad2 = np.zeros(len(image_data[n]))
	#  Находим максимальное значеие для центровки по нему
	for i in range(len(image_data[n])): 
		if maxx_value < image_data[n][i]:
			maxx_value = image_data[n][i]
			maxx = i
	# Центрируем разрез 
	for i in range(len(image_data[n])):
		axe[i] = i - maxx


	# maxx_value, maxx = 0, 0
	# for i in range(len(mu_real)): 
	# 	if maxx_value < mu_real[i]:
	# 		maxx_value = mu_real[i]
	# 		maxx = i
	# for i in range(len(mu_real)):
	# 	axe2[i] = i - maxx

	# plt.plot(axe2, mu_real, color = 'green')
	# plt.plot(axe, image_data[n], color = 'red')
	# plt.show()

	# Отражаем кривую относительно центра в положительную (правую) сторону        
	for i in range(maxx):
		axe[i] = sqrt(axe[i]**2)  
	# Строим в указанных координатах
	for i in range(len(axe)):
		rad[i] = axe[i] **(1/4)
		# if rad[i] > 0:
		# 	rad2[i] = log(rad[i])
		# # else:
			# rad2[i] = -3        # Не помню зачем такое писала, пусть будет

	# Можно посмотреть на разрез, устредненый и отраженный
	plt.plot(rad, image_data[n], color = 'green')
	plt.show()
	# Тот же разрез в логарифмическом масштабе
	# plt.plot(rad2, image_data[n], color = 'green')
	# plt.show()
	# plt.clf()
	# Находим звездные величины
	# mag = np.zeros(len(mu_real))
	mag = np.zeros(len(image_data[n]))
	for i in range(len(image_data[n])):
		# mag2 = (image_data[n][i] + image_data[n+1][i])/2   + (image_data[n][i+1] + image_data[n+1][i+1])/2
		# mag[i] = 22.08 - 2.5*log10(mag2) - 2.5*log10((1 + z)**3)
		mag[i] = 22.08 - 2.5*log10(image_data[n][i]) - 2.5*log10((1 + z)**3) -2.5*log10(1/(0.04*0.04))
		# mag[i] = 22.08 - 2.5*log10(mu_real[i]) - 2.5*log10((1 + z)**3)


	# Давайте выделим область для аппроксимации, выберу наиболее прямой промежуток на всем срезе
	rad1 = rad[rad < 2]
	radd = rad1[rad1 > 1.1]
	magg = mag[rad < 2]
	magg = magg[rad1 > 1.1]
	# rad1 = rad2[rad2 < 0.33]
	# radd = rad1[rad1 > 0.1]
	# magg = mag[rad2 < 0.33]
	# magg = magg[rad1 > 0.1]


	# Если хотите увидеть разрез целиком, то раскоменьте эту строку, а следующую заокмментируете
	# c, eps = Approximation(rad2, mag, file_name)
	c, eps = Approximation(radd, magg, file_name)
	# Находим значения радиуса и звездных величин из формулы (12) в ссылке пдф
	mu_e = c[1] + 8.3268
	r_e = (8.3268/c[0])**4
	# Домножаем на расчитанный ранее масштаб на сайте, указанного в пдф
	r_e = r_e * scale
	return r_e, mu_e

def main_for_all(file_name, z, scale, vvmin, vvmax, center):
	# Тут собираем все данные и строим итогувую картиночку
	r_e, mu_e = [], []
	all_par=np.zeros((len(file_name),2))
	for gal in range(len(file_name)):
		# cutss, n = cuts(file_name[gal], z[gal], vvmin, vvmax)   # Для обрезки кадра онлайн без смс и регистрации
		cutss = 'cuts_' + str(file_name[gal])                     #  ДЛя уже готовых обрезанных кадров
		all_par[gal] = main(cutss, z[gal], scale[gal], int(center[gal]), vvmin, vvmax, center[gal])
		r_e.append(log10(all_par[gal][0]))
		mu_e.append(all_par[gal][1])
	print('r_e = ', r_e, 'mu = ', mu_e)

	# mu_e2 = np.zeros(len(r_e))
	# for j in range(int(max(r_e)) ):
	# 	mu_e2[j] = 3* log10(r_e[j]) + 20.1
	# x = np.linspace(-3, 3,100)
	# mu_e2 = lambda x: 3* log10(x) + 20.1
	# mu = np.vectorize(mu_e2)


	# fig, ax = plt.subplots()
	# popt, pcov = curve_fit(func, r_e, mu_e)
	# perr = np.sqrt(np.diag(pcov))   
	# ax.scatter(r_e, mu_e, label='Галлактики глубокого поля Хаббла', color='purple')
	# ax.plot(r_e, func(np.asarray(r_e), *popt), 'g--', label='fit: c1=%5.5f, c2=%5.5f' % tuple(popt))
	# ax.set_ylim(max(mu_e)+0.5, min(mu_e)-0.5)
	# ax.set_xlabel(r'$\lg r_e$  [kpc]', fontsize=12)
	# ax.set_ylabel(r'$\mu_e$ ',fontsize=12)
	# ax.plot(r_e,mu(r_e), color='blue')
	# # print(r_e, mu_e2)
	# ax.legend()
	# plt.show()

	# c, eps = Approximation(r_e, mu_e, '_____for all__________')
	# plt.scatter(r_e, mu_e)
	# plt.gca().invert_yaxis()
	# plt.show()

def image_save_cutout(filename, position, size):
    # Load the image and the WCS
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)
    # Make the cutout, including the WCS
    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)
    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data
    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())
    # Write the cutout to a new FITS file
    cutout_filename = 'cuts_' + str(filename)
    hdu.writeto(cutout_filename, overwrite=True)
    return str(cutout_filename)

def cuts(filename, z, vvmin, vvmax):
	coordinate_x = []
	coordinate_y = []
	# Будем тыкать на кадр для обрезки. Выбираем сначала левый верхний угол, потом правый нижний и закрываем окно
	def onclick(event):
		print('x=%f, y=%f'%(event.xdata, event.ydata))
		coordinate_x.append(event.xdata)
		coordinate_y.append(event.ydata)
	image_data = fits.getdata(filename)
	# image_data[image_data == 'nan'] = vvmin
	image_data[image_data <= 0 ] = vvmin
	fig, ax = plt.subplots()
	plt.title('choose right region for ' + filename + '  with z=' + str(z))
	ax.imshow(image_data, cmap='gray', norm=LogNorm(vmin = vvmin, vmax = vvmax))
	# plt.colorbar()
	fig.canvas.mpl_connect('button_press_event', onclick)
	plt.show()
	plt.clf()

	# Вычисляем размеры получившегося кадра для его сохранения и его центр для будующего среза
	# Либо можно вручную выбрать центр и ввести в переменную center
	position = [(coordinate_x[1] - coordinate_x[0])/2 + coordinate_x[0], (coordinate_y[1] - coordinate_y[0])/2 + coordinate_y[0]]
	size = [coordinate_y[1] - coordinate_y[0], coordinate_x[1] - coordinate_x[0]]
	n = size[0]/2
	file_name = image_save_cutout(filename, position, size)
	return file_name, n

if __name__ == "__main__":
	# Вводим все нащи данные, подготовленные заранее
	filename = '17_true.fits', '122_true.fits', '124_true.fits', '125_true.fits', '273_true.fits', '303_true.fits', '495_true.fits', '524_true.fits', '619_true.fits', '653_true.fits'
	z =     [1.013, 0.764, 0.504, 0.562, 0.680, 1.000, 0.880, 0.678, 0.370, 0.600]
	scale = [8.147, 7.482, 6.196, 6.553, 7.148, 8.122, 7.844, 7.139, 5.168, 6.762]
	vvmax, vvmin = 0.005, 9e-06
	center = [48, 50, 28, 43, 19, 22, 13, 41, 30, 28]
	main_for_all(filename, z, scale, vvmin, vvmax, center)