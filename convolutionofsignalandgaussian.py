#Convolution combines two functions, one should be a signal in a time domain and the other in our case is a gaussian distribution centered at 0. The result is a continuous signal

import numpy as np
import matplotlib.pyplot as plt

filename1 = "gaussgen.dat" # should have a time column and an amplitude column. 
filename2 = 'signal.dat' # should have the same number of indicies as gaussian file and same time entries. 

data1 = np.loadtxt(filename1,skiprows=0,usecols=(0,1))
data2 = np.loadtxt(filename2,skiprows=0,usecols=(0,1))

time1  = data1[:,0]
time2  = data2[:,0]
gauss  = data1[:,1]
signal  = data2[:,1]



ymax2 = 1.1*np.max(np.abs(signal))
ymin2 = -ymax2
xmax2 = np.max(time2)
xmin2 = np.min(time2)

signal_hat = np.fft.fft(signal)
gauss_hat = np.fft.fft(gauss)

signal_hat_shift = np.fft.fftshift(signal_hat)
gauss_hat_shift = np.fft.fftshift(gauss_hat)

y_hat = signal_hat_shift * gauss_hat_shift
y= np.fft.ifft(y_hat)

y_shift = np.fft.fftshift(y)

ymax_fft = 1.1*np.max(np.abs(y_shift))
ymin_fft = -ymax_fft


# Plotting the original signal and Gaussian
plt.figure(figsize=(12, 6))
plt.subplot(3, 1, 1)
plt.plot(time2, signal, label='Original Signal')
plt.title('Original Signal')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(time1, gauss, label='Gaussian')
plt.title('Gaussian Function')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

# Plotting the Convolution Result
plt.subplot(3, 1, 3)
plt.plot(time2, np.abs(np.real(y_shift)), label='Convolution Result')
plt.title('Convolution Result')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

plt.tight_layout()
plt.show()
