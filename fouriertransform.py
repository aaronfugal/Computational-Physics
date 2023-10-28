import numpy as np
import matplotlib.pyplot as plt

filename = '.txt' #insert text file with data where first column in time and column 2 is intensity of signal
data = np.loadtxt(filename,skiprows=0,usecols=(0,1,2))

# time (s)
# mean intensity of time
# standard deviation of the mean

# Construct workable data here


time = data[:,0]
intensity = data[:,1]


#14854


zeros = np.zeros(14854)

intensity_ = np.append(intensity,zeros)

N=len(intensity_)
dt = time[1]-time[0] #sampling interval

transform = np.fft.fft(intensity_)
freq = np.fft.fftfreq(N,dt)
df = freq[1]-freq[0]


transform_shift = np.fft.fftshift(transform)
freq_shift = np.fft.fftshift(freq)

# Visually inspect data here


fig = plt.figure(figsize=(6,6))

fig.add_subplot(311)
plt.plot(time,intensity)
fig.add_subplot(212)
plt.plot(freq_shift,transform_shift)
plt.show()

# fig = plt.figure(figsize = (6,6))

# fig.add_subplot(311)
# plt.plot(freq_shifted,transform_shift.real, c = 'b', label = 'real FFT')