import math
import numpy as np
import numpy.matlib

def ramp_filter(sinogram, scale=0.1, alpha=0.001):
	""" Ram-Lak filter with raised-cosine for CT reconstruction

	fs = ramp_filter(sinogram, scale) filters the input in sinogram (angles x samples)
	using a Ram-Lak filter.

	fs = ramp_filter(sinogram, scale, alpha) can be used to modify the Ram-Lak filter by a
	cosine raised to the power given by alpha."""

	# get input dimensions
	angles = sinogram.shape[0]
	n = sinogram.shape[1]

	#Set up filter to be at least twice as long as input
	m = np.ceil(np.log(2*n-1) / np.log(2)) 
	m = int(2 ** m)

	# Frequency axis, normalized to [-1, 1]
	freqs = np.fft.fftfreq(m, d=scale)
	norm_freqs = freqs / np.max(np.abs(freqs))

	# Create Ram-Lak filter
	filter_ramlak = np.abs(freqs)

	# Apply raised-cosine taper
	if alpha > 0:
		filter_ramlak *= np.cos(np.pi * norm_freqs / 2) ** alpha

	print('Ramp filtering')

	# Allocate output array
	fs = np.zeros_like(sinogram)

	# Apply filter to each angle
	for i in range(angles):
		# Zero-pad projection to length m
		proj = np.zeros(m)
		proj[:n] = sinogram[i, :]

		# FFT, filter, IFFT
		proj_fft = np.fft.fft(proj)
		filtered_fft = proj_fft * filter_ramlak
		filtered = np.real(np.fft.ifft(filtered_fft))

		# Truncate back to original length
		fs[i, :] = filtered[:n]

	return fs
