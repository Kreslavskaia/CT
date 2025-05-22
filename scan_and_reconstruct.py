from ct_scan import *
from ct_calibrate import *
from ct_lib import *
from ramp_filter import *
from back_project import *
from hu import *

def scan_and_reconstruct(
    photons,
    material,
    ph,
    scale=0.1,
    angles=256,
    mas=1000,
    alpha=0.01):

	""" Simulation of the CT scanning process
		reconstruction = scan_and_reconstruct(photons, material, phantom, scale, angles, mas, alpha)
		takes the phantom data in phantom (samples x samples), scans it using the
		source photons and material information given, as well as the scale (in cm),
		number of angles, time-current product in mas, and raised-cosine power
		alpha for filtering. The output reconstruction is the same size as phantom."""
	print("I am workign correct")

	# convert source (photons per (mas, cm^2)) to photons
	photons = photons*mas*scale**2

	# create sinogram from phantom data, with received detector values
	sino_raw = ct_scan(photons*mas*scale**2, material, ph, scale, angles)
	draw(np.log1p(sino_raw))

	# convert detector values into calibrated attenuation values
	sino_cal = ct_calibrate(photons*mas*scale**2, material, sino_raw, scale)
	draw(sino_cal)

	# Ram-Lak
	sino_filt = ramp_filter(sino_cal, scale, alpha)
	draw(sino_filt)

	# Back-projection
	reco = back_project(sino_filt)
	#Clip negative values
	reco[reco<0] = 0
	draw(reco)

	# convert to Hounsfield Units
	#Not yet

	return reco
