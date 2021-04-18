import numpy as np
from astropy.io import fits
import glob
import sep
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.table import Table
import extinction

calibrated_dir = "../phys391-data/RRLeo_calibrated/"
sky_dir = "../phys391-data/RRLeo_sky/"
nosky_dir = "../phys391-data/RRLeo_nosky/"
test_sky_subtraction_dir = "../phys391-data/RRLeo_test_sky_subtraction/"


def calibration():
    print("Calibrating...")
    raw_dir_a = "../phys391-data/BKZ_group/RR_Leo_a_11.03.21/"
    raw_dir_b = "../phys391-data/BKZ_group/RR_Leo_b_11.03.21/"

    hdu = fits.open("../phys391-data/MasterBias.fit")  # open the master bias file
    MasterBias = hdu[0].data.astype("float")  # take the data, convert to floats, and assign it to MasterBias
    hdu = fits.open("../phys391-data/MasterDark.fit")  # open the master dark file
    MasterDark = hdu[0].data.astype("float")  # take the data, convert to floats, and assign it to MasterDark
    hdu = fits.open("../phys391-data/2021.03.10-FLAT-MASTERS-bin2/MasterFlat_PhotV.fit")  # open FlatFileName
    MasterFlat = hdu[0].data.astype("float")  # take the data, convert to floats, and assign it to MasterFlat

    filelist_a = glob.glob(raw_dir_a + "*")
    filelist_b = glob.glob(raw_dir_b + "*")
    file_index = 1

    for f in filelist_a:
        hdu = fits.open(f)
        rawdata = hdu[0].data.astype("float")  # get the data, convert to float, and assign to rawdata
        texp = hdu[0].header["EXPTIME"]  # get the exposure time

        # Use the equation in lecture 4 to calibrate the data
        calibrated_data = (rawdata - texp * MasterDark - MasterBias) / MasterFlat

        # modify the hdu file, but save it in the calibrated directory
        hdu[0].data = calibrated_data
        hdu[0].header["CLBRATED"] = (True, "Calibrated from file " + raw_dir_a + f)
        hdu[0].header["BZERO"] = 0
        hdu[0].header["BSCALE"] = 1.0

        filename = "RR_Leo-" + str(file_index).zfill(4) + "_PhotV.fit"
        hdu.writeto(calibrated_dir + filename, overwrite=True)
        file_index += 1

    for f in filelist_b:
        hdu = fits.open(f)
        rawdata = hdu[0].data.astype("float")  # get the data, convert to float, and assign to rawdata
        texp = hdu[0].header["EXPTIME"]  # get the exposure time

        # Use the equation in lecture 4 to calibrate the data
        calibrated_data = (rawdata - texp * MasterDark - MasterBias) / MasterFlat

        # modify the hdu file, but save it in the calibrated directory
        hdu[0].data = calibrated_data
        hdu[0].header["CLBRATED"] = (True, "Calibrated from file " + raw_dir_a + f)
        hdu[0].header["BZERO"] = 0
        hdu[0].header["BSCALE"] = 1.0

        filename = "RR_Leo-" + str(file_index).zfill(4) + "_PhotV.fit"
        hdu.writeto(calibrated_dir + filename, overwrite=True)
        file_index += 1

    print("Calibration done.")


def subtractSky(image_data, file_index, filter_size=3, box_size=64):
    def makeplots(image_data, sky_data, image_data_nosky, file_index):
        def colourbar(sc, ax):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(sc, cax=cax, orientation="vertical")

        # Plot the 3 images side by side to compare them (the original image, the sky image, and the image minus the sky)
        fig, ax = plt.subplots(1, 3, figsize=[12, 4])

        norm = ImageNormalize(image_data, interval=ZScaleInterval())  # scales the image by same 'zscale' algorithm as ds9
        im0 = ax[0].imshow(image_data, origin="lower", cmap="gray", norm=norm)
        ax[0].set_title("Original")

        im1 = ax[1].imshow(sky_data, origin="lower", cmap="gray")  # linear.  no need for zscale
        ax[1].set_title("Sky")

        norm = ImageNormalize(image_data_nosky, interval=ZScaleInterval())  # scales the image by same 'zscale' algorithm as ds9
        im2 = ax[2].imshow(image_data_nosky, origin="lower", cmap="gray", norm=norm)
        ax[2].set_title("After sky subtraction")

        # Remove the ticks and tick labels
        for a in ax:
            a.xaxis.set_visible(False)
            a.yaxis.set_visible(False)

        # Add colour bars to all three panels (not as simple when using subplots; calls function below)
        colourbar(im0, ax[0])
        colourbar(im1, ax[1])
        colourbar(im2, ax[2])
        fig.tight_layout()
        fig.savefig(test_sky_subtraction_dir + "test_sky_subtraction_" + str(file_index).zfill(4) + ".png")
        plt.close(fig)

    # Determine the sky level as a function of position
    sky = sep.Background(image_data, fw=filter_size, fh=filter_size, bh=box_size, bw=box_size)  # use SEP to determine the background sky level
    sky_data = sky.back()  # This is the 2D array containing the sky level in each pixel (ADUs)

    image_data_nosky = image_data - sky_data  # Subtract the sky data from the image data and assign the result to image_data_sub

    makeplots(image_data, sky_data, image_data_nosky, file_index)  # Call the function makeplots defined below

    return sky_data, image_data_nosky, sky.globalrms


def computeNoise(image_data_nosky, sky_data, texp):
    # 1. Declare the gain and read noise (in electrons)
    gain = 1.5
    readnoise = 19.0

    # 2. Read the master dark data from Assignment 4 and compute the Poision noise (in electrons)
    hdu = fits.open("../phys391-data/MasterDark.fit")
    MasterDark = hdu[0].data.astype("float")
    noise_dark = np.mean((np.absolute(gain * MasterDark)) ** 0.5)

    # 3. Compute the Poisson noise of the sky in electrons
    noise_sky = (np.absolute(gain * sky_data)) ** 0.5

    # 4. Compute the Poisson noise of the sky-subtracted image in electrons
    noise_data = (np.absolute(gain * image_data_nosky)) ** 0.5

    # 5. Total noise: add the image noise, sky noise, dark noise and read noise in quadrature.  Note that the dark noise defined above
    # is the dark noise per second.  Make sure to scale it up by the exposure time.
    noise = (noise_data ** 2 + noise_sky ** 2 + (noise_dark) ** 2 * texp + readnoise ** 2) ** 0.5

    # 6. Turn the noise back into units of ADU
    noise /= gain

    return noise  # Return the noise in units of ADU


def sourceExtraction(image_data_nosky, skyrms):
    # Use sep.extract to find all the objects in image_data_nosky that are 2-sigma above the background
    objects = sep.extract(image_data_nosky, 2, err=skyrms)

    # Get the dimensions of image_data_nosky
    ny, nx = image_data_nosky.shape

    mask = (objects['x'] > 750) & (objects['x'] < nx-750) & (objects['y'] > 750) & (objects['y'] < ny-750)  # pixel location of RR Leo
    objects = objects[mask]  # Overwritten "objects" with a subset that fulfill the above criteria

    return objects  # Return the structured array of objects


def referenceSourceExtraction(image_data_nosky, skyrms):
    # Use sep.extract to find all the objects in image_data_nosky that are 2-sigma above the background
    objects = sep.extract(image_data_nosky, 2, err=skyrms)

    # Get the dimensions of image_data_nosky
    ny, nx = image_data_nosky.shape

    mask = (objects['x'] > 450) & (objects['x'] < 550) & (objects['y'] > 1025) & (objects['y'] < 1125)  # pixel location of refernce object
    objects = objects[mask]  # Overwritten "objects" with a subset that fulfill the above criteria

    return objects  # Return the structured array of objects


def getMags(image_data_nosky, objects, noise, texp):
    # This function first computes the kron radius and then finds the flux within the 2.5*kron_radius.
    # Follow the examples in the slides to fill this out.
    kronrad, krflag = sep.kron_radius(image_data_nosky, objects['x'], objects['y'], objects['a'], objects['b'], objects['theta'], 6.0)
    flux, fluxerr, flag = sep.sum_ellipse(image_data_nosky, objects['x'], objects['y'], objects['a'], objects['b'], objects['theta'], 2.5*kronrad, subpix=5, err=noise)

    index = np.argmax(flux)
    flux, fluxerr, flag = flux[index], fluxerr[index], flag[index]

    # Compute the S/N as the flux divided by the error on the flux
    sn = flux / fluxerr

    # Now compute the magnitudes and magnitude error
    x = flux / texp
    xerr = fluxerr / texp
    mag = -2.5 * np.log10(x) + 0
    mag_err = np.absolute(-2.5 * (1 / (x * np.log(10)))) * xerr

    # Correct the magnitudes for atmospheric extinction
    # mag -= ext_coeff * airmass

    return mag, mag_err, sn, flag  # Return the magnitudes, magnitude errors, S/N and the flag


ref_apparent_mag = 11.02  # GSC 1968:912

mags = []
mags_err = []
times = []

calibration()
ext_coeff = extinction.extinction()
print("ext_coeff: " + str(ext_coeff))
filelist = glob.glob(calibrated_dir + "*")
file_index = 1
for f in filelist:
    print("Processing fit file: " + str(file_index).zfill(4))
    hdu = fits.open(f)
    image_data = hdu[0].data  # get the data, and assign to image_data (after calibration they're already floats)
    image_data = image_data.byteswap(inplace=True).newbyteorder()  # This line is need for SEP; otherwise it produces an error
    texp = hdu[0].header["EXPTIME"]
    airmass = hdu[0].header["AIRMASS"]
    sky_data, image_data_nosky, skyrms = subtractSky(image_data, file_index, filter_size=10, box_size=64)

    hdu[0].header["GLOBLRMS"] = skyrms
    hdu[0].data = sky_data  # Re-use the header information from the old hdu, but overwrite the data with the sky data
    hdu.writeto(sky_dir + "sky_" + str(file_index).zfill(4) + ".fit", overwrite=True)  # save the hdu as new file

    hdu[0].data = image_data_nosky  # Re-use the header information from the old hdu, but overwrite the data with the sky-subtracted image
    hdu.writeto(nosky_dir + "nosky_" + str(file_index).zfill(4) + ".fit", overwrite=True)  # save the hdu as new file

    noise = computeNoise(image_data_nosky, sky_data, texp)
    objects = sourceExtraction(image_data_nosky, skyrms)
    print(Table(objects))
    mag, mag_err, sn, flag = getMags(image_data_nosky, objects, noise, texp)
    referenceObjects = referenceSourceExtraction(image_data_nosky, skyrms)
    print(Table(referenceObjects))
    ref_mag, ref_mag_err, ref_sn, ref_flag = getMags(image_data_nosky, referenceObjects, noise, texp)
    mag_offset = ref_mag - ref_apparent_mag
    apparent_mag = mag - mag_offset
    mags.append(apparent_mag)
    mags_err.append(mag_err)
    times.append(texp * file_index)

    file_index += 1
    print()

# Amplitude of collected data i.e. the half the difference between the max and min
amp = (max(mags) - min(mags)) / 2

plt.errorbar(times, mags, yerr=mags_err, fmt=".")
plt.title("RR Leo luminosity period relationship")
plt.ylabel("apparent magnitude")
plt.gca().invert_yaxis()  # invert the y-axis because we want "bright" to be at the top
plt.xlabel("time (seconds)")
text = "amplitude: " + str(np.round(amp, 5))
text += "\n"
text += "ext_coeff: " + str(np.round(ext_coeff, 5))
text += "\n"
period = file_index / 60.0 / 24
text += "observed period: " + str(np.round(period, 7)) + " d"
text += "\n"
text += "published period: 0.4524021 d"
plt.text(0.1, 0.07, text, transform=plt.gca().transAxes)
# plt.show()
plt.savefig("RRLeo_period.eps")

print("Maximum apparent magnitude: " + str(np.round(max(mags), 5)))
print("Minimum apparent magnitude: " + str(np.round(min(mags), 5)))
