# => makes the plot interactive
# %matplotlib widget 
# inline makes the plots static
#%matplotlib inline 

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from numpy.polynomial import chebyshev
import warnings
# filter astropy warning on fits headers
warnings.filterwarnings('ignore', category=UserWarning, append=True)
from spectro_tools import *

IN_DIR = "./Input/"
OUT_DIR = "./Output/"
LIGHT_DIR = IN_DIR + "Light/"
DARK_DIR = IN_DIR + "Dark/"
FLAT_DIR = IN_DIR + "Flat/"
BIAS_DIR = IN_DIR + "Bias/"
THAR_DIR = IN_DIR + "ThAr/"

def main(ref_file: str = "reduced_master_debiased_ThAr.fits",
         ref_directory: str = THAR_DIR):

    fits_file = ref_directory+ref_file
    xaxis,data=read_raw_spectrum(fits_file)
    spectrum=data

    # Find peaks in the spectrum

    # TODO THRESHOLD

    peaks = find_peaks(data, height=5000.)[0]  # You can adjust the 'height' threshold
    # NB: 'fiducial value': height=5000

    # Get the centroid (x-value) of each peak
    centroid_x_values = peaks
    # Positions in pixels of the peaks

    # Plot the spectrum and mark the centroids
    if False:
        plt.plot(data)
        plt.plot(centroid_x_values, data[peaks], 'ro', label='Max')
        plt.hlines(5000., 0, len(data), "r")
        plt.hlines(np.quantile(data,0.95), 0, len(data), "g")
        plt.xlabel('pixel')
        plt.ylabel('Flux')
        plt.title('Spectrum with peaks and lines')
        plt.grid(True)

        plt.show()
        # Nice, now in order to improve the precision on the line centers,



    # let's fit each detected peak with a gaussian to get a better centroid position
    # generate first_guess for the fitting routine
    # The method below just makes up credible values for a triplet (intensity, centre, width) for each line
    # (~credible) using the peaks detected 
    # and concatenates all that into a large vector first_guess
    first_guess=generate_first_guess(peaks)
    #print(first_guess)

    # fit the lamp spectrum as a sum of gaussian lines using curve_fit and our first guess
    params, covariance = curve_fit(gaussian, xaxis, data, p0=first_guess)
    #print(np.shape(covariance))
    # Reshape params into a 2D array (N, 3) for readability
    num_peaks = len(params) // 3
    params = np.array(params).reshape((num_peaks, 3))
    allamps=params[:,0]
    allcens=params[:,1] # => THIS ARRAY HAS THE FITTED GAUSSIAN CENTROILDS OF THE LINES
    allwids=params[:,2]

    if(0):
        # remove the huge saturaed line at pixel 1987  & 6965 Angstrom
        # well not 100% needed it seems we throw it away later
        print(len(allcens))
        ibad=np.argmin(np.abs(allcens-1987.))
        print(ibad)
        allcens=np.delete(allcens,ibad)
        print(len(allcens))
        allamps=np.delete(allamps,ibad)
        allwids=np.delete(allwids,ibad)
        print(allcens)

        
    # Now plot the spectrum again
    if False:
        plt.plot(data)
        plt.plot(centroid_x_values, data[peaks], 'ro', label='Max')
        plt.xlabel('Pixel')
        plt.ylabel('Flux')
        plt.title('Spectrum with peaks and lines')
        # plot individual gaussian fit for each line, for check
        for i in range(num_peaks):
            fit_params = params[i]  # Extract parameters for each Gaussian
            gau=gaussian(xaxis, *fit_params)
            plt.plot(xaxis, gau)#, label=f'Gaussian {i+1}')
            plt.text(allcens[i], np.max(gau)+3000, str(i), fontsize=12, ha='center', va='center', color='blue')

        plt.legend()
        plt.show()

    fig, ax0 = plt.subplots(1)

    #ax0 = axs[0]
    #ax1 = axs[1]

    ax0.plot(xaxis, data)
    for i in range(num_peaks):
        fit_params = params[i]  # Extract parameters for each Gaussian
        gau=gaussian(xaxis, *fit_params)
        ax0.text(allcens[i], np.max(gau)+3000, str(i), fontsize=10, 
                 ha='center', 
                 va='center', 
                 color='C3')
    #atlas = np.genfromtxt("Source/atlas_linelist.csv", usecols=(0,4), delimiter=",")
    #atlas_val = atlas[:,1]
    #atlas_l = atlas[:,0]
    #w = np.argwhere(atlas_val > 0)
    #atlas_val = atlas_val[w]
    #atlas_l = atlas_l[w]
    #ax1.plot(atlas_l, atlas_val)
    plt.show(block=False)
    ans = "!"
    print("use: https://github.com/pocvirk/astronomical_data_reduction/blob/main/doc/line_atlas_ThAr.pdf")
    print("[number] [wavelength (nm)] ; enter ok when done")
    numbers = []
    lambdas = []
    while not ans == "ok":
        ans = input("> ")
        if ans not in ["ok", "", " "]:
            try:
                n, l = ans.split(" ")
                numbers.append(int(n)), lambdas.append(float(l))
                print("ok !")
            except:
                print("error")
    pixel_lambda = []
    for i in range(len(numbers)):
        pixel_lambda.append([allcens[numbers[i]], lambdas[i]])
    pixel_lambda = np.array(pixel_lambda)
    plt.close()
    # Now derive the full dispersion law as a polynomial fit through the points above
    # Fit a Chebyshev polynomial of degree 1 (linear)
    degree = 1
    coeffs = chebyshev.chebfit(pixel_lambda[:,0], pixel_lambda[:,1], degree)
    # Evaluate the Chebyshev polynomial across xaxis
    y_fit = chebyshev.chebval(xaxis, coeffs)

    print("{},{}".format(coeffs[0], coeffs[1]))
    with open("./Input/values.csv", "w+") as values_file:
        values_files.write("{},{}".format(coeffs[0], coeffs[1]))
    # plot the fit with our calibration points:
    plt.figure(figsize=(5,5))
    plt.scatter(pixel_lambda[:,0],pixel_lambda[:,1])
    plt.xlabel('pixel')
    plt.ylabel('Angstrom')
    plt.plot(xaxis, y_fit, label=f'Chebyshev Polynomial (Degree {degree})', color='red')
    plt.show()

    # thats a pretty good fit.
    # to see how good it is, we will check the residuals in the next cell

    return None

if __name__ == "__main__":
    main()
