import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


IN_DIR = "./Input/"
OUT_DIR = "./Output/"
LIGHT_DIR = IN_DIR + "Light/"
DARK_DIR = IN_DIR + "Dark/"
FLAT_DIR = IN_DIR + "Flat/"
BIAS_DIR = IN_DIR + "Bias/"
THAR_DIR = IN_DIR + "ThAr/"

PLT_STYLE = "YII_light_1"

def list_fits(directory: str = LIGHT_DIR):
    """
    Returns the list of files in the directory
    @params:    
        - directory: the directory containing the fits
    @output:
        - list of fits files in the directory
    """
    filelist = [fname for fname in os.listdir(directory) if ".fit" in fname] 
    # works for fit and fits
    return filelist

def stack(filelist: list, 
          in_directory: str, 
          out_name: str,
          out_directory: str = None,
          method:str = "median"):
    in_data = []
    if out_directory == None:
        out_directory = in_directory
    N = len(filelist)
    for file in filelist:
        if not "master_" in file:
            in_data.append(fits.getdata(in_directory + file))
    in_data = np.array(in_data)
    if method == "mean":
        out_data = np.mean(in_data, axis=0)
    if method == "median":
        out_data = np.median(in_data, axis=0)
    header = fits.getheader(in_directory+file)
    header["history"] = "stacking with {} files".format(N)
    fits.writeto(out_directory \
            + out_name, out_data, header, overwrite=True)
    return 0

def debias(filelist: list, 
           in_directory: str, 
           master_bias: str = "master_bias.fits",
           bias_directory: str = BIAS_DIR,
           out_directory: str = None):
    if out_directory == None:
        out_directory = in_directory
    for file in filelist:
        if not "debiased_" in file:
            in_data = fits.getdata(in_directory + file)
            header = fits.getheader(in_directory + file)
            bias = fits.getdata(bias_directory+master_bias)
            out_data = in_data - bias
            header["history"] = "debiased with master bias"
            fits.writeto(out_directory + "debiased_" + file, \
                         out_data, header, \
                         overwrite=True)
    return 0

def normalize_image_flat(filelist: str, 
                         in_directory: str = FLAT_DIR,
                         out_directory: str = None):
    if out_directory == None:
        out_directory = in_directory

    for file in filelist:
        if not "normalized_" in file:
            data = fits.getdata(in_directory + file)
            head = fits.getheader(in_directory + file)
            normalized_data = data/np.median(data)
            fits.writeto(out_directory + "normalized_" + file, 
                         normalized_data, overwrite=True)
    return 0


def normalize_spectrum_flat(filelist: str, 
                            in_directory: str = FLAT_DIR,
                            out_directory: str = None,
                            fit_degree: int = 3, 
                            d_lim: int = 10000):
    if out_directory == None:
        out_directory = in_directory

    poly = np.polynomial.chebyshev
    for file in filelist:
        if not "normalized_" in file:
            data = fits.getdata(in_directory + file).flatten()
            head = fits.getheader(in_directory + file)
            pixel = list(range(len(data)))
            mask = np.ones_like(data)
            mask[np.argwhere(data < d_lim)] = 0
            fit_coefs = poly.chebfit(pixel, data, fit_degree, w=mask)
            fit_data = poly.chebval(pixel, fit_coefs)
            normalized_data = data/fit_data
            normalized_data[np.argwhere(normalized_data<0.5)] = np.nan
            fits.writeto(out_directory + "normalized_" + file, 
                         normalized_data, overwrite=True)
    return 0

def reduction_operation(target: str,
                        master_light: list,
                        light_directory: str = LIGHT_DIR,
                        master_dark: str = None,
                        dark_directory: str = DARK_DIR,
                        master_flat: str = None,
                        flat_directory: str = FLAT_DIR,
                        out_directory: str = OUT_DIR):

    data = fits.getdata(light_directory + master_light)
    head = fits.getheader(light_directory + master_light)
    if master_dark != None:
        dark = fits.getdata(dark_directory + master_dark)
        data = data - dark
        head["history"] = "removed dark with {}".format(master_dark)
    if master_flat != None:
        flat = fits.getdata(flat_directory + master_flat)
        data = data / flat
        head["history"] = "flatten with {}".format(master_flat)
    fits.writeto(out_directory + "reduced_" + target, 
                 data, head, overwrite=True)


    return 0

def plot_spectrum(reduced: str,
                  directory: str = OUT_DIR):
    if PLT_STYLE in plt.style.available: plt.style.use(PLT_STYLE)
    data = fits.getdata(directory + reduced).flatten()
    head = fits.getheader(directory + reduced)
    pixel = list(range(len(data)))

    plt.plot(pixel, data)
    plt.show(block=True)
    return 0

def find_peaks_highest(data: list,
                       N: int = 10):
    """Find the N highest peaks in the data"""
    peaks = []
    mask = np.
    for i in range(N):
        None
    return 0





