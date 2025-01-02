import os
import numpy as np
import matplotlib.pyplot as plt
import astroalign as aa
from astropy.io import fits
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

IN_DIR = "./Input/"
OUT_DIR = "./Output/"
LIGHT_DIR = IN_DIR + "Light/"
DARK_DIR = IN_DIR + "Dark/"
FLAT_DIR = IN_DIR + "Flat/"
BIAS_DIR = IN_DIR + "Bias/"
THAR_DIR = IN_DIR + "ThAr/"

PLT_STYLE = "YII_light_1"
if PLT_STYLE in plt.style.available: plt.style.use(PLT_STYLE)

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
    if N > 0:
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
                            d_lim: int = 100):
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
            #normalized_data[np.argwhere(normalized_data<0.5)] = np.nan
            fits.writeto(out_directory + "normalized_" + file, 
                         normalized_data, overwrite=True)
    return 0

def reduction_individual(lightlist: list,
                         light_directory: str = LIGHT_DIR,
                         master_dark: str = None,
                         dark_directory: str = DARK_DIR,
                         master_flat: str = None,
                         flat_directory: str = FLAT_DIR,
                         out_directory: str = None):
    if out_directory == None:
        out_directory = light_directory
    for light in lightlist:
        if not "reduced_" in light:
            data = fits.getdata(light_directory + light)
            head = fits.getheader(light_directory + light)
            if master_dark != None:
                dark_data = fits.getdata(dark_directory + master_dark)
                dark_head = fits.getheader(dark_directory + master_dark)
                exp = head["EXPTIME"]
                dark_exp = dark_head["EXPTIME"]
                data = data - dark_data * exp/dark_exp
                head["history"] = "removed dark with {}".format(master_dark)
            if master_flat != None:
                flat = fits.getdata(flat_directory + master_flat)
                data = data / flat
                head["history"] = "flatten with {}".format(master_flat)
            fits.writeto(out_directory + "reduced_" + light, 
                         data, head, overwrite=True)
    return 0
def plot_spectrum(reduced: str,
                  directory: str = OUT_DIR):
    data = fits.getdata(directory + reduced).flatten()
    head = fits.getheader(directory + reduced)
    pixel = list(range(len(data)))

    plt.plot(pixel, data)
    plt.show(block=True)
    return 0

def align(filelist: str,
          ref_file: str,
          directory: str = LIGHT_DIR,
          ref_dir: str = LIGHT_DIR):
    data_ref = fits.getdata(ref_dir + ref_file)
    good = []
    if len(filelist) > 0:
        for file in filelist:
            data = fits.getdata(directory + file)
            head = fits.getheader(directory + file)
            try:
                aligned_data, fp = aa.register(np.float32(data), np.float32(data_ref))
                good.append(file)
                head["history"] = "aligned with {}".format(ref_file)
                fits.writeto(directory + "aligned_" + file,
                         aligned_data, head,
                         overwrite = True)
            except:
                print("alignement failed with {}".format(file))
    return good

def wavelength_calibration(file: str,
                           out_name: str,
                           directory: str = LIGHT_DIR,
                           reference: str = "values.csv",
                           ref_directory: str = IN_DIR,
                           out_directory: str = OUT_DIR):
    with open(ref_directory+reference) as ref_file:
        vals = ref_file.readline().replace(" ", "").split(",")
    a0 = float(vals[0]) 
    a1 = float(vals[1])
    in_data = fits.getdata(directory + file).flatten()
    head = fits.getheader(directory + file)
    pixel = np.array(range(len(in_data)))
    lambda_ = a1 * pixel + a0
    out_data = np.array([lambda_, in_data])
    np.savetxt(out_directory + out_name, out_data)    
    head["history"] = ("calibrated wavelength: "
                       "lambda = {} + {} * pixel").format(a0, a1)
    fits.writeto(out_directory + out_name.replace(".csv", ".fits"), out_data, head, overwrite=True)
    return 0
