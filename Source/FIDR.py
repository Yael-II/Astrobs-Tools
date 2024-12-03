import os
import numpy as np
import data_reduction_tools as tools
import astroalign as aa

TYPE = ["image", "spectroscopy"][1]
target = "NGC_40"
IN_DIR = "./Input/"
OUT_DIR = "./Output/"
LIGHT_DIR = IN_DIR + "Light/"
DARK_DIR = IN_DIR + "Dark/"
FLAT_DIR = IN_DIR + "Flat/"
BIAS_DIR = IN_DIR + "Bias/"
THAR_DIR = IN_DIR + "ThAr/"

def spectro_reduction():
    # BIAS
    filelist = tools.list_fits(BIAS_DIR)
    tools.stack(filelist, BIAS_DIR, "master_bias.fits")

    # DARK
    #filelist = tools.list_fits(DARK_DIR)
    #tools.debias(filelist, DARK_DIR)
    #filelist = ["debiased_" + fname for fname in filelist \
    #        if not "debiased_" in fname]
    #tools.stack(filelist, DARK_DIR, "master_debiased_dark.fits")

    # FLAT
    filelist = tools.list_fits(FLAT_DIR)
    tools.debias(filelist, FLAT_DIR)
    filelist = ["debiased_" + fname for fname in filelist \
            if not "debiased_" in fname]
    tools.normalize_spectrum_flat(filelist)
    filelist = ["normalized_" + fname for fname in filelist \
            if not "normalized_" in fname]
    tools.stack(filelist, FLAT_DIR, "master_normalized_debiased_flat.fits")

    # LIGHT
    filelist = tools.list_fits(LIGHT_DIR)
    tools.debias(filelist, LIGHT_DIR)
    filelist = ["debiased_" + fname for fname in filelist \
            if not "debiased_" in fname]
    tools.stack(filelist, LIGHT_DIR, "master_debiased_light.fits")
    tools.reduction_operation(target + ".fits", 
                        master_light = "master_debiased_light.fits",
                        master_flat = "master_normalized_debiased_flat.fits")
    # MAIN
    # tools.plot_spectrum("master_debiased_flat.fits", FLAT_DIR)
    # tools.plot_spectrum("reduced_" + target + ".fits")
    return 0
        
def spectro_calibration():
    # BIAS
    filelist = tools.list_fits(BIAS_DIR)
    tools.stack(filelist, BIAS_DIR, "master_bias.fits")

    # DARK
    #filelist = tools.list_fits(DARK_DIR)
    #tools.debias(filelist, DARK_DIR)
    #filelist = ["debiased_" + fname for fname in filelist \
    #        if not "debiased_" in fname]
    #tools.stack(filelist, DARK_DIR, "master_debiased_dark.fits")

    # FLAT
    filelist = tools.list_fits(FLAT_DIR)
    tools.debias(filelist, FLAT_DIR)
    filelist = ["debiased_" + fname for fname in filelist \
            if not "debiased_" in fname]
    tools.normalize_spectrum_flat(filelist)
    filelist = ["normalized_" + fname for fname in filelist \
            if not "normalized_" in fname]
    tools.stack(filelist, FLAT_DIR, "master_normalized_debiased_flat.fits")

    # ThAr
    filelist = tools.list_fits(THAR_DIR)
    tools.debias(filelist, THAR_DIR)
    filelist = ["debiased_" + fname for fname in filelist \
            if not "debiased_" in fname]
    tools.stack(filelist, THAR_DIR, "master_debiased_ThAr.fits")
    tools.reduction_operation("ThAr.fits", 
                        light_directory = THAR_DIR,
                        master_light = "master_debiased_ThAr.fits",
                        master_flat = "master_normalized_debiased_flat.fits")
    # MAIN
    tools.plot_spectrum("reduced_ThAr.fits")
    return 0

def image_reduction():
    # BIAS
    filelist = tools.list_fits(BIAS_DIR)
    tools.stack(filelist, BIAS_DIR, "master_bias.fits")

    # DARK
    filelist = tools.list_fits(DARK_DIR)
    tools.debias(filelist, DARK_DIR)
    filelist = ["debiased_" + fname for fname in filelist \
            if not "debiased_" in fname]
    tools.stack(filelist, DARK_DIR, "master_debiased_dark.fits")

    # FLAT
    filter_list = [d for d in os.listdir(LIGHT_DIR) if not d[0] == "."]
    for filter_ in filter_list:
        filter_dir = filter_ + "/"
        filelist = tools.list_fits(FLAT_DIR + filter_dir)
        tools.debias(filelist, FLAT_DIR + filter_dir)
        filelist = ["debiased_" + fname for fname in filelist \
                if not "debiased_" in fname]
        tools.normalize_image_flat(filelist, FLAT_DIR + filter_dir)
        filelist = ["normalized_" + fname for fname in filelist \
                if not "normalized_" in fname]
        tools.stack(filelist, FLAT_DIR + filter_dir, 
                    "master_normalized_debiased_flat_{}.fits".format(filter_))

    # LIGHT
    for filter_ in filter_list:
        filter_dir = filter_ + "/"
        filelist = tools.list_fits(LIGHT_DIR + filter_dir)
        tools.debias(filelist, LIGHT_DIR + filter_dir)
        filelist = ["debiased_" + fname for fname in filelist \
                if not "debiased_" in fname]
        tools.stack(filelist, 
                    LIGHT_DIR + filter_dir, 
                    "master_debiased_light_{}.fits".format(filter_))
        tools.reduction_operation(target + "_" + filter_ + ".fits", 
            master_light = "master_debiased_light_{}.fits".format(filter_),
            light_directory = LIGHT_DIR + filter_dir,
            master_dark = "master_debiased_dark.fits",
            master_flat = "master_normalized_" \
                          + "debiased_flat_{}.fits".format(filter_),
            flat_directory = FLAT_DIR + filter_dir)


    return 0

def main():
    if TYPE == "image":
        image_reduction()
    if TYPE == "spectroscopy":
        spectro_reduction()
        spectro_calibration()
    return 0

if __name__ == "__main__":
    main()
