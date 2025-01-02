import os
import numpy as np
import data_reduction_tools as tools
import astroalign as aa

TYPE = ["image", "spectroscopy"][1]
IN_DIR = "./Input/"
OUT_DIR = "./Output/"
LIGHT_DIR = IN_DIR + "Light/"
DARK_DIR = IN_DIR + "Dark/"
FLAT_DIR = IN_DIR + "Flat/"
BIAS_DIR = IN_DIR + "Bias/"
THAR_DIR = IN_DIR + "ThAr/"

def spectro_reduction(target: str):
    # BIAS
    filelist = tools.list_fits(BIAS_DIR)
    tools.stack(filelist, BIAS_DIR, "master_bias.fits")

    # DARK
    filelist = tools.list_fits(DARK_DIR)
    has_dark = False
    if len(filelist) > 0:
        has_dark = True
        tools.debias(filelist, DARK_DIR)
        filelist = ["debiased_" + fname for fname in filelist \
                if not "debiased_" in fname]
        tools.stack(filelist, DARK_DIR, "master_debiased_dark.fits")

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
    if has_dark:
        tools.reduction_individual(["master_debiased_light.fits"],
                        master_dark = "master_debiased_dark.fits",
                        master_flat = "master_normalized_debiased_flat.fits")
    else:
        tools.reduction_individual(["master_debiased_light.fits"],
                        master_flat = "master_normalized_debiased_flat.fits")

    # MAIN
    # tools.plot_spectrum("master_debiased_flat.fits", FLAT_DIR)
    # tools.plot_spectrum("reduced_" + target + ".fits")
    return 0
        
def spectro_calibration(target: str):
    # BIAS
    filelist = tools.list_fits(BIAS_DIR)
    tools.stack(filelist, BIAS_DIR, "master_bias.fits")

    # DARK
    filelist = tools.list_fits(DARK_DIR)
    has_dark = False
    if len(filelist) > 0:
        has_dark = True
        tools.debias(filelist, DARK_DIR)
        filelist = ["debiased_" + fname for fname in filelist \
                if not "debiased_" in fname]
        tools.stack(filelist, DARK_DIR, "master_debiased_dark.fits")

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
    if has_dark:
        tools.reduction_individual(["master_debiased_ThAr.fits"], 
                        light_directory = THAR_DIR,
                        master_dark = "master_debiased_dark.fits",
                        master_flat = "master_normalized_debiased_flat.fits")
    else:
        tools.reduction_individual(["master_debiased_ThAr.fits"], 
                        light_directory = THAR_DIR,
                        master_flat = "master_normalized_debiased_flat.fits")
    # MAIN
    filename = "reduced_master_debiased_light.fits"
    tools.wavelength_calibration(filename,
                                 "calibrated_" + target + ".csv")
    return 0

def image_reduction(target: str):
    filter_list = [d for d in os.listdir(LIGHT_DIR) if not d[0] == "."]
    # BIAS
    filelist = tools.list_fits(BIAS_DIR)
    tools.stack(filelist, BIAS_DIR, "master_bias.fits")

    # DARK
    dark_list = [d for d in os.listdir(DARK_DIR) if not d[0] == "."]
    filelist = tools.list_fits(DARK_DIR)
    has_dark = False
    has_filter_dark = False
    if np.array_equal(np.sort(dark_list), np.sort(filter_list)):
        has_dark = True
        has_filter_dark = True
        for filter_ in filter_list:
            filter_dir = filter_ + "/"
            filelist = tools.list_fits(DARK_DIR + filter_dir)
            tools.debias(filelist, DARK_DIR + filter_dir)
            filelist = ["debiased_" + fname for fname in filelist \
                    if not "debiased_" in fname]
            tools.stack(filelist, DARK_DIR + filter_dir,
                        "master_debiased_dark_{}.fits".format(filter_))

    elif len(filelist) > 0:
        has_dark = True
        tools.debias(filelist, DARK_DIR)
        filelist = ["debiased_" + fname for fname in filelist \
                if not "debiased_" in fname]
        tools.stack(filelist, DARK_DIR, "master_debiased_dark.fits")




    # FLAT
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
    ref_dir = filter_list[0] + "/"
    ref_file = tools.list_fits(LIGHT_DIR + ref_dir)[0]
    ref_file = "reduced_debiased_" + ref_file
    good = []
    bad = []
    for filter_ in filter_list:
        filter_dir = filter_ + "/"
        filelist = tools.list_fits(LIGHT_DIR + filter_dir)
        tools.debias(filelist, LIGHT_DIR + filter_dir)
        filelist = ["debiased_" + fname for fname in filelist \
                if not "debiased_" in fname]
        if has_dark:
            if has_filter_dark:
                tools.reduction_individual(filelist, 
                        light_directory = LIGHT_DIR + filter_dir,
                        master_dark = "master_debiased_dark_{}.fits".format(filter_),
                        master_flat = "master_normalized_" \
                                      + "debiased_flat_{}.fits".format(filter_),
                        dark_directory = DARK_DIR + filter_dir,
                        flat_directory = FLAT_DIR + filter_dir)
            else:
                tools.reduction_individual(filelist, 
                        light_directory = LIGHT_DIR + filter_dir,
                        master_dark = "master_debiased_dark.fits",
                        master_flat = "master_normalized_" \
                                      + "debiased_flat_{}.fits".format(filter_),
                        flat_directory = FLAT_DIR + filter_dir)
        else: 
            tools.reduction_individual(filelist, 
                    light_directory = LIGHT_DIR + filter_dir,
                    #master_dark = "master_debiased_dark.fits",
                    master_flat = "master_normalized_" \
                                  + "debiased_flat_{}.fits".format(filter_),
                    flat_directory = FLAT_DIR + filter_dir)

        filelist = ["reduced_" + fname for fname in filelist \
                if not "reduced_" in fname]
        align = tools.align(filelist, 
                           ref_file, 
                           LIGHT_DIR + filter_dir, 
                           LIGHT_DIR + ref_dir)
        if len(align) > 1:
            good.append(filter_)
            filelist = ["aligned_" + fname for fname in align \
                    if not "aligned_" in fname]
            tools.stack(filelist, 
                    LIGHT_DIR + filter_dir, 
                    "master_aligned_reduced_debiased_light_{}_{}.fits".format(
                        target, filter_), out_directory = OUT_DIR)
        else:
            bad.append(filter_)
            tools.stack(filelist, 
                    LIGHT_DIR + filter_dir, 
                    "master_reduced_debiased_light_{}_{}.fits".format(
                        target, filter_), out_directory = OUT_DIR)
    if len(good) > 0 and len(bad) > 0:
        filter_ref = good[0]
        filter_ref_file = [fname for fname in tools.list_fits(OUT_DIR)
                           if "aligned_" in fname][0]
        print(filter_ref, filter_ref_file)
        print(good)
        print(bad)
        for filter_ in bad:
            tools.align(["master_reduced_debiased_light_{}_{}.fits".format(
                            target, filter_)],
                        filter_ref_file, 
                        OUT_DIR,
                        OUT_DIR)
    return 0

def main():
    print("Select the type of reduction")
    print("\t1. Images")
    print("\t2. Spectra")
    i = int(input("Choice [1 or 2]: "))
    if i == 1:
        print(("IMPORTANT: Please put raw files in the Input/ directory "
               "(with any name)"))
        print("Input/")
        print("  ├ Bias/")
        print("  │  └ bias_1.fits, ...")
        print("  ├ Dark/")
        print("  │  └ dark_1.fits, ...")
        print("  ├ Flat/")
        print("  │  ├ B/")
        print("  │  │ └ flat_B_1.fits, ...")
        print("  │  └ V/...")
        print("  └ Light/")
        print("     ├ B/")
        print("     │ └ raw_B_1.fits")
        print("     └ V/ ...")
    if i == 2:
        print(("IMPORTANT: Please put raw files in the Input/ directory "
               "(with any name)"))
        print("Input/")
        print("  ├ Bias/")
        print("  │  └ bias1.fits, ...")
        print("  ├ Dark/ (optional)")
        print("  │  └ dark_1.fits, ...")
        print("  ├ Flat/")
        print("  │  └ flat_1.fits")
        print("  ├ Light/")
        print("  │  └ raw_1.fits")
        print("  └ ThAr/")
        print("     └ thar_raw_1.fits")
    name = input("Target name: ")
    if i == 1:
        image_reduction(name)
    if i == 2:
        spectro_reduction(name)
        spectro_calibration(name)
    return 0

if __name__ == "__main__":
    main()
