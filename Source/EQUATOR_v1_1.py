"""
* EQUATOR: Equator Queries simbAd to create Tables of Objects
* Version 1 - 2024-09-21
@ Yaël Moussouni
@ Observatory of Strasbourg 
"""
import astrobs_v1_toolbox as toolbox
import numpy as np
import os

YES = ["YES", "Y", "1"]
QUIT = ["QUIT", "EXIT", "Q"]
OUT_DIR = "./Output/"


# TODO objects --> targets!
# TODO https://astroquery.readthedocs.io/en/latest/imcce/imcce.html#miriade-ephemeris-service

def main():
    stop = False
    targets = toolbox.create_table()
    try:
        files = os.listdir(OUT_DIR)
        files = np.sort([file[:-4] for file in files if file[-4:] == ".cfg"])
        print("\033[32m"
              + "Please select a configurations:"
              + "\033[0m")
        for i in range(len(files)):
            print("\033[32m"
                  + "\t{}. ".format(i+1)
                  + "\033[0m"
                  + files[i])
        selection = int(input("\033[32m"
              + "Choice: "
              + "\033[0m"))-1
        if selection > len(files):
            selection = -1
            print("\033[93m"
                  + ("Warning: Configuration unavailable, "
                     "selected the last one instead.")
                  + "\033[0m")
        config_name = files[selection]
        config = toolbox.read_cfg(config_name)
    except Exception as e:
        config = toolbox.read_cfg(" ")
    while not stop:
        ans = input("\033[32m"
                    + "Select an operation: "
                    + "\033[0m")
        if ans.upper() in QUIT:
            stop = True
        toolbox.resolve_input(ans, targets, config)
        targets.pprint_all()
    return 0

if __name__ == "__main__":
    main()
