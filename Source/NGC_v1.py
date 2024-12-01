import os
from astroquery.vizier import Vizier
from astrobs_v1_toolbox import read_cfg, create_table
import astropy.coordinates as coord
import astropy.units as u

OUT_DIR = "./Output/"
QUIT = ["QUIT", "EXIT", "Q"]

obj_type_list = """\t\033[32m1. \033[0mOpen Cluster
\t\033[32m2. \033[0mGlobular Cluster
\t\033[32m3. \033[0mDiffuse Nebula
\t\033[32m4. \033[0mPlanetary Nebula 9 Object in Small Magellanic Cloud
\t\033[32m5. \033[0mGalaxy
\t\033[32m6. \033[0mCluster associated with nebulosity
\t\033[32m7. \033[0mNon existent
\t\033[32m8. \033[0mObject in Large Magellanic Cloud
\t\033[32m0. \033[0mUnverified southern object
\t\033[32m(blank). \033[0mAny"""


def get_NGC(filename: str,
            obj_type: int = None,
            directory: str = OUT_DIR,
            extension: str = ".cfg"):
    if obj_type == None:
        obj_type = ">0"

    Vizier.clear_cache()
    config = read_cfg(filename, directory, extension)
    date = config["SUN_SET"][:10]
    loc = config["LOCATION"]

    NGC = Vizier(catalog="VII/1B", 
                 columns=["NGC", "Type", "_RAJ2000", "_DEJ2000", "+Mag"])
    objects = NGC.query_constraints(RA1975 = config["RA_CONST"],
                                    DE1975 = config["DE_CONST"],
                                    Mag = "<999",
                                    Type = str(obj_type))[0]
    table = create_table()
    for i in range(len(objects)):
        mag = "Mag: {}".format(objects["Mag"][i])
        ra = "{}".format(coord.Angle(objects["_RAJ2000"][i], 
                         unit=u.degree).to_string(u.hourangle,
                                                  sep=":", 
                                                  pad=True))
        dec =  "{}".format(coord.Angle(objects["_DEJ2000"][i], 
                           unit=u.degree).to_string(u.degree, 
                                                    alwayssign=True, 
                                                    sep=":",
                                                    pad=True))
        name = "NGC {}".format(objects["NGC"][i])
        table.add_row({"seq": 0, 
                       "name": name, 
                       "main_id": name,
                       "ra": ra, 
                       "dec": dec,
                       "notes": mag})
    table.write("{}{}_{}_ngc.xml".format(directory, date, loc), 
                  format="votable",
                  overwrite=True)
    table.pprint()
    return None

def main(directory: str = OUT_DIR,
         extension: str = ".cfg") -> tuple:
    filelist = [fname for fname in os.listdir(directory) if extension in fname]
    print("\033[32m"
          + "Select a configuration"
          + "\033[0m")
    for i in range(len(filelist)):
        print("\033[32m"
              + "\t{}. ".format(i+1)
              + "\033[0m"
              + "{}".format(filelist[i][:-4]))
    try: 
        answer = input("\033[32m" + "Choice: " + "\033[0m")
        if answer in QUIT:
            return 0
        selection = int(answer) - 1
        filename = filelist[selection][:-4]
    except Exception as e:
        print("\033[91m"
              + ("Error! Input value not recognized, "
                 "selected the last location instead")
              + "\033[0m")
        selection = -1
        filename = filelist[selction][:-4]
    print("\033[32m"
          + "Select a type of object"
          + "\033[0m")
    print("\033[32m" + obj_type_list + "\033[32m")
    try:
        answer = input("\033[32m" + "Choice: " + "\033[0m")
        if answer in QUIT:
            return 0
        if answer in [" ", ""]:
            obj_type = None
        obj_type = int(answer)
    except Exception as e:
        obj_type = None

    get_NGC(filename, obj_type)
    return 0

if __name__ == "__main__":
    main()
