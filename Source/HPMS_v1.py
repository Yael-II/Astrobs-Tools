import os
from astroquery.simbad import Simbad
from astrobs_v1_toolbox import read_cfg

OUT_DIR = "./Output/"
QUIT = ["QUIT", "EXIT", "Q"]

def get_HPMS(filename: str,
             directory: str = OUT_DIR,
             extension: str = ".cfg",
             pm_max: float = 800.0, 
             pm_min: float = 30.0,
             mag_max: float = 8.5,
             mag_min: float = 13.0) -> tuple:
    config = read_cfg(filename, directory, extension)
    date = config["SUN_SET"][:10]
    constraint = config["CONSTRAINT"]
    loc = config["LOCATION"]

    query_hpms_high = """SELECT TOP 20
           main_id,
           ra,
           dec,
           pmra,
           pmdec,
           SQRT(pmra*pmra + pmdec*pmdec) as "PM",
           filter,
           flux
    FROM basic JOIN flux ON oid=oidref
    WHERE {where}
            AND SQRT(pmra*pmra + pmdec*pmdec) > {pm}
            AND flux > {mag_max} AND flux < {mag_min}
            AND filter='V'
    ORDER BY "PM" DESC;
    """.format(where = constraint, 
               pm = pm_max, 
               mag_max = mag_max, 
               mag_min = mag_min)

    query_hpms_low = """SELECT TOP 20
           main_id,
           ra,
           dec,
           pmra,
           pmdec,
           SQRT(pmra*pmra + pmdec*pmdec) as "PM",
           filter,
           flux
    FROM basic JOIN flux ON oid=oidref
    WHERE {where}
            AND SQRT(pmra*pmra + pmdec*pmdec) > {pm}
            AND flux > {mag_max} AND flux < {mag_min}
            AND filter='V'
    ORDER BY "PM" ASC;
    """.format(where = constraint, 
               pm = pm_min, 
               mag_max = mag_max, 
               mag_min = mag_min)
    
    hpms_high = Simbad.query_tap(query_hpms_high)
    hpms_low = Simbad.query_tap(query_hpms_low)

    hpms_high.write("{}{}_{}_hpms_high.xml".format(directory, date, loc), 
                    format="votable", 
                    overwrite=True)
    hpms_low.write("{}{}_{}_hpms_low.xml".format(directory, date, loc), 
                   format="votable", 
                   overwrite=True)
    print("\033[36m" + "High Proper Motion Stars: Highest 20" + "\033[0m")
    hpms_high.pprint()
 
    print("\033[36m" + "High Proper Motion Stars: Limiting 20" + "\033[0m")
    hpms_low.pprint()
    return 0

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
    get_HPMS(filename)
    return None
    
if __name__ == "__main__":
    main()
