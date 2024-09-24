"""
* Astrobs Toolbox
* Version 1 - 2024-09-21
@ Yaël Moussouni
@ Observatory of Strasbourg 
"""
import numpy as np
import astropy.time as time
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table
from astroquery.simbad import Simbad

SPECIAL = ["CALIBRATION", "IGNORE"]
COLS = ["SEQ", 
        "NAME", 
        "MAIN_ID", 
        "RA",
        "DEC", 
        "N_EXP", 
        "T_EXP", 
        "ST_BEGIN", 
        "ST_END",
        "NOTES"]

UNITS = ["", # SEQ
         "", # NAME
         "", # MAIN_ID
         "hms", # RA
         "dms", # DEC
         "", # N_EXP
         "s", # T_EXP
         "hms", # ST_BEGIN
         "hms", # ST_END
         ""]

TYPES = [np.int32, # SEQ
		 np.str_, # NAME
         np.str_, # MAIN_ID
         coord.angles.core.Angle, # RA
         coord.angles.core.Angle, # DEC
         np.int32, # N_EXP
         u.quantity.Quantity, # T_EXP
         coord.angles.core.Angle, # ST_BEGIN
         coord.angles.core.Angle, # ST_END
         np.str_]

def check_window(table, config, i=None):
    def check_operation_(line, config):
        if line["MAIN_ID"] in SPECIAL:
            return True
        else:
            ra = coord.Angle(line["RA"])
            dec = coord.Angle(line["DEC"])
            east = coord.Angle(config["WINDOW_EAST"])
            west = coord.Angle(config["WINDOW_WEST"])
            upper = coord.Angle(config["WINDOW_UPPER"])
            lower = coord.Angle(config["WINDOW_LOWER"])
            if east <= west:
                if ra < east or ra > west:
                    print("\033[93m"
                          +"Warning: The target {} is outside the observation window!".format(line["NAME"])
                          +"\033[0m")
                    return False
            else:
                if ra < east and ra > west:
                    print("\033[93m"
                          +"Warning: The target {} is outside the observation window!".format(line["NAME"])
                          +"\033[0m")
                    return False
            if dec < lower or dec > upper:
                print("\033[93m"
                      +"Warning: The target {} is outside the observation window!".format(line["NAME"])
                      +"\033[0m")
                return False
        return True
    if i == None:
        output = True
        for i in range(len(table)):
            val = check_operation_(table[i], config)
            if val == False:
                output = False
        return output
    else:
        return check_operation_(table[i], config)

def create_table(cols: list = COLS, 
                 units: list = UNITS,
                 types: list = TYPES) -> list:
    return Table(names=cols, dtype=types, units=units)

def add_manual(table,
               config,
               seq=-1, 
               name="", 
               main_id="",
               ra=None,
               dec=None,
               n_exp=0, 
               t_exp=0*u.s,
               obj_begin=None,
               obj_end=None,
               notes=""):

    while name in table["NAME"]:
        name = "*" + name

    while seq >= 0 and seq in table["SEQ"]:
        seq += 1

    ra_str = ""
    dec_str = ""
    obj_begin_str = ""
    obj_end_str = ""

    if type(ra) == coord.angles.core.Angle:
        ra_str = ra.to_string()
    elif type(ra) == str:
        ra_str = ra

    if type(dec) == coord.angles.core.Angle:
        dec_str = dec.to_string()
    elif type(dec) == str:
        dec_str = dec

    if type(obj_begin) == coord.angles.core.Angle:
        obj_begin_str = obj_begin.to_string()
    elif type(obj_begin) == str:
        obj_begin_str = obj_begin

    if type(obj_end) == coord.angles.core.Angle:
        obj_end_str = obj_end.to_string()
    elif type(obj_end) == str:
        obj_end_str = obj_end

    table.add_row([seq,
                   name,
                   main_id,
                   ra_str,
                   dec_str,
                   n_exp,
                   t_exp.to_value(u.s),
                   obj_begin_str,
                   obj_end_str,
                   notes])
    check_window(table, config, -1)
    return table

def add_calib(table,
              config=None,
              seq=-1,
              name="CALIB",
              n_exp=0,
              t_exp=0*u.s,
              notes=""):

    while name in table["NAME"]:
        name = "*" + name

    while seq >= 0 and seq in table["SEQ"]:
        seq += 1

    table.add_row([seq,
                   name,
                   "CALIBRATION",
                   "",
                   "",
                   n_exp,
                   t_exp.to_value(u.s),
                   "",
                   "",
                   notes])
    check_window(table, config, -1)
    return table

def add_simbad(table, 
               obj, 
               config,
               seq=-1, 
               name="", 
               n_exp=0, 
               t_exp=0*u.s,
               obj_begin = None,
               obj_end = None,
               notes=""):

    main_id = obj["MAIN_ID"][0]

    if name == "":
        name = main_id

    while name in table["NAME"]:
        name = "*" + name

    while seq >= 0 and seq in table["SEQ"]:
        seq += 1

    ra = coord.Angle(obj["RA"][0], unit=u.hour)
    dec = coord.Angle(obj["DEC"][0], unit=u.degree)
    ra_str = ra.to_string()
    dec_str = dec.to_string()
    obj_begin_str = ""
    obj_end_str = ""

    if type(obj_begin) == coord.angles.core.Angle:
        obj_begin_str = obj_begin.to_string()
    elif type(obj_begin) == str:
        obj_begin_str = obj_begin

    if type(obj_end) == coord.angles.core.Angle:
        obj_end_str = obj_end.to_string()
    elif type(obj_end) == str:
        obj_end_str = obj_end

    table.add_row([seq,
                   name,
                   main_id,
                   ra_str,
                   dec_str,
                   n_exp,
                   t_exp.to_value(u.s),
                   obj_begin_str,
                   obj_end_str,
                   notes])
    check_window(table, config, -1)
    return table

def read_cfg(filename:str,
             directory:str = "output/") -> dict:
    with open(directory+filename, "r") as file:
        params_index = []
        params_value = []
        for line in file.readlines():
            if line[0] != "#":
                params = line.strip().split(": ")
                params_index.append(params[0])
                params_value.append(params[1])
        config = {params_index[i]:params_value[i] 
                  for i in range(len(params_index))}
    return config

def set_st_window(objects, config):
    before = coord.Angle(config["ST_BEGIN"]) \
            - coord.Angle(config["WINDOW_EAST"])
    after = coord.Angle(config["WINDOW_WEST"]) \
            - coord.Angle(config["ST_END"])
    for i in range(len(objects)):
        if not objects[i]["MAIN_ID"] in SPECIAL:
            objects[i]["ST_BEGIN"] = (coord.Angle(objects[i]["RA"]) \
                    - before).wrap_at(24*u.h).to_string()
            objects[i]["ST_END"] = (coord.Angle(objects[i]["RA"]) \
                    + after).wrap_at(24*u.h).to_string()
        else:
            continue
    return objects

def set_seq(objects, config):
    # TODO
    return None
    
