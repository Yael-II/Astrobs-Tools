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
YES = ["YES", "Y", "1"]
ALL = ["ALL", "*"]
DONE = ["DONE", "OK"]
QUIT = ["QUIT", "EXIT", "Q"]
CANCEL = ["CANCEL", "BACK"]
WRITE = ["WRITE", "SAVE"]
READ = ["READ", "OPEN"]
CALIB = ["CALIB", "CALIBRATION"]
SIMBAD = ["SIMBAD", "OBJECT"]
REGION = ["SEARCH", "REGION"]
ST = ["SIDEREAL", "ST"]
SEQ = ["SEQUENCE", "SEQ"]
CHECK = ["CHECK"]
OUT_DIR = "./Output/"
DEFAULT_CONFIG = """# Default config
LOCATION: None
LAT: 0d00m00s
LON: 0d00m00s
SUN_SET: 1970-01-01 12:00:00.000
SUN_SET_CIVIL: 1970-01-01 12:00:00.000
SUN_SET_NAUTICAL: 1970-01-01 12:00:00.000
SUN_SET_ASTRONOMICAL: 1970-01-01 12:00:00.000
SUN_RISE_ASTRONOMICAL: 1970-01-01 12:00:00.000
SUN_RISE_NAUTICAL: 1970-01-01 12:00:00.000
SUN_RISE_CIVIL: 1970-01-01 12:00:00.000
SUN_RISE: 1970-01-01 12:00:00.000
OBS_BEGIN: 1970-01-01 12:00:00.000
OBS_END: 1970-01-01 12:00:00.000
ST_BEGIN: 12h00m0s
ST_END: 12h00m00s
WINDOW_EAST: 0h00m00s
WINDOW_WEST: 23h59m59s
WINDOW_UPPER: 90d00m00s
WINDOW_LOWER: -90d00m00s
"""

COLS = ["seq", 
        "name", 
        "main_id", 
        "ra",
        "dec", 
        "n_exp", 
        "t_exp", 
        "st_begin", 
        "st_end",
        "notes"]

UNITS = ["", # SEQ
         "", # NAME
         "", # MAIN_ID
         "h", # RA
         "deg", # DEC
         "", # N_EXP
         "s", # T_EXP
         "h", # ST_BEGIN
         "deg", # ST_END
         ""]

TYPES = [np.int_, # SEQ
		 np.str_("<U64"), # NAME
         np.str_("<U64"), # MAIN_ID
         coord.angles.core.Angle, # RA
         coord.angles.core.Angle, # DEC
         np.int_, # N_EXP
         u.quantity.Quantity, # T_EXP
         coord.angles.core.Angle, # ST_BEGIN
         coord.angles.core.Angle, # ST_END
         np.str_("<U512")]

def check_window(table: Table, 
                 config: dict, 
                 i: np.int_ = None) -> np.bool_:
    def check_operation_(line, config):
        if line["main_id"] in SPECIAL:
            return True
        if line["ra"] == "":
            return True
        if line["dec"] == "":
            return True
        ra = coord.Angle(line["ra"])
        dec = coord.Angle(line["dec"])
        east = coord.Angle(config["WINDOW_EAST"])
        west = coord.Angle(config["WINDOW_WEST"])
        upper = coord.Angle(config["WINDOW_UPPER"])
        lower = coord.Angle(config["WINDOW_LOWER"])
        if east <= west:
            if ra < east or ra > west:
                print("\033[93m"
                      +("Warning: The target {} "
                        "is outside the observation window.").format(
                         line["name"])
                      +"\033[0m")
                return False
        else:
            if ra < east and ra > west:
                print("\033[93m"
                      +("Warning: The target {} "
                        "is outside the observation window.").format(
                         line["name"])
                      +"\033[0m")
                return False
            if dec < lower or dec > upper:
                print("\033[93m"
                      +("Warning: The target {} "
                        "is outside the observation window.").format(
                         line["name"])
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
                 types: list = TYPES) -> Table:
    return Table(names=cols, dtype=types, units=units)

def select_obj(table: Table, 
               line: Table,
               config: dict) -> Table:
    table_old = table.copy()
    for i in range(len(line)):
        while line[i]["name"] in line[:i]["name"]: 
            line[i]["name"] = "*"+line[i]["name"]
    line.pprint_all()
    answer = "$"
    print("\033[32m"
          + ("Select the name of the objects "
             "you want to observe (column \"NAME\"):")
          + "\033[0m")
    while answer.upper() not in DONE + [""]:
        answer = input("\033[32m"
                       + "> "
                       + "\033[0m")
        print("\033[A"
              + "\033[32m" 
              + "> " 
              + "\033[0m"
              + "\033[C"*len(answer),
              end="")
        if answer.upper() in DONE + [""]:
            print("\033[34m"
                  + " - ok"
                  + "\033[0m")
        elif answer.upper() in CANCEL or answer.upper() in QUIT:
            i = 0
            while i < len(table):
                if table[i]["name"] not in table_old["name"]:
                    table.remove_row(i)
                    print("remove {}".format(i))
                else: i+= 1
            print("")
            print("\033[34m"
                  + "Operation canceled"
                  + "\033[0m")
            return table
        elif answer in ALL:
            for i in range(len(line)):
                name = line[i]["name"]
                while name in table["name"]: 
                    name = "*" + name
                table.add_row(line[i])
                table[-1]["name"] = name
            print("\033[34m"
                  + " - ok"
                  + "\033[0m")

        elif answer in line["name"]:
            i = np.argwhere(line["name"] == answer)[0][0]
            name = line[i]["name"]
            while name in table["name"]: 
                name = "*" + name
            table.add_row(line[i])
            table[-1]["name"] = name
            print("\033[34m"
                  + " - ok"
                  + "\033[0m")
        else: 
            print("\033[34m"
                  + " - not found"
                  + "\033[0m")
    return table
    
def add_manual(table: Table,
               config: dict,
               seq: np.int_ = 0, 
               name: np.str_ = "", 
               main_id: np.str_ = "",
               ra: coord.angles.core.Angle = None,
               dec: coord.angles.core.Angle = None,
               n_exp: np.int_ = 0, 
               t_exp: u.quantity.Quantity = 0*u.s,
               obj_begin: coord.angles.core.Angle = None,
               obj_end: coord.angles.core.Angle = None,
               notes: np.str_ = "") -> Table:

    while name in table["name"]: name = "*" + name
    while seq > 0 and seq in table["seq"]: seq += 1

    ra_str = ""
    dec_str = ""
    obj_begin_str = ""
    obj_end_str = ""

    if type(ra) == coord.angles.core.Angle:
        ra_str = ra.to_string(unit=u.hour)
    elif type(ra) == str:
        ra_str = ra

    if type(dec) == coord.angles.core.Angle:
        dec_str = dec.to_string(unit=u.degree)
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

def add_calib(table: Table,
              config: dict = None,
              seq: np.int_ = 0,
              name: np.str_ = "CALIB",
              n_exp: np.int_ = 0,
              t_exp: u.quantity.Quantity = 0*u.s,
              notes: np.str_ = "") -> Table:

    while name in table["name"]:
        name = "*" + name

    while seq > 0 and seq in table["seq"]:
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

def add_simbad(table: Table, 
               obj: Table, 
               config: dict,
               seq: np.int_ = 0, 
               name: np.str_ = "", 
               n_exp: np.int_ = 0, 
               t_exp: u.quantity.Quantity = 0*u.s,
               obj_begin: coord.angles.core.Angle = None,
               obj_end: coord.angles.core.Angle = None,
               notes: np.str_ = "") -> Table:
    for i in range(len(obj)):
        main_id = obj["main_id"][i]

        if name == "":
            name = main_id

        while name in table["name"]:
            name = "*" + name

        while seq > 0 and seq in table["seq"]:
            seq += 1

        ra = coord.Angle(obj["ra"][i], unit=u.degree)
        dec = coord.Angle(obj["dec"][i], unit=u.degree)
        ra_str = ra.to_string(unit=u.hour)
        dec_str = dec.to_string(unit=u.degree)
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
        name = ""
    check_window(table, config)
    return table

def read_cfg(filename:str,
             directory:str = OUT_DIR,
             extension:str = ".cfg") -> dict:
    if directory[-1] != "/": directory += "/"
    if extension[0] != ".": extension = "." + extension
    try:
        with open(directory+filename+extension, "r") as file:
            params_index = []
            params_value = []
            for line in file.readlines():
                if len(line) > 0:
                    if line[0] != "#":
                        params = line.strip().split(": ")
                        params_index.append(params[0])
                        params_value.append(params[1])
            config = {params_index[i]:params_value[i] 
                      for i in range(len(params_index))}
    except:
        print("\033[93m"
              + "Warning: No configuration file, using default."
              + "\033[0m")
        if True:
            file = DEFAULT_CONFIG # file is a string
            params_index = []
            params_value = []
            for line in file.split("\n"):
                if len(line) > 0:
                    if line[0] != "#":
                        params = line.strip().split(": ")
                        params_index.append(params[0])
                        params_value.append(params[1])
            config = {params_index[i]:params_value[i] 
                      for i in range(len(params_index))}




    return config

def set_st_window(table, config):
    before = coord.Angle(config["ST_BEGIN"]) \
            - coord.Angle(config["WINDOW_EAST"])
    after = coord.Angle(config["WINDOW_WEST"]) \
            - coord.Angle(config["ST_END"])
    for i in range(len(table)):
        if not table[i]["main_id"] in SPECIAL:
            table[i]["st_begin"] = (coord.Angle(table[i]["ra"]) \
                    - before).wrap_at(24*u.h).to_string()
            table[i]["st_end"] = (coord.Angle(table[i]["ra"]) \
                    + after).wrap_at(24*u.h).to_string()
        else:
            continue
    return table

def set_seq(table, config, stack=False):
    window_east = coord.Angle(config["WINDOW_EAST"]).degree
    window_west = coord.Angle(config["WINDOW_WEST"]).degree
    dupl_seq = np.argwhere(
            np.unique(table["seq"], return_counts=True)[1] > 1).flatten()
    index_dupl = np.argwhere(np.any(
            [table["seq"] == dupl_seq[i] 
             for i in range(len(dupl_seq))], axis=0)).flatten()
    index_zero = np.argwhere(
            table["seq"] == 0).flatten()
    index_special = np.argwhere(np.any(
            [table["main_id"] == SPECIAL[i] 
             for i in range(len(SPECIAL))], axis=0)).flatten()
    index_nora = np.argwhere(table["ra"] == "").flatten()
    index_yes = np.union1d(index_dupl, index_zero)
    index_no = np.union1d(index_special, index_nora)
    index = np.setdiff1d(index_yes, index_no)
    order_table = table[index]
    table["seq"][index] = np.full(len(index), 0)
    order_table.add_columns(
            [coord.Angle(order_table["ra"]).degree], 
            names=["RA_DEG"])
    if window_east <= window_west:
        order_table.sort("RA_DEG")
    else :
        w1 = np.argwhere(
                order_table["RA_DEG"] > (window_east + window_west)/2)
        w2 = np.argwhere(
                order_table["RA_DEG"] <= (window_east + window_west)/2)
        order = np.zeros(len(order_table))
        order[w1] = 1
        order[w2] = 2
        order_table.add_columns(
                [order],
                names=["ORDER"])
        order_table.sort("RA_DEG")
        order_table.sort("ORDER")
    i = 0
    seq = 0
    while i < len(order_table):
        while seq in table["seq"]:
            seq += 1
        order_table["seq"][i] = seq
        i += 1
        seq += 1
    table[index] = order_table
    table.sort("seq")
    while table[0]["seq"] < 0:
        table.add_row(table[0])
        table.remove_row(0)
    if stack:
        i = 0
        while table[i]["seq"] >= 0:
            table[i]["seq"] = i + 1
            i += 1
        i = -1
        while table[i]["seq"] < 0:
            table[i]["seq"] = i
            i -= 1
    return table

def swap(table, name_1, name_2):
    if name_1 not in table["name"] or name_2 not in table["name"]:
        print("\033[93m"
              +"Warning: {} or {} cannot be found in the table.".format(
                  name_1,
                  name_2)
              +"\033[0m")
        return table
    i = 0
    seq_1 = 0
    seq_2 = 0
    while table[i]["name"] not in [name_1, name_2]:
        i += 1
    seq_1 = table[i]["seq"]
    table.add_row(table[i])
    j = i + 1
    while table[j]["name"] not in [name_1, name_2]:
        j += 1
    seq_2 = table[j]["seq"]
    table[i] = table[j]
    table[j] = table[-1]
    table[i]["seq"] = seq_1
    table[j]["seq"] = seq_2
    table.remove_row(-1)
    return table

def write_table(table, 
                config, 
                filename=None, 
                extension=".xml", 
                format_="votable", 
                directory=OUT_DIR):
    data = table.copy()
    if directory[-1] != "/": directory += "/"
    if extension[0] != ".": extension = "." + extension
    if filename == None:
        filename = "{}_{}".format(config["OBS_BEGIN"][:10],
                                  config["LOCATION"])
    try:
        data.write(directory + filename + extension, 
                      format=format_, 
                      overwrite=True)
        print("\033[34m"
              + "File saved in the {} directory ".format(directory)
              + "with the name {}".format(filename+extension)
              + "\033[0m")
    except Exception as e:
        print("\033[91m"
              + ("Error! File cannot be saved "
                 "in the {} directory").format(directory)
              + "\033[0m")

def read_table(filename,
               table,
               config,
               extension=".xml",
               format_="votable",
               directory=OUT_DIR):
    if directory[-1] != "/": directory += "/"
    if extension[0] != ".": extension = "." + extension
    if filename == None:
        filename = "{}_{}".format(config["OBS_BEGIN"][:10],
                                  config["LOCATION"])
    try:
        data = Table.read(directory + filename + extension,
                          format=format_)
        for i in range(len(data)):
            line = data[i]
            if "name" in line.columns: name = line["name"]
            else: name = "{}_{}".format(filename, i)
            if "main_id" in line.columns: main_id = line["main_id"]
            else: main_id = name
            if "ra" in line.columns: 
                if "h" not in str(line["ra"]):
                    ra = coord.Angle("{}d".format(line["ra"])).to_string(unit=u.hour)
                else: ra = line["ra"]
            else: ra = ""
            if "dec" in line.columns: 
                if "d" not in str(line["dec"]):
                    dec = coord.Angle("{}d".format(line["dec"])).to_string(unit=u.degree)
                else: dec = line["dec"]
            else: dec = ""
            if "n_exp" in line.columns: n_exp = line["n_exp"]
            else: n_exp = 0
            if "t_exp" in line.columns: t_exp = line["t_exp"]
            else: t_exp = 0.
            if "st_begin" in line.columns: st_begin = line["st_begin"]
            else: st_begin = ""
            if "st_end" in line.columns: st_end = line["st_end"]
            else: st_end = ""
            if "notes" in line.columns: notes = line["notes"]
            else: notes = ""
            while name in table["name"]:
                name = "*" + name

            for col in line.columns:
                if col not in COLS:
                    notes += " & {}: {}".format(col, line[col])

            table.add_row({"seq": 0,
                           "name": name,
                           "main_id": main_id,
                           "ra": ra,
                           "dec": dec,
                           "n_exp": n_exp,
                           "t_exp": t_exp,
                           "st_begin": st_begin,
                           "st_end": st_end,
                           "notes": notes})

    except Exception as e:
        print("\033[91m"
              + ("Error! File cannot be loaded "
                 "from the {} directory").format(directory)
              + "\033[0m")
        return create_table()
    return table

def resolve_input(text, table, config):
    swap_table = create_table()
    args = text.split(" ")
    try:
        if args[0].upper() in DONE + ["", " "]:
            None
        elif args[0].upper() in WRITE: # TODO ajouter args
            write_table(table, config)
        elif args[0].upper() in READ: # TODO ajouter args
            read_table(args[1], swap_table, config)
            select_obj(table, swap_table, config)
        elif args[0].upper() in SIMBAD:
            name = " ".join(args[1:])
            obj = Simbad.query_object(name)
            swap_table = add_simbad(swap_table, obj, config, name=name)
            select_obj(table, swap_table, config)
        elif args[0].upper() in REGION:
            name = " ".join(args[1]+" "+args[2])
            region = coord.SkyCoord(ra=args[1], dec=args[2])
            radius = coord.Angle(args[3])
            obj = Simbad.query_region(region, radius)
            swap_table = add_simbad(swap_table, obj, config)
            select_obj(table, swap_table, config)
        elif args[0].upper() in CALIB: # TODO ajouter args
            add_calib(table, config)
        elif args[0].upper() in ST:
            set_st_window(table, config)
        elif args[0].upper() in SEQ: # FIXME marche bien une fois mais pas deux ??
            set_seq(table, config)
        elif args[0].upper() in CHECK:
            check_window(table, config) 
        elif args[0].upper() in QUIT:
            None
        else:
            print("\033[93m"
                  + ("Warning: {} is not recognized "
                     "as a valid keyword.").format(args[0])
                  + "\033[0m")
    except Exception as e:
        print("\033[91m"
              + "Error: Your operation is not recognized".format(args[0])
              + "\033[0m")
    return table



