"""
* Astrobs Toolbox
* Version 1 - 2024-09-21
@ Yaël Moussouni
@ Observatory of Strasbourg 
"""
import numpy as np
import matplotlib.pyplot as plt
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
CANCEL = ["CANCEL", "BACK", "UNDO"]
WRITE = ["WRITE", "SAVE"]
READ = ["READ", "OPEN", "LOAD"]
CALIB = ["CALIB", "CALIBRATION"]
SIMBAD = ["SIMBAD", "OBJECT"]
REGION = ["SEARCH", "REGION"]
MANUAL = ["MANUAL", "ADD"]
ST = ["SIDEREAL", "ST"]
SEQ = ["SEQUENCE", "SEQ"]
CHECK = ["CHECK"]
PLOT = ["PLOT", "GRAPH"]
HELP = ["HELP", "H", "?"]

FIGURE_DPI = 100
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
WINDOW_EAST: 00:00:00.0
WINDOW_WEST: 23:59:59.9
WINDOW_UPPER: +90:00:00.0
WINDOW_LOWER: -90:00:00.0
CONSTRAINT: ' ' LIKE ' '
"""

COLS = ["seq", 
        "name", 
        "main_id", 
        "ra",
        "dec", 
        "notes"]

UNITS = ["", # SEQ
         "", # NAME
         "", # MAIN_ID
         "h", # RA
         "deg", # DEC
         ""]

TYPES = [np.int_, # SEQ
		 np.str_("<U64"), # NAME
         np.str_("<U64"), # MAIN_ID
         coord.angles.core.Angle, # RA
         coord.angles.core.Angle, # DEC
         np.str_("<U512")]

def check_window(table: Table, 
                 config: dict, 
                 i: np.int_ = None) -> np.bool_:
    """
    This function returns False if any target is outside the observation 
    window.
    @params:
        - table: the table of targets to check
        - config: the configuration dictionary
        - i: if not None, check only the index i of the table
    @returns:
        - True if all targets are in the observation window
        - False if at least one target is outside the observation window.
    """
    def check_operation_(line, config):
        """Repeated operation to check all the elements of the table."""
        if line["main_id"] in SPECIAL:
            return True
        if line["ra"] == "":
            return True
        if line["dec"] == "":
            return True
        ra = coord.Angle(line["ra"], unit=u.hourangle)
        dec = coord.Angle(line["dec"], unit=u.degree)
        east = coord.Angle(config["WINDOW_EAST"], unit=u.hourangle)
        west = coord.Angle(config["WINDOW_WEST"], unit=u.hourangle)
        upper = coord.Angle(config["WINDOW_UPPER"], unit=u.degree)
        lower = coord.Angle(config["WINDOW_LOWER"], unit=u.degree)
        if east <= west:
            if ra < east or ra > west:
                print("\033[93m"
                      +("Warning: The target {} "
                        "is outside the observation window.").format(
                         line["name"])
                      +"\033[0m")
                return False
        else:
            if ra < east and ra > west:
                print("\033[93m"
                      +("Warning: The target {} "
                        "is outside the observation window.").format(
                         line["name"])
                      +"\033[0m")
                return False
            if dec < lower or dec > upper:
                print("\033[93m"
                      +("Warning: The target {} "
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
    """
    This function creates the target table.
    @params:
        - cals: list of columns
        - units: list of units
        - type: lust of types
    @returns:
        - a Table object
    """
    return Table(names=cols, dtype=types, units=units)

def select_obj(table: Table, 
               line: Table,
               config: dict) -> Table:
    """
    This function interacts with the user to select targets from one table 
    and append them in another table.
    @params:
        - table: the target table to modify
        - line: a table of targets to choose from
        - config: the configuration dictionary
    @returns:
        - table: the modified table
    """
    table_old = table.copy()
    for i in range(len(line)):
        while line[i]["name"] in line[:i]["name"]: 
            line[i]["name"] = "*"+line[i]["name"]
    line.pprint_all()
    answer = "$"
    print("\033[32m"
          + ("Select the name of the objects "
             "you want to observe (column \"NAME\"):"
             "\n(available commands: [name], all, cancel, done)")
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
        elif answer in HELP:
            print("")
            print_help()
        elif answer.upper() in ALL:
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
               notes: np.str_ = "") -> Table:
    """
    This function allows the user to add manually an entry to the table.
    @params
        - table: the target table
        - config: the configuration dictionary:
        - *params: all the parameters of the target to add
    @returns
        - table: the modified table
    """
    while name in table["name"]: name = "*" + name
    while seq > 0 and seq in table["seq"]: seq += 1

    ra_str = ""
    dec_str = ""
    obj_begin_str = ""
    obj_end_str = ""

    if type(ra) == coord.angles.core.Angle:
        ra_str = ra.to_string(unit=u.hourangle,
                              sep=":",
                              pad=True)
    elif type(ra) == str:
        ra_str = ra

    if type(dec) == coord.angles.core.Angle:
        dec_str = dec.to_string(unit=u.degree,
                                sep=":",
                                pad=True,
                                alwayssign=True)
    elif type(dec) == str:
        dec_str = dec

    table.add_row([seq,
                   name,
                   main_id,
                   ra_str,
                   dec_str,
                   notes])
    check_window(table, config, -1)
    return table

def add_calib(table: Table,
              config: dict = None,
              seq: np.int_ = 0,
              name: np.str_ = "CALIB",
              notes: np.str_ = "") -> Table:
    """
    Add a calibration line to the table.
    @params
        - table: the target table
        - config: the configuration dictionary:
        - *params: all the parameters of the calibration to add
    @returns
        - table: the modified table
    """
    while name in table["name"]:
        name = "*" + name

    while seq > 0 and seq in table["seq"]:
        seq += 1

    table.add_row([seq,
                   name,
                   "CALIBRATION",
                   "",
                   "",
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
               notes: np.str_ = "") -> Table:
    """
    Add a target imported from Simbad using astroquery.Simbad:
    @params
        - table: the target table
        - obj: the table of objects from astroquery.Simbad
        - config: the configuration dictionary:
        - *params: all the parameters of the target to add
    @returns
        - table: the modified table
    """
    for i in range(len(obj)):
        try:
            main_id = obj["main_id"][i]
        except KeyError:
            main_id = obj["MAIN_ID"][i]

        if name == "":
            name = main_id

        while name in table["name"]:
            name = "*" + name

        while seq > 0 and seq in table["seq"]:
            seq += 1
        try:
            ra = coord.Angle(obj["ra"][i], unit=u.degree)
            dec = coord.Angle(obj["dec"][i], unit=u.degree)
        except KeyError:
            ra = coord.Angle(obj["RA"][i], unit=u.hourangle)
            dec = coord.Angle(obj["DEC"][i], unit=u.degree)
        ra_str = ra.to_string(unit=u.hourangle, 
                              sep=":", 
                              pad=True)
        dec_str = dec.to_string(unit=u.degree, 
                                sep=":", 
                                pad=True, 
                                alwayssign=True)
        table.add_row([seq,
                       name,
                       main_id,
                       ra_str,
                       dec_str,
                       notes])
        name = ""
    check_window(table, config)
    return table

def read_cfg(filename: str,
             directory: str = OUT_DIR,
             extension: str = ".cfg") -> dict:
    """
    Read the configuration file
    @params:
        - filename: the configuration file name (without extension)
        - directory: the directory in which the configuration file is located
        - extension: the extension of the configuration file
    @returns:
        - config: the configuration dictionary
    """
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

def set_seq(table: str, 
            config: dict, 
            stack: bool = False):
    """
    Sets the sequence number of all targets in the table, if the sequence 
    number has a duplicate or is zero, excluding special id (calib, ...), and 
    targets with no coordinates.
    @params:
        - table: the target table
        - config: the configuration dictionary
        - stack: if True, the first target is indexed with 1 
        (the last target is indexed with -1, if necessary)
    @returns:
        - table: the updated table
    """
    window_east = coord.Angle(config["WINDOW_EAST"], unit=u.hourangle).degree
    window_west = coord.Angle(config["WINDOW_WEST"], unit=u.hourangle).degree
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
            [coord.Angle(order_table["ra"], unit=u.hourangle).degree], 
            names=["ra_deg"])
    w_before = np.argwhere(order_table["ra_deg"] < window_east)

    order_table["ra_deg"][w_before] += 360
    order_table.sort("ra_deg")

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

def swap(table: Table, 
         name_1: str, 
         name_2: str) -> Table:
    """
    This function swaps to lines in the table
    @params:
        - table: the table of targets
        - name_1: the name of the first target to swap
        - name_2: the name of the second target to swap
    @returns:
        - table: the updated table
    """
    if name_1 not in table["name"] or name_2 not in table["name"]:
        print("\033[93m"
              +"Warning: {} or {} cannot be found in the table.".format(
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

def write_table(table: Table, 
                config: dict, 
                filename=None, 
                extension=".xml", 
                format_="votable", 
                directory=OUT_DIR):
    """
    Writes the table in a file.
    @params:
        - table: the table of targets
        - config: the configuration dictionary
        - filename: the name of the file (if None, YYYY-MM-DD_Location)
        - extension: the file extension (default: .xml)
        - format_: the format of the file (default: votable)
        - directory: the output directory
    @returns:
        - 0 if the code exited with no issues
    """
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
        return 0
    except Exception as e:
        print("\033[91m"
              + ("Error! File cannot be saved "
                 "in the {} directory").format(directory)
              + "\033[0m")
        return 1

def read_table(filename: str,
               table: Table,
               config: dict,
               extension=".xml",
               format_="votable",
               directory=OUT_DIR):
    """
    Read a table previously stored in a file.
    @params:
        - filename: the name of the file
        - table: the (empty) table in which the file is loaded
        - config: the configuration dictionary
        - extension: the file extension (default: .xml)
        - format_: the format of the file (default: votable)
        - directory: the output directory
    @returns:
        - table: the table with the file loaded
    """
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
            if "seq" in line.columns: seq = int(line["seq"])
            else: seq=0
            if "name" in line.columns: name = line["name"]
            else: name = "{}_{}".format(filename, i)
            if "main_id" in line.columns: main_id = line["main_id"]
            else: main_id = name
            if "ra" in line.columns: 
                if ":" not in str(line["ra"]) \
                        and str(line["ra"]).upper() not in ["", " ", "NONE"]:
                    ra = coord.Angle("{}d".format(line["ra"])).to_string(
                            unit=u.hourangle,
                            sep=":",
                            pad=True
                        )
                else: ra = line["ra"]
            else: ra = ""
            if "dec" in line.columns: 
                if ":" not in str(line["dec"]) \
                        and str(line["dec"]).upper() not in ["", " ", "NONE"]:
                    dec = coord.Angle("{}d".format(line["dec"])).to_string(
                            unit=u.degree,
                            sep=":",
                            pad=True,
                            alwayssign=True
                        )
                else: dec = line["dec"]
            else: dec = ""
            if "notes" in line.columns: notes = line["notes"]
            else: notes = ""
            while name in table["name"]:
                name = "*" + name

            for col in line.columns:
                if col not in COLS:
                    notes += " & {}: {}".format(col, line[col])

            table.add_row({"seq": seq,
                           "name": name,
                           "main_id": main_id,
                           "ra": ra,
                           "dec": dec,
                           "notes": notes})

    except Exception as e:
        print("\033[91m"
              + ("Error! File cannot be loaded "
                 "from the {} directory").format(directory)
              + "\033[0m")
        print(e)
        return create_table()
    return table

def make_plot(table: Table, 
              config: dict) -> int:
    """
    Make an observation plot of the targets
    @params:
        - table: the target table
        - config: the configuration dictionary
    @returns: 
        - 0
    """
    def get_timestamp(datetime, ref):
        return (datetime - ref)/np.timedelta64(1, 's')
    def get_time(timestamp, ref):
        return timestamp*np.timedelta64(1, 's') + ref

    if "YII_light_1" in plt.style.available: plt.style.use("YII_light_1")
    plt.rcParams["figure.dpi"] = FIGURE_DPI
    plt.rcParams["figure.facecolor"] = "#121212"
    plt.rcParams['axes.facecolor'] = '#216576'
    COLOR = 'DDDDDD'
    plt.rcParams['text.color'] = COLOR
    plt.rcParams['axes.labelcolor'] = COLOR
    plt.rcParams['xtick.color'] = COLOR
    plt.rcParams['ytick.color'] = COLOR

    fig, ax = plt.subplots(1)
    ax.set_title("Location: {}".format(config["LOCATION"]))

    sun_set = np.datetime64(config["SUN_SET"])
    sun_civil = np.datetime64(config["SUN_SET_CIVIL"])
    sun_nautical = np.datetime64(config["SUN_SET_NAUTICAL"])
    sun_astronomical = np.datetime64(config["SUN_SET_ASTRONOMICAL"])
    sun_rastronomical = np.datetime64(config["SUN_RISE_ASTRONOMICAL"])
    sun_rnautical = np.datetime64(config["SUN_RISE_NAUTICAL"])
    sun_rcivil = np.datetime64(config["SUN_RISE_CIVIL"])
    sun_rise = np.datetime64(config["SUN_RISE"])
    
    obs_begin = np.datetime64(config["OBS_BEGIN"])
    obs_end = np.datetime64(config["OBS_END"])

    timestamp_set = get_timestamp(sun_set, obs_begin)
    timestamp_civil = get_timestamp(sun_civil, obs_begin)
    timestamp_nautical = get_timestamp(sun_nautical, obs_begin)
    timestamp_astronomical = get_timestamp(sun_astronomical, obs_begin)
    timestamp_rastronomical = get_timestamp(sun_rastronomical, obs_begin)
    timestamp_rnautical = get_timestamp(sun_rnautical, obs_begin)
    timestamp_rcivil = get_timestamp(sun_rcivil, obs_begin)
    timestamp_rise = get_timestamp(sun_rise, obs_begin)

    timestamp_begin = get_timestamp(obs_begin, obs_begin)
    timestamp_end = get_timestamp(obs_end, obs_begin)

    st_begin = coord.Angle(config["ST_BEGIN"], unit=u.hourangle).degree
    st_end = coord.Angle(config["ST_END"], unit=u.hourangle).degree
    window_east = coord.Angle(config["WINDOW_EAST"], unit=u.hourangle).degree
    window_west = coord.Angle(config["WINDOW_WEST"], unit=u.hourangle).degree

    if st_end < st_begin: 
        st_end += 360

    slope = (timestamp_end - timestamp_begin)/(st_end - st_begin)
    
    def compute_timestamp(st, a=slope, sb=st_begin, hb=timestamp_begin):
        return (st-sb) * a + hb
    def compute_st(timestamp, a=slope, sb=st_begin, hb=timestamp_begin):
        return (timestamp-hb) / a + sb

    N = len(table)

    # Night times
    ax.axvspan(sun_set, sun_rise, color="k", alpha=0.2)
    ax.axvspan(sun_civil, sun_rcivil, color="k", alpha=0.4)
    ax.axvspan(sun_nautical, sun_rnautical, color="k", alpha=0.6)
    ax.axvspan(sun_astronomical, sun_rastronomical, color="k", alpha=0.8)

    # Observation time
    ax.axvspan(obs_begin, obs_end, color="C0", alpha=0.3)

    # Observation range
    if st_begin < window_east:
        st_begin += 360
    if window_west < st_end:
        window_west += 360
    a_before = st_begin - window_east
    a_after =  window_west - st_end
    
    plot_colors = ["#ED1C24", "#E8BD0F"]
    N_colors = len(plot_colors)
    N_same = 4

    # Targets
    for i in range(N):
        ra = coord.Angle(table["ra"][i], unit=u.hourangle).degree
        if ra > window_west:
            ra -= 360
        if ra < window_east:
            ra += 360
        ra_before = ra + a_before
        ra_after = ra - a_after

        h = compute_timestamp(ra)
        h_before = compute_timestamp(ra_before)
        h_after = compute_timestamp(ra_after)

        time_ra = get_time(h, obs_begin)
        time_before = get_time(h_before, obs_begin)
        time_after = get_time(h_after, obs_begin)
        
        color_index = (i // N_same) % N_colors
        color = plot_colors[color_index]
        
        ax.plot([time_before, time_after], [i,i], color=color, marker="|")
        ax.scatter([time_ra], [i], s=5, color=color)

        ax.text(time_after, i, table["name"][i] + "  ", 
                horizontalalignment="right", verticalalignment="center")
        ra_text = table["ra"][i]
        dec_text = table["dec"][i]
        ax.text(time_before, i, "  {}{}".format(ra_text, dec_text),
                horizontalalignment="left", verticalalignment="center")

    ax.set_xlabel("Observation date and time (UTC)")
    ax.set_yticks(range(N), table["seq"])

    plt.show(block=True)

    # obs_time = np.arange(sun_set, sun_rise, dtype="datetime64[m]")
    return 0

def print_help():
    hilfe = ("\033[36m"
             "Help will be always given to those who ask for it [1].\n"
             "Available commands (not case sensitive):\n"
             "\t- help, h, ?: show this page\n"
             "\t- quit, exit, q: quit the current code "
             "(WARNING: this does not save the current state!)\n"
             "\t- write, save: write the current table in a file "
             "(no options available yet)\n"
             "\t- read [filename], open [filename], load [filename]: "
             "loads the file \"filename\" "
             "in the current table (no additional options available yet)\n"
             "\t- calibration, calib: adds a calibration in the target list\n"
             "\t- simbad [object name], object [object name]: "
             "add an object from simbad\n"
             "\t- search [ra] [dec] [radius], region [ra] [dec] [radius]: "
             "search a region centred on the ra/dec coordinates, "
             "with a given radius (ra is given in hour, dec in degree "
             "and the radius in any specified unit, "
             "e.g. search 01:03:40 +35:40:20 30\')\n"
             "\t- manual [name] -s [seq (optional)], "
             "add [name] -s [seq (optional)]: "
             "manually add a target (only the name the sequence are "
             "available for now)\n"
             "\t- sequence, seq: computes the sequence order for each "
             "target\n"
             "\t- check: check if all targets are in the observation field\n"
             "\t- plot, graph: creates a graph with all the targets "
             "observation date and time"
             "\n"
             "General actions: \n"
             "\t- cancel, back: cancels the current action\n"
             "\t- yes, y, 1: yes\n"
             "\t- no, n, 0, or anything else: no\n"
             "\t- all, *: select all\n"
             "\t- done, ok: confirm, save the current state and quit\n"
             "\n"
             "[1] Pr. Dumbledore, Albus Percival Wulfric Brian, 1992, Hogwarts School of Witchcraft and Wizardry"
             "\033[0m")
    print(hilfe)
    return None

def resolve_input(text: str, 
                  table: Table, 
                  config: dict) -> Table:
    """
    Resolve the input to execute the expected function in an 
    interactive way
    @params:
        - text: the input string entered by the user
        - table: the table of targets
        - config: the configuration dictionary
    @returns:
        - table: the updated table
    """
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
            ra = coord.Angle(args[1], unit=u.hourangle)
            dec = coord.Angle(args[2], unit=u.degree)
            region = coord.SkyCoord(ra=ra, dec=dec)
            radius = coord.Angle(args[3])
            obj = Simbad.query_region(region, radius)
            swap_table = add_simbad(swap_table, obj, config)
            select_obj(table, swap_table, config)
        elif args[0].upper() in MANUAL: # TODO ajouter args
            if "-s" in args:
                name = args[1]
                i = 2
                while args[i] != "-s":
                    name += (" " + args[i]) 
                    i+=1
                seq = int(args[i+1])
            else:
                name = " ".join(args[1:])
                seq = 0
            add_manual(table, config, seq, name)
        elif args[0].upper() in CALIB: # TODO ajouter args
            add_calib(table, config)
        elif args[0].upper() in SEQ: # FIXME marche bien une fois mais pas deux ??
            set_seq(table, config)
        elif args[0].upper() in CHECK:
            check_window(table, config) 
        elif args[0].upper() in PLOT:
            make_plot(table, config)
        elif args[0].upper() in HELP:
            print_help()
        elif args[0].upper() in QUIT:
            None
        else:
            print("\033[93m"
                  + ("Warning: {} is not recognized "
                     "as a valid keyword.").format(args[0])
                  + "\033[0m")
    except Exception as e:
        print("\033[91m"
              + "Error: Your operation is not recognized".format(args[0])
              + "\033[0m")
        print(e)
    return table



