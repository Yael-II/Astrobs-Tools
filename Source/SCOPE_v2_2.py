"""
* SCOPE: Sky Coordinates for Observations Python Estimator
* Version 2.1 - 2024-09-21
@ Yaël Moussouni
@ Observatory of Strasbourg 
"""
import astropy.time as time
import astropy.coordinates as coord
import astropy.units as u

YES = ["YES", "Y", "1"]
OUT_DIR = "./Output/"
# Definitions
def get_sun(sun_time, sun_frame, angle, reverse=False):
    i = 0
    while coord.get_sun(sun_time).transform_to(sun_frame).alt.value > angle:
        if reverse: 
            sun_time -= 1 * u.min
        else: 
            sun_time += 1 * u.min
        i += 1
        if i > 60*24:
            print("\033[93m"
                  + "Warning: Sun did not reached the expected angle."
                  + "\033[0m")
            return time.Time("1970-01-01T12:00") # Error date
    return sun_time

# Locations
available_obs = ["Observatoire astronomique de Strasbourg (ObAS)", 
                 "Observatoire de Haute-Provence (OHP)"]
available_obs_short = ["ObAS",
                       "OHP"]
available_obs_lat = [coord.Latitude("""48d34m59s"""), 
                     coord.Latitude("""43d55m53s""")]  
available_obs_lon = [coord.Longitude("""07d46m07s"""), 
                     coord.Longitude("""5d42m45s""")]


def main():
    try:
        print("\033[32m"+ "Please select a location:"+"\033[0m")
        for i in range(len(available_obs)):
            print("\033[32m"+ "\t{}. ".format(i+1)+"\033[0m"+available_obs[i]+".")
        selection = int(input("\033[32m"+ "Choice: "+"\033[0m"))-1
        if selection >= len(available_obs):
            selection = -1
            print("\033[93m"
                  + ("Warning: Location unavailable, "
                     "selected the last one instead.")
                  + "\033[0m")
    except:
        selection = -1
        print("\033[91m"
              + ("Error! Input value not recognized, "
                 "selected the last location instead")
              + "\033[0m")

    obs_name = available_obs[selection]
    obs_short = available_obs_short[selection]
    obs_lat = available_obs_lat[selection]
    obs_lon = available_obs_lon[selection]
    obs_loc = coord.EarthLocation(lon=obs_lon, lat=obs_lat)

    # Date and time
    try:
        obs_date = input("\033[32m"
                         + "Date of observation [YYYY-MM-DD]: "
                         + "\033[0m")
        test_date = time.Time(obs_date)
    except:
        obs_date = time.Time("1970-01-01", format="iso").now()
        obs_date = obs_date.to_string()[0:10]
        print("\033[91m"
              + "Error! Input value not recognized, selected today instead"
              + "\033[0m")
    else:
        del test_date
    print("\033[34m" + "Computing Sun position..." + "\033[0m")
    sun_frame = coord.AltAz(location=obs_loc)
    sun_time = time.Time(obs_date + " 12:00:00", scale="utc")

    sun_set = get_sun(sun_time, sun_frame, 0)
    sun_civil = get_sun(sun_set, sun_frame, -6)
    sun_nautical = get_sun(sun_civil, sun_frame, -12)
    sun_astronomical = get_sun(sun_nautical, sun_frame, -18)

    sun_time = time.Time(obs_date + " 12:00:00", scale="utc") + 1 * u.day

    sun_rise = get_sun(sun_time, sun_frame, 0, reverse=True)
    sun_rcivil = get_sun(sun_rise, sun_frame, -6, reverse=True)
    sun_rnautical = get_sun(sun_rcivil, sun_frame, -12, reverse=True)
    sun_rastronomical = get_sun(sun_rnautical, sun_frame, -18, reverse=True)
    print("\033[34m" + "done!" + "\033[0m")

    print("\033[36m"
          + "Night time [UTC]:"
          + "\033[0m")
    print("\033[36m"
          + "\t    === Evening ===    "
          + "\033[0m")
    print("\033[36m"
          + "\t               sunset: "
          + "\033[0m"
          +sun_set.to_string())
    print("\033[36m"
          + "\t       civil twilight: "
          + "\033[0m"
          +sun_civil.to_string())
    print("\033[36m"
          + "\t    nautical twilight: "
          + "\033[0m"
          +sun_nautical.to_string())
    print("\033[36m"
          + "\tastronomical twilight: "
          + "\033[0m"
          +sun_astronomical.to_string())
    print("\033[36m"
          + "\t    === Morning ===    "
          + "\033[0m")
    print("\033[36m"
          + "\tastronomical twilight: "
          + "\033[0m"
          +sun_rastronomical.to_string())
    print("\033[36m"
          + "\t    nautical twilight: "
          + "\033[0m"
          +sun_rnautical.to_string())
    print("\033[36m"
          + "\t       civil twilight: "
          + "\033[0m"
          +sun_rcivil.to_string())
    print("\033[36m"
          + "\t              sunrise: "
          + "\033[0m"
          +sun_rise.to_string())


    obs_begin = input("\033[32m"
                      + "Begin time of observation [hh:mm format; UTC]: "
                      + "\033[0m")
    obs_end =   input("\033[32m"
                      + "  End time of observation [hh:mm format; UTC]: "
                      + "\033[0m")

    if obs_begin == "set" or obs_begin == "rise":
        obs_begin_time = sun_set
    elif obs_begin == "civil":
        obs_begin_time = sun_civil
    elif obs_begin == "nautical":
        obs_begin_time = sun_nautical
    elif obs_begin == "astronomical" or obs_begin == "auto" or obs_begin == "":
        obs_begin_time = sun_astronomical
    else: 
        try:
            obs_begin_time = time.Time(obs_date + " " + obs_begin, scale="utc")
        except:
            obs_begin_time = sun_astronomical
            print("\033[91m"
                  + ("Error! Input value not recognized, "
                     "selected astronomical twilight instead.")
                  + "\033[0m")


    if obs_end == "set" or obs_end == "rise":
        obs_end_time = sun_rise
    elif obs_end == "civil":
        obs_end_time = sun_rcivil
    elif obs_end == "nautical":
        obs_end_time = sun_rnautical
    elif obs_end == "astronomical" or obs_end == "auto" or obs_end == "":
        obs_end_time = sun_rastronomical
    else:
        try:
            obs_end_time = time.Time(obs_date + " " + obs_end, scale="utc")
        except:
            obs_end_time = sun_rastronomical
            print("\033[91m"
                  + ("Error! Input value not recognized, "
                     "selected astronomical twilight instead.")
                  + "\033[0m")

    if obs_end_time.value < obs_begin_time.value:
        obs_end_time += 1*u.day

    obs_duration = obs_end_time - obs_begin_time

    # Coordinates

    print("\033[36m"
          + "Observation constraints"
          + "\033[0m")
    try:
        min_alt = coord.Angle(
            input("\033[32m"
                  + "           Minimum altitude above horizon [deg]: "
                  + "\033[0m")
            + "d")
    except:
        min_alt = coord.Angle(50 * u.deg)
        print("\033[91m"
              + "Error! Input value not recognized, selected 50° instead."
              + "\033[0m")
    try:
        min_dec = coord.Angle(
            input("\033[32m"
                  + "     Minimum declination around the north [deg]: "
                  + "\033[0m")
            + "d")
    except:
        min_dec = coord.Angle(30 * u.deg)
        print("\033[91m"
              + "Error! Input value not recognized, selected 30° instead."
              + "\033[0m")
    try:
        window_east = coord.Angle(
            input("\033[32m"
                  + "Observation window before crossing meridian [h]: "
                  + "\033[0m")
            + " hours")
    except:
        window_east = coord.Angle(4 * u.h)
        print("\033[91m"
              + "Error! Input value not recognized, selected 4h instead."
              + "\033[0m")
    try:
        window_west = coord.Angle(
            input("\033[32m"
                  + " Observation window after crossing meridian [h]: "
                  + "\033[0m")
            + " hours")
    except:
        window_west = coord.Angle(2 * u.h)
        print("\033[91m"
              + "Error! Input value not recognized, selected 2h instead."
              + "\033[0m")

    obs_sky_rotation = coord.Angle(
            "{} hours".format(obs_duration.to_value("h")))\
            * (1*u.sday/u.day).decompose()
       
    # * Date/time/location objects
    obs = time.Time(obs_begin_time, location=obs_loc)

    # * Local sidereal time
    st = obs.sidereal_time("apparent") # ? apparent or absolute

    # * Borders and corners coordinates
    dec_upper = (90*u.deg-min_dec).wrap_at(360*u.deg)
    dec_lower = (min_alt + obs_loc.lat - 90*u.deg).wrap_at(180*u.deg)
    a_west = (st + obs_sky_rotation + window_west).wrap_at(24*u.h)
    a_east = (st - window_east).wrap_at(24*u.h)

    if a_east.degree <= a_west.degree:
        comp_symb = "&&"
        comp_text = "AND"
    else: 
        comp_symb = "||"
        comp_text = "OR"
    # * Output
    print("")
    print("\033[36m"
          + "Location: "
          + "\033[0m"
          + obs_name)
    print("\033[36m"
          + "\tlat: "
          + "\033[0m" 
          + obs_lat.to_string())
    print("\033[36m"
          + "\tlon: "
          + "\033[0m" 
          + obs_lon.to_string())
    print("\033[36m"
          + "Observation"
          + "\033[0m")
    print("\033[36m"
          + "\tbegin date and time: "
          + "\033[0m" 
          + obs_begin_time.to_string())
    print("\033[36m"
          + "\t  end date and time: "
          + "\033[0m" 
          + obs_end_time.to_string())
    print("\033[36m"
          + "\tbegin sidereal time: "
          + "\033[0m" 
          + st.to_string(unit=u.hour))
    print("\033[36m"
          + "\t  end sidereal time: "
          + "\033[0m" 
          + (st + obs_sky_rotation).wrap_at(24*u.h).to_string(unit=u.hour))
    print("\033[36m"
          + "\t       sky rotation: "
          + "\033[0m" 
          + obs_sky_rotation.to_string(unit=u.hour))
    print("\033[36m"
          + "Sky coordinates window"
          + "\033[0m")
    print("\033[36m"
          + "\t  east ra: "
          + "\033[0m" 
          + a_east.to_string(unit=u.hour))
    print("\033[36m"
          + "\t  west ra: "
          + "\033[0m" 
          + a_west.to_string(unit=u.hour))
    print("\033[36m"
          + "\tupper dec: "
          + "\033[0m" 
          + dec_upper.to_string(unit=u.degree))
    print("\033[36m"
          + "\tlower dec: "
          + "\033[0m" 
          + dec_lower.to_string(unit=u.degree))
    print("\033[36m"
          + "Simbad query constraints: "
          + "\033[0m"
          + "WHERE (ra < {} {} ra > {}) AND (dec < {} AND dec > {}))".format(
              a_west.degree,
              comp_text,
              a_east.degree,
              dec_upper.degree,
              dec_lower.degree))
    print("\033[36m"
          + "Vizier query constraints:"
          + "\033[0m")
    print("\033[36m"
          + "\tRA: "
          + "\033[0m" 
          + "< " 
          + a_west.to_string(unit=u.hour, sep=":") 
          + " {} ".format(comp_symb) 
          + "> " 
          + a_east.to_string(unit=u.hour, sep=":"))
    print("\033[36m"
          + "\tDE: "
          + "\033[0m" 
          + "< " 
          + dec_upper.to_string(unit=u.degree, sep=":") 
          + " && " 
          + "> " 
          + dec_lower.to_string(unit=u.degree, sep=":"))

    print("")
    print("\033[32m" 
          + "Write output to file ? [yes/no]" 
          + "\033[0m")
    write_output = input("\033[32m"
                         + "Choice: "
                         + "\033[0m").upper() in YES
    if write_output:
        try:
            with open("{}{}_{}.cfg".format(OUT_DIR, obs_date, obs_short),
                      "w+") as file:
                file.write("# File generated with SCOPE_v2\n")
                file.write("LOCATION: ")
                file.write(obs_short + "\n")
                file.write("LAT: ")
                file.write(obs_lat.to_string() + "\n")
                file.write("LON: ")
                file.write(obs_lon.to_string() + "\n")
                file.write("SUN_SET: ")
                file.write(sun_set.to_string() + "\n")
                file.write("SUN_SET_CIVIL: ")
                file.write(sun_civil.to_string() + "\n")
                file.write("SUN_SET_NAUTICAL: ")
                file.write(sun_nautical.to_string() + "\n")
                file.write("SUN_SET_ASTRONOMICAL: ")
                file.write(sun_astronomical.to_string() + "\n")
                file.write("SUN_RISE_ASTRONOMICAL: ")
                file.write(sun_rastronomical.to_string() + "\n")
                file.write("SUN_RISE_NAUTICAL: ")
                file.write(sun_rnautical.to_string() + "\n")
                file.write("SUN_RISE_CIVIL: ")
                file.write(sun_rcivil.to_string() + "\n")
                file.write("SUN_RISE: ")
                file.write(sun_rise.to_string() + "\n")
                file.write("OBS_BEGIN: ")
                file.write(obs_begin_time.to_string() + "\n")
                file.write("OBS_END: ")
                file.write(obs_end_time.to_string() + "\n")
                file.write("ST_BEGIN: ")
                file.write(st.to_string(unit=u.hour) + "\n")
                file.write("ST_END: ")
                file.write((st + obs_sky_rotation).wrap_at(
                    24*u.h).to_string(unit=u.hour) + "\n")
                file.write("WINDOW_EAST: ")
                file.write(a_east.to_string(unit=u.hour) + "\n")
                file.write("WINDOW_WEST: ")
                file.write(a_west.to_string(unit=u.hour) + "\n")
                file.write("WINDOW_UPPER: ")
                file.write(dec_upper.to_string(unit=u.degree) + "\n")
                file.write("WINDOW_LOWER: ")
                file.write(dec_lower.to_string(unit=u.degree) + "\n")
            print("\033[34m" 
                  + "File saved in the {} directory ".format(OUT_DIR)
                  + "with the name {}_{}.cfg".format(obs_date, obs_short)
                  + "\033[0m")
            print("\033[34m" + "done!" + "\033[0m")
        except:
            print("\033[91m"
                  + "Error! File cannot be save in the output directory."
                  + "\033[0m")
    return 0
if __name__ == "__main__":
    main()
