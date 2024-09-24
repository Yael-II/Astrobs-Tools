"""
* EQUATOR: Equator Queries simbAd to create Tables of Objects
* Version 1 - 2024-09-21
@ Yaël Moussouni
@ Observatory of Strasbourg 
"""
import astropy.time as time
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table
from astroquery.simbad import Simbad
import astrobs_v1_toolbox as toolbox

YES = ["YES", "Y", "1"]

config = toolbox.read_cfg("2024-09-25_ObAS.cfg")

objects = toolbox.create_table()
name = "NGC 2496"
obj = Simbad.query_object(name)
toolbox.add_simbad(objects, obj, config, name=name)
name = "NGC 792"
obj = Simbad.query_object(name)
toolbox.add_simbad(objects, obj, config, name=name)
toolbox.add_calib(objects)
toolbox.add_calib(objects)
toolbox.add_manual(objects, config, name="test", main_id="IGNORE", ra=coord.Angle("12h30m30s"))

toolbox.set_st_window(objects, config)

print(objects)
