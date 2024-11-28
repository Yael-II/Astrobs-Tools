# Astronomical Observation Tools

Astronomical Observation Tools (Astrobs Tools) is composed of:
- SCOPE: a code to get the sky coordinate limits from observation constraints;
- EQUATOR: a software to organize targets before observations.


## Requirements

Astrobs Tools requires `Python 3.10` or newer

## Setup

First step is to initialize a virtual environment for the software with:
```
python3 -m venv ./venv
```

Then install the required packages with:
```
source activate.sh && pip install -r requirements.txt && deactivate
```

You may also require to authorize the execution of Astrobs Tools with:
```
chmod u+x Astrobs-Tools.sh 
```

## Usage

To launch Astrobs tools, simply execute:
```
./Astrobs-Tools.sh
```
Then select the desired code

### SCOPE

1. Select a location (available location: ObAS, OHP) — default: (last one)
2. Date of observation `[YYYY-MM-DD]` - default: today
3. Begin time of observation `[hh:mm format, UTC]` — default: astronomical twilight
4. End time of observation `[hh:mm format, UTC]` — default:  astronomical twilight
5. Minimum altitude above horizon `[deg]` — default: 50°
6. Minimum declination around the north `[deg]` — default: 30°
7. Observation window before crossing the meridian `[h]` — default: 4h
8. Observation window after crossing the meridian `[h]` - default: 2h

Observation parameters are then printed

- Write output to file `[yes/no]` — default: no
By default, the file is saved in the Output directory, with the name `YYYY-MM-DD_Location.cfg`

### EQUATOR

First, select a configuration file (if none, default is selected).
Then, select the action you wish to perform. 

Available commands (not case sensitive):
- `help`, `h`, `?`: show this page
- `quit`, `exit`, `q`: quit the current code 
(WARNING: this does not save the current state!)
- `write`, `save`: write the current table in a file 
(no options available yet)
- `read [filename]`, `open [filename]`: loads the file "filename" 
in the current table (no additional options available yet)
- `calibration`, `calib`: adds a calibration in the target list
- `simbad [object name]`, `object [object name]`: 
add an object from simbad
- `search [ra] [dec] [radius]`, `region [ra] [dec] [radius]`: 
search a region centred on the ra/dec coordinates, 
with a given radius (coordinates should be expressed as 
`12h30m30s`, `90d30m30s` or `90.555d`)
- `manual [name] -s [seq]`, `add [name] -s [seq]`: manually add a target (only the name and the sequence are available for now)
- `sidereal`, `st`: computes the sidereal time for each target
- `sequence`, `seq`: computes the sequence order for each 
target
- `check`: check if all targets are in the observation field

General actions:
- `cancel`, `back`: cancels the current action
- `yes`, `y`, `1`: yes
- `no`, `n`, `0`, or anything else: no
- `all`, `*`: select all
- `done`, `ok`: confirm, save the current state and quit
