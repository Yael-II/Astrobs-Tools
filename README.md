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

...

### EQUATOR

First, select a configuration file (if none, default is selected).
Then, select the action you wish to perform. 

Available commands (not case sensitive):
- help, h, ?: show this page
- quit, exit, q: quit the current code
- cancel, back: cancels the current action
- write, save: write the current table in a file 
(no options available yet)
- read, open [filename]: loads the file "filename" 
in the current table (no additional options available yet)
- calibration, calib: adds a calibration in the target list
- simbad [object name], object [object name]: 
add an object from simbad
- search [ra] [dec] [radius], region [ra] [dec] [radius]: 
search a region centred on the ra/dec coordinates, 
with a given radius (coordinates should be expressed as 
12h30m30s, 90d30m30s or 90.555d)
- sidereal, st: computes the sidereal time for each target
- sequence, seq: computes the sequence order for each 
target
- check: check if all targets are in the observation field

