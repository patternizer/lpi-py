![image](https://github.com/patternizer/lpi-py/blob/main/8-out-of-10-bat-stripes.png)
![image](https://github.com/patternizer/lpi-py/blob/main/nbis-2023-bat-10-biodiversity-stripes.png)
![image](https://github.com/patternizer/lpi-py/blob/main/nbis-2023-bat-10-gam.png)
![image](https://github.com/patternizer/lpi-py/blob/main/lpd-2020-biodiversity-class-temporal-distribution.png)
![image](https://github.com/patternizer/lpi-py/blob/main/linear-interpolation-gap-filling.png)

# lpi-py

Python code translation from R to re-calculate the Living Planet Index for the LPD 2020 database and merge in local data from the Norfolk Biodiversity Information Service (NBIS).

## Contents

* `gap-finder.py` - python algorithm to find and infill gaps using a combination of linear regression interpolation and a generalised additive model
* `plot-species-lpd-2020.py` - python code to ingest biodiversity abundance counts 1950-2020 from the Living Planet Index (LPI) dataset (LPD2020) and plot species total abundnace timeseries and source locations on a map
* `plot-species-nbis-2023.py` - python code to ingest biodiversity abundance counts 1950-2023 from the Norfolk Biodiversity Information Service (NBIS) 2023 dataset and plot species total abundnace timeseries and source locations on a map
* `nbis-2023-lpi.py` - python code to ingest biodiversity abundance counts 1950-2023 from the Norfolk Biodiversity Information Service (NBIS) 2023 dataset and compute the LPI for selected species using a new implementation of the generalised additive model.
* `8-out-of-10-bat-stripes.png` - NBIS 2023 LPI for 8 out of 10 Norfolk bat species tribute to Chris Packham and team's amazing "8 Out of 10 Bats"

The first step is to clone the latest lpi-py code and step into the check out directory: 

    $ git clone https://github.com/patternizer/lpi-py.git
    $ cd lpi-py

### Usage

The code was tested locally in a Python 3.8.11 virtual environment. 
Run with:

    $ python nbis-2023-lpi.py
    $ python plot-species-nbis-2023.py (research code)
    $ python plot-species-lpd-2020.py (research code)
    $ python gap-finder.py (optional)
    
## License

The code is distributed under terms and conditions of the [MIT license](https://opensource.org/licenses/MIT).

## Contact information

* [Michael Taylor](patternizer@proton.me)






