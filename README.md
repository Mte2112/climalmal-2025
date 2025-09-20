## Project Description

This repo contains the scripts needed to reproduce the results in Elling et al., 2025 *Tropical oceans drive Malawi’s malaria risk*   
Code formatted with [Black](https://pypi.org/project/black/)

## Directory structure

```
├── README.md                           
├── analysis_main.py                    # produces the main results (all but CMIP6 analysis) 
├── analysis_cmip6.ipynb                # produces CMIP6 soil moisture results
├── load_preprocessed_data.py           # data loading utilities 
├── analysis_tools.py                   # helper functions for main analysis
├── format_data.py                      # data processing/formatting formatting functions 
├── figures.ipynb                       # jupyter notebook for figure creation 
├── data/
│   └── shape/                         # shapefiles for Malawi boundaries
│── preprocessed/                      # preprocessed data from analysis.py
└── figures/                           
    └── manuscript/                    # output figures, manuscript
    └── supp/                          # output figures, supplementary
```

## How to run  
Clone the repository  
```bash
git clone https://github.com/Mte2112/climalmal-2025  
```  

### 1. Data setup
- Unpack shapefile tarball  
```bash
gunzip -c data/shapefiles.tar.gz | tar -xf - -C data
```  
- Optional: To dramatically speed up the process, [download preprocessed inputs](https://zenodo.org/records/17166386) from Zenodo. If using preprocessed inputs, skip step 3-Main analysis (data processing script)  
- To run the full analysis, download public climate datasets. See step 3. 

### 2. Set up the environment
```bash
pip install -r requirements.txt  
conda install -c conda-forge esmpy
```

### 3. Run the analysis

### Main analysis   
#### Public data inputs (if running main analysis)
- [Malaria incidence](https://zenodo.org/records/17161438) (district-level, monthly)
- [SST](https://climatedataguide.ucar.edu/climate-data/sst-data-noaa-optimal-interpolation-oi-sst-analysis-version-2-oisstv2-1x1) (NOAA OI SST)
- [ERA5 reanalysis](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview) (300/850 mb wind, 300/850 mb geopotential, 2m temperature, vimfd)
- [CHIRPS precipitation](https://www.chc.ucsb.edu/data/chirps)    
- [GRACE-DA-DM soil moisture](https://disc.gsfc.nasa.gov/datasets/GRACEDADM_CLSM025GL_7D_3.0/summary?keywords=grace%20soil%20moisture) 
- [Elevation](https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/tiled/) (ETOPO1)
- Malawi national and district shapefiles (in repo)  

Update the file paths at the top of analysis_main.py to reflect data locations  

Run the analysis
```bash
python analysis_main.py
```
*Inludes all aspects of core analysis. Saves preprocessed output for easy import & plotting in figures.ipynb. Does not include CMIP6 analysis from Figure 7 in manuscript.   

### CMIP6 projections 
`analysis_cmip6.ipynb`  

CMIP6 analysis via Pangeo for efficient, reproducible analysis. Completely different pipeline from the main analysis, so kept separate.  

### 4. Create figures seen in manuscript
`figures.ipynb` 


## Citation

Elling, M. T., Karnauskas, K. B., Kowalcyk, M., Mategula, D., Chirombo, J., Livneh, B., McCann, R., & Buchwald, A. G. (2025). Tropical oceans drive Malawi’s malaria risk. GitHub. https://github.com/Mte2112/climalmal-2025

