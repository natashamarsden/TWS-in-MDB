# TWS-in-MDB
This repository contains code used for an Engineering Capstone project at the University of Melbourne, entitled "Assessment of Hydrological Model Performance using GRACE-derived Terrestrial Water Storage in the Murray-Darling Basin." The authors are Natasha Marsden, Arfaa Feezanul Islam, and Jun Liu.

The notebooks have been saved with all outputs (e.g., figures, maps, tables) already rendered, so they can be viewed directly without requiring re-execution. Note that the raw datasets are not included in the repository due to file size limitations. However, all datasets are publicly available and can be accessed via the sources detailed in Section 3.2 of the report. Minimal preprocessing was performed outside of the Python environment.

- functions.py: Python script with all necessary functions
- MDB_Regions.ipynb: Sub-basin delineation
- P_and_ET_Data.ipynb: Loading and processing of all P and ET datasets (except MODIS)
- MODIS_Pre-Processing.ipynb: Loading and processing of MODIS data (this had to be done separately given the complexity of the dataset)
- Runoff_Data_Pre-Processing.ipynb: Gap filling the missing runoff data
- Runoff_Data.ipynb: Processing the runoff data
- Data_Processing.ipynb: Loading and processing (including downscaling and spatial averaging) of GRACE; Final processing of all other datasets
- Main_Code_File.ipynb: Main code file which was used to generate all graphs and outputs
