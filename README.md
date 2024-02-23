# Mammals under pressure: Presence data for assessing extinction of endemic, threatened and mammals subject to use in Colombia

This is the code used to obtain the graphics, statistics, and tables provided in the paper. The results were obtained by applying the dataset, which includes 97,944 records corresponding to 136 species, of which 38 are endemic, 92 identified as species subject to use by humans in the literature, and 56 categorized as either Data Deficient or threatened.

## Contents

The dataset contains

- 00-rawData: The original Data used in the publication
- 01-scripts Contains the file "Script_connectivity_andean_bear.R" with the steps to create the resources for the manuscript
- 02-outData: It contains the data obtained after applying the analysis [it is developing]
- 03 Figures: It has the figures used in the manuscript 
- 04-Manuscript: Contains the R Markdown with the manuscript in the "manuscript.Rmd" file. Also, have the bibliography in "mybibifile.bib", and the Cumulative Layout Shift in "elsarticle.cls" file.
- renv: It has the local storage created to reproducibility the data. See instructions after.

## Installing Required R Packages

To reproduce the information described in the publication you need to the "00RawData" folder with the dataset  descrited in the publication, and the file "DPMC_Analysis.R" in the "01Script" folder

To reproduce the analysis and visualizations in this repository, you will need to install several R packages. You can install these packages using the install.packages() function in R. Here are the packages you need to install:

terra: Install the terra package, which provides advanced capabilities for spatial data manipulation and analysis.
''' 
install.packages("terra")
'''

tmap: Install the tmap package, which offers powerful mapping functionalities for spatial data visualization.

'''
install.packages("tmap")
'''

tidyverse: Install the tidyverse package, which includes a collection of R packages for data manipulation, visualization, and analysis.
'''
install.packages("tidyverse")
'''

sp: Install the sp package, which provides classes and methods for spatial data representation and manipulation.
'''
install.packages("sp")
'''

geodata: Install the geodata package, which offers functions for working with geographical data.
'''
install.packages("geodata")
'''

treemap: Install the treemap package, which allows you to create treemaps for hierarchical data visualization.
'''
install.packages("treemap")
'''

forcats: Install the forcats package, which provides tools for working with categorical variables (factors) in R.
'''
install.packages("forcats")
'''

Once you have installed these packages, you can load them into your R session using the library() function:

'''
library(terra)
library(tmap)
library(tidyverse)
library(sp)
library(geodata)
library(treemap)
library(forcats)
'''

Now you're ready to run the code and explore the analyses and visualizations provided in this repository. If you encounter any issues during installation or usage, please refer to the package documentation or consult the R community for assistance.

## Repository Structure

The following tree represents the files that are organized in this repository

├── 00RawData
├── 01Scripts
│   └── DPMC_Analysis.R
├── 02Manuscript
│   └── Datapaper endemic, threatened and mammals subject to use final version.pdf
├── 03ProcessedData
├── 04Plots
│   ├── recods_BoR.jpeg
│   ├── recods_per_Family.jpeg
│   └── recods_per_order.jpeg
├── 05OutData
│   ├── Table2.csv
│   └── Tables.xls
├── LICENSE
└── README.md
