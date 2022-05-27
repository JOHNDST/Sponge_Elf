# Sponge_Elf
A LIDs allocation decision-aid program
![001_final-cs6-bd-01](https://user-images.githubusercontent.com/62609466/170732784-88fa9edb-25ff-4509-bba0-c17e65a2236c.jpg)
# Description
Sponge Elf is an integrated decision-aid program developed by us to guide the spatial allocation of LIDs. Combining hydrological models and a multi-objective optimization algorithm, Sponge Elf enables automatic sub-catchment modeling and gives direct spatial suggestions for the most cost-effective LID allocation. 

# How to Use?
## 1. Install dependent libraries
This algorithm is written in Python language and equipped with SWMM5 (https://pypi.org/project/SWMM5/) and Pymoo (https://pymoo.org/), please complete the installation of the dependent libraries first. Commonly used libraries including Matplotlib, Pandas, GDAL, bisect, and Numpy are also needed. EPA SWMM 5.1.012 (Only this version is currently supported) (https://www.epa.gov/water-research/storm-water-management-model-swmm) is recommended to use when you want to localize you own model.

## 2. Pre-treatment 
This section of code is used to construct sub-catchments model for calculation based on given data. Sub-catchments will be generated and added into the INP template for simulation. Currently the program can't modify the pre-set hydrological data in the template, if you want to change something, you can change it mannually by using EPA SWMM. 

## 3. Set decision variables and objective functions
We set 4 types of LIDs for the algorithm to allocate, you can change thier occupying propotion of each sub-catchment cells or change their costs in our program. Then, you can set your hydrological indicators for optimization. 

## 4. Multi-objective Spatial Optimization
After setting some basic parameters you can start the multi-objective spatial optimization, the NSGA-II algorithm will generate LIDs layouts by putting LID types in each sub-catchment cell, and evaluate their hydrological indicators and costs. Finally a pareto solutions set will be obtained.

## 5. Visualization and Analysis
In this part you are able to see your optimal LID allocations and analyse them. You can also export them as tif or other formats for further use.
