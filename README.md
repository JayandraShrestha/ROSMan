# ROSMan
# A new hydropower reservoir operation routine of SWAT for improved management of flows, energy production and sediment.


# Description
The Reservoir Operation and Sediment MANagement (ROSMan) routine (Shrestha et al. 2020; Shrestha et al. 2021) is a new reservoir routine implemented in SWAT (Soil and Water Assessment Tool). The ROSMan has fundamentally two capabilities: 1) hydropower reservoir operations without considering sedimentation (HydROR) (Shrestha et al. 2020) and 2) accumulation and removal of sediment under hydropower reservoir operations and sediment management techniques (ResSMan). Thus, the user can choose between the HydROR and ResSMan. If the main purpose of simulation is to predict energy generation and impacts on the hydrologic regime of a river due to operation of hydropower reservoirs under different policies at the river basin scale, then it is suggested to simulate with the HydROR. The HydROR calculates the water balance of a reservoir and energy generation of a hydropower plant using predefined rule curves and plant efficiency without considering the impacts of sedimentation on the storage capacity. On the other hand, if the main objective of the study is to assess the accumulation of sediment and its impact on reservoir storage capacity, then the user must simulate with the ResSMan. The ResSMan routine has capabilities to predict the accumulation of trapped sediment, its impacts on the storage-capacity of a reservoir, and losses in hydropower generation under user-specified operation policies. Furthermore, it allows to compute the restoration of storage volume due to the removal of sediment by flushing (removal of sediment from a reservoir by passing water and sediment through flush gates located at the low level of a dam) and sluicing (passing sediment before suspended sediment solids have settled down in reservoirs).


# Quick start guide
Set up a SWAT model for the study area with assigning a reservoir/s.

Load all the required input files for the ROSMan, as described in the user’s manual, in the SWAT “TxtInOut” folder. 

Placed the executable file of the ROSMan (rosman.exe) also in the “TxtInOut” folder. 

Run the model by clicking the rosman.exe file.


# Contact
Jayandra P Shrestha
shresthajayandra@gmail.com

Citation: Shrestha, J. P. (2021). "Development and application of a SWAT hydropower routine for flow, sediment, and energy management." PhD Thesis, University of Canterbury, Christchurch, New Zealand.
DOI: http://dx.doi.org/10.26021/11342

Related publication: Shrestha, J. P., Pahlow, M., and Cochrane, T. A. (2020). "Development of a SWAT Hydropower Operation Routine and Its Application to Assessing Hydrological Alterations in the Mekong." Water, 12(8), 2193.
DOI: https://doi.org/10.3390/w12082193

