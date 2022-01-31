Aeroelastic Lumped Parameter (ALP) Modeling Codes
Created by David W. Fellows and Daniel J. Bodony, UIUC

applications -- contains directories including:
    - testing of fluid and structural utilities used to compute various quantities needed to construct reduced-order models
    - files used to obtain reduced-order models of panel flutter scenarios in order to test/partially validate reduced-order modeling framework
    - files used to constuct reduced-order models of turbochargers
    --- the directories contained here include driver files needed to perform stability computations and compute damping ratios

struct_energies_utils -- contains subroutines that compute structural energy coefficients

fluid_function_utils -- contains subroutines that are needed to compute quanities related to fluid forcing functions

mesh_interpolation_utils -- contains subroutines that are needed to perform mesh-to-mesh interpolation of fluid mesh onto structural mesh

hidden directories: -- not included in this repository, included for author's completeness
    .fluid_model/: the above directories are re-organizations of this directory, in order to make the code appaear cleaner
                   this directory is maintained only to ensure no files have been lost

