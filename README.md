Investigating mesh attributes in fitting LGCPs using the INLA-SPDE
approach.

    ## deldir 0.2-10      Nickname: "Morpheus and Euripides"

    ## 
    ##      Note 1: As of version 0.2-1, error handling in this 
    ##      package was amended to conform to the usual R protocol. 
    ##      The deldir() function now actually throws an error 
    ##      when one occurs, rather than displaying an error number 
    ##      and returning a NULL.
    ##  
    ##      Note 2:  As of version 0.1-29 the arguments "col" 
    ##      and "lty" of plot.deldir() had their names changed to 
    ##      "cmpnt_col" and "cmpnt_lty" respectively basically 
    ##      to allow "col" and and "lty" to be passed as "..." 
    ##      arguments.
    ##  
    ##      Note 3: As of version 0.1-29 the "plotit" argument 
    ##      of deldir() was changed to (simply) "plot".
    ##  
    ##      See the help for deldir() and plot.deldir().

![Illustration of the different triangle attributes refered to
below](README_files/figure-markdown_strict/triang_properties-1.png)

## Overview

## Files

-   `functions.r` contains a number of function to calculate different
    metrics measures and attributes. The main function
    `get_triag_attributes()` takes a `inla.mesh.2d()` object and returns
    a named list of length three:

    -   `edges` a dataframe of mesh segments  
    -   `triangles` a dataframe with, for each triangle in the
        triangulation, the circumcircle and incircle radii, the
        assocoated centroid locations, and the radius-edge ratio &
        radius ratio,
    -   angles\` is a 3 x nvert dataframe of all triangle interior
        angles
