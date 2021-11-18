## Overview

Investigating mesh attributes in fitting LGCPs using the INLA-SPDE
approach.

![Illustration of the different triangle attributes refered to
below](README_files/figure-markdown_strict/triang_properties-1.png)

**Note**: The radius-edge ratio of a plolygon *τ* is given by
$re = \\frac{R}{l\_{\\text{min}}}$, where R is the circumradius and
*l*<sub>min</sub> is the length of the shortest edge. The radius ratio
$rr = \\frac{r}{R}$, where *r* is the inradius.

## Files

-   `functions.r` contains a number of function to calculate different
    metrics measures and attributes. The main function
    `get_triag_attributes()` takes a `inla.mesh.2d()` object and returns
    a named list of length three:

    -   `edges` a dataframe of mesh segments  
    -   `triangles` a dataframe with, for each triangle in the
        triangulation, the circumcircle and incircle radii, the
        associated centroid locations, and the radius-edge ratio &
        radius ratio,
    -   `angles` is a 3 x nvert dataframe of all triangle interior
        angles
