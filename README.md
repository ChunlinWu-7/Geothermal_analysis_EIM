# Geothermal engineering -- spherical water tanks under time-harmonic heat transfer

## Description
This software package serves as supplemental information for Reviewers and Editorial board members to review/verify results/figures/findings in the manuscript "The Green's Function Based Analysis for Seasonal Heat Storage and Utilization in a Solar-Geothermal Building". 

This software package is not open-access currently, which will be available after acceptance of the above manuscript. 

## Supplemental library
Eigen-3.4.0 and boost library are applied to conduct analysis. 

Please install these two libraries before compiling the code. 

## Intput files:
(1) soil_properties: specify the thermal conductivity and heat capcity of the soil;

(2) particle.txt: specify number of water tanks; positions of tanks; effective thermal conductivity; effective heat capacity; 

(3) post_info.txt: specify post-process points.

## Output files:
(1) post_temp.txt: postprocess temperature; 

(2) post_flux.txt: postprocess heat flux; 

(3) post_temp_ori.txt: annual temperature profile.

## Contributor

Chunlin Wu, Assistant Professor, Shanghai Institute of Applied Mathematics and Mechanics; 

Huiming Yin, Associate Professor, Columbia University. 

© 2025 Chunlin Wu. All rights reserved.  

## To cite this work:
Chunlin Wu, Tengxiang Wang, Huiming Yin, The Green's Function Based Analysis for Seasonal Heat Storage and Utilization in a Solar-Geothermal Building, Journal of Energy Storage. 
