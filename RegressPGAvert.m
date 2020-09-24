% AN EXAMPLE - Multiple linear regression analysis for the vertical PGA values
%
% Inputs:
% 1) MRSGSLdata.mat - containing a variable MRSGSL with 203 rows 
%                   (one for each strong motion record) and 6 columns: 
%                               1st column - Magnitudes, M
%                               2nd column - Hypocentral Depths, h
%                               3rd column - Epicentral distances, E
%                               4th column - Dummy variable for Local Soil conditions SL1
%                               5th column - Dummy variable for Deep Geology conditions SG1
%                               6th column - Dummy variable for Deep Geology conditions SG2
%   The dummy variables are defined in the following manner:
%                               “Rock” soil sites:	SL1 = 0	
%                               Stiff soil sites: SL1 = 1
%                               Basement (Geological) rock: SG1 = SG2 = 0
%                               Intermediate sites: SG1 = 1 and SG2 = 0
%                               Sediments: SG1 = 0 and SG2 = 1
%   The empirical scaling equation is used to estimate decadic logarithm of
%   the horizontal PGA value (in g), as the following sum: c1 + c2*M +
%   + c3*log[(R^2+R0^2)^0.5] + c4*SL1 + c6*SG1 + c7*SG2 + sigma*eps
%   where c1, c2, c3, c4, c6, and c7 are the unknown coefficients to be obtained by 
%   the multiple linear regression. The R0 values were iteratively adjusted to 
%   maximize the R2 statistics of the PGA prediction. For all regression analyses 
%   it is assumed that the analyzed data follow the log-normal distribution. 
%   Hence, sigma is the standard deviation for the common logarithm of PGA, 
%   where eps = 0 for the median estimates and eps = 1 for the median +1*sigma. 
%   R can be either the epicentral distance, E, or the hypocentral distance (h^2+E^2)^0.5
% 2) VertPGAdata.mat - containing a variable VertPGA with 203 rows and 1 column. 
%   The 203 rows consist of the PGA values recorded in the vertical direction.
%
% Outputs:
% 1) ScallingCoefficientsVERT.mat - containing the input variables MRSGSL and HorizPGA, 
% as well as the four variables:
% - ScallingCoefficientsEPICallR - scaling coefficients for the equation with R as the 
%   epicentral distance, based on the data recorded at all source to site distances
% - ScallingCoefficientsEPICsmallerthan30km - scaling coefficients for the equation with R as the 
%   epicentral distance, based on the data recorded at the epicentral distances smaller than 30 km
% - ScallingCoefficientsHYPOCallR - scaling coefficients for the equation with R as the 
%   hypocentral distance, based on the data recorded at all source to site distances
% - ScallingCoefficientsHYPOCsmallerthan30km - scaling coefficients for the equation with R as the 
%   hypocentral distance, based on the data recorded at the epicentral distances smaller than 30 km
% Each of the four variables contains the following values:
% [1] c1, [2] c2, [3] c3, [4] R0, [5] c4, [6] zero value (this value is to be determined 
% in the second phase of the analysis as the scaling coefficient for the deep soil sites), 
% [7] c6, [8] c7, [9] sigma
%
% Copyright (C) Borko Bulajic 2020 -


clear all 
clc

load MRSGSLdata.mat
load VertPGAdata.mat

% EPICENTRAL DISTANCE, THE DATA FOR ALL DISTANCES
DistanceBound=1000;
SourceToSiteDistance=1;
ScallingCoefficientsEPICallR=regressionBB2020vert(MRSGSL,VertPGA,DistanceBound,SourceToSiteDistance,1);

% HYPOCENTRAL DISTANCE, THE DATA FOR ALL DISTANCES
DistanceBound=1000;
SourceToSiteDistance=2;
ScallingCoefficientsHYPOCallR=regressionBB2020vert(MRSGSL,VertPGA,DistanceBound,SourceToSiteDistance,1);

% EPICENTRAL DISTANCE, THE DATA FOR EPICENTRAL DISTANCES SMALLER THAN 30km
DistanceBound=30;
SourceToSiteDistance=1;
ScallingCoefficientsEPICsmallerthan30km=regressionBB2020vert(MRSGSL,VertPGA,DistanceBound,SourceToSiteDistance,1);

% HYPOCENTRAL DISTANCE, THE DATA FOR EPICENTRAL DISTANCES SMALLER THAN 30km
DistanceBound=30;
SourceToSiteDistance=2;
ScallingCoefficientsHYPOCsmallerthan30km=regressionBB2020vert(MRSGSL,VertPGA,DistanceBound,SourceToSiteDistance,1);

clear DistanceBound SourceToSiteDistance

save ScallingCoefficientsVERT.mat