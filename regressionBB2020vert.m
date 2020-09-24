function [ScallingCoefficients] = regressionBB2020vert(MRSGSL,VertPGA,DistanceBound,SourceToSiteDistance,saNo)

% regressionBB2020vert  - Multiple linear regression analysis for the vertical PGA values
%
% Inputs:
% 1) variable MRSGSL with 203 rows (one for each strong motion record) and 6 columns: 
%    1st column - Magnitudes, M
%    2nd column - Hypocentral Depths, h
%    3rd column - Epicentral distances, E
%    4th column - Dummy variable for Local Soil conditions SL1
%    5th column - Dummy variable for Deep Geology conditions SG1
%    6th column - Dummy variable for Deep Geology conditions SG2
%    The dummy variables are defined in the following manner:
%    “Rock” soil sites:	SL1 = 0	(the local soil conditions)
%    Stiff soil sites: SL1 = 1 (the local soil conditions)
%    Basement (Geological) rock: SG1 = SG2 = 0 (the deep geology conditions)
%    Intermediate sites: SG1 = 1 and SG2 = 0 (the deep geology conditions)
%    Sediments: SG1 = 0 and SG2 = 1 (the deep geology conditions)
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
% 2) variable VertPGA with 203 rows and 1 column. The 203 rows consist of the 
%   PGA values recorded in the vertical direction.
% 3) variable DistanceBound - defines the distance range (in km) of the data used for 
%   the regression analyses. For example, if DistanceBound = 1000, all the data will be used.
% 4) variable SourceToSiteDistance - defines the type of the source to site 
%    distance to be used for R:
%   - SourceToSiteDistance = 1 means that the epicentral distances E will be used.
%   - SourceToSiteDistance = 2 means that the hypocentral distances (h^2+E^2)^0.5 will be used.
% 5) variable saNo - if only empirical scaling equations for the PGA values
%   are to be calculates then saNo = 1 (the default value). If scaling equations 
%   for more strong motion intensity measures are to be calculated then
%   saNo is the number of the equations. In that case the variable VertPGA
%   should have n = saNo number of columns, each with the values of the
%   intensity measure for which the scaling equation is to be defined.
%
% Outputs:
% 1) variable ScallingCoefficients - containing the following values:
% [1] c1, [2] c2, [3] c3, [4] R0, [5] c4, [6] zero value (this value is to be determined 
% in the second phase of the analysis as the scaling coefficient for the deep soil sites), 
% [7] c6, [8] c7, [9] sigma
%
% Copyright (C) Borko Bulajic 2020 -


thedata=MRSGSL;
sa=VertPGA;
thedata=thedata(find(abs(MRSGSL(:,1))<=DistanceBound),:);
sa=[sa(find(abs(MRSGSL(:,1))<=DistanceBound),:)];

% saNo - number of spectral amplitudes: saNo=1 - scaling coefficients are calculated only for PGA
n_sa=length(saNo);

RESULTSvert=[];
RESIDUALSvert=[];
ScallingCoefficients=[];

for j=1:n_sa

	R0values=[0.1:0.1:40]; % Range of values for R0 to be considered in iterative procedure
	n_R0values=length(R0values);

	bcoeffS=[];
	statcoeffS=[];
	RESIDUALS=[];
	Scalingcoeff=[];
	for i=1:n_R0values
		
		if SourceToSiteDistance==1
			%EPICENTRAL DISTANCE
			Cs=[ones(size(thedata,1),1) thedata(:,3) log10((thedata(:,1).^2+(R0values(i))^2).^0.5) thedata(:,5) thedata(:,6) thedata(:,4)];
		elseif SourceToSiteDistance==2
			%HYPOCENTRAL DISTANCE
			Cs=[ones(size(thedata,1),1) thedata(:,3) log10((thedata(:,1).^2+(R0values(i))^2+thedata(:,2).^2).^0.5) thedata(:,5) thedata(:,6) thedata(:,4)];
		end    
		ds=[log10(abs(sa(:,j))./980.7)];
		[b,bint,r,rint,stats] = regress(ds,Cs);
		   
		statcoeffS = [statcoeffS; stats];
		[minS,Is]=min(statcoeffS(:,4));
		optCs=R0values(Is);
		bcoeffS = [bcoeffS; b' (stats(4))^0.5 R0values(i) stats(1:3)];
		Scalingcoeff=[Scalingcoeff;b(1),b(2),b(3),R0values(i),b(6),0,b(4),b(5),(stats(4))^0.5];
		RESIDUALS = [RESIDUALS; r'];

	end

    RESULTSvert1=bcoeffS(Is,:);
    RESULTSvert=[RESULTSvert; RESULTSvert1];
    RESIDUALSvert1=RESIDUALS(Is,:);
    RESIDUALSvert=[RESIDUALSvert; RESIDUALSvert1];
    ScallingCoefficients1=Scalingcoeff(Is,:);
    ScallingCoefficients=[ScallingCoefficients;ScallingCoefficients1];
end



