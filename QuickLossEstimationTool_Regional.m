clear
close all

% Initialize values for user inputs (these can be modified as needed)
BuildingOccupancy     = 'RES1';            % Example: Selected building occupancy % Building Occupancy RES1,RES2,RES3A,RES3B,RES3C,RES3D,RES3E,RES3F,RES4,RES5,RES6,COM1,COM2,COM3,COM4A,COM4B,COM4C,COM5,COM6,COM7A,COM7B,COM7C,COM8,COM9,COM10,IND1,IND2,IND3,IND4,IND5,IND6,AGR1,REL1,GOV1,GOV2,EDU1,EDU2
BuildingType          = 'W1';              % Example: Selected building type % Bilding Type W1,W2,W3,W4,C1H,C1M,C1L,C2H,C2M,C2L,C3H,C3M,C3L,MH,PC1,PC2H,PC2M,PC1L,RM1M,RM1L,RM2H,RM2M,RM2L,S1H,S1M,S1L,S2H,S2M,S2L,S3,S4H,S4M,S4L,URMM,URML
DesignLevel           = 'PC';              % Example: Selected design level % Design Level PC,LC,MC,HC 
site_of_interest      = 'Toronto';        % Victoria, Vancouver, Calgary, Toronto, Ottawa, Montreal, Quebec, LaMalbaie
Vs30msPick            = 'Constant';        % 'Global  ', 'Constant', 'Regional' % it shoudl be 8 characters
Vs30msValue           = 360;               % Constant Vs30 in m/s
StructuralValueCAD    = 120000;            % Example: User input for structural value in CAD
NonStructuralValueCAD = 360000;            % Example: User input for non-structural value in CAD
ContentValueCAD       = 220000;            % Example: User input for content value in CAD
Return_period         = 2500;              % 100, 200, 500, 1000, 2500, 5000, 10000
NumberofSites         = '121 (11*11)    '; % '121 (11*11)    ', '441 (21*11)    ', '961 (31*31)    ', '1681 (41*41)   ', '2601 (51*51)   ', '10201 (101*101)' % it should be 15 characters
NumberSample          = 100000;            % Number of samples above the minimum im level

load('coastline_Canada_GEODAS.mat');
load('Vulnerability_NRC2019.mat');
load('fsauid_data.mat','Exposure'); 
load('VS30_SLL.mat','VS30_SLL');
load('MetroVancouver_Vs30Data.mat','MetroVancouver_Vs30Data');

[x_global_vs30,y_global_vs30,z_global_vs30] = grdread2("global_vs30.grd"); % Global Vs30 data

load('GSC2020_HazardMap_Total_Classified_RegionalVsAdded.mat','PE','Site', ...
    'Hazard_GSC2020_140' ,'Hazard_GSC2020_160' ,'Hazard_GSC2020_180' ,'Hazard_GSC2020_250' ,'Hazard_GSC2020_300', ...
    'Hazard_GSC2020_360' ,'Hazard_GSC2020_450' ,'Hazard_GSC2020_580' ,'Hazard_GSC2020_760' ,'Hazard_GSC2020_910', ...
    'Hazard_GSC2020_1100','Hazard_GSC2020_1500','Hazard_GSC2020_1600','Hazard_GSC2020_2000','Hazard_GSC2020_3000');

Vs30_vector = [140 160 180 250 300 360 450 580 760 910 1100 1500 1600 2000 3000]';

disp('                                                                                             ');
disp('--- Seismic Risk Analysis for Canadian Buildings Using OpenQuake-Canada Model Information ---');
disp('---                                    Regional Option                                    ---');
disp('                                                                                             ');

rng(0)

%% 

sites_of_interest = site_of_interest;

% Define coordinates for each city [left longitude, right longitude, bottom latitude, top latitude]
switch char(sites_of_interest)
    case 'Victoria'
        region_coordinates = [-123.59, -123.14, 48.21, 48.66];
    case 'Vancouver'
        region_coordinates = [-123.33, -122.88, 49.03, 49.48]; % Original
        region_coordinates = [-123.30, -122.75, 48.95, 49.50]; % Modified
    case 'Calgary'
        region_coordinates = [-114.29, -113.84, 50.79, 51.24];
    case 'Toronto'
        region_coordinates = [-79.60, -79.15, 43.48, 43.93];
    case 'Ottawa'
        region_coordinates = [-76.05, -75.60, 45.08, 45.53];
    case 'Montreal'
        region_coordinates = [-73.93, -73.48, 45.32, 45.77]; % Original
        region_coordinates = [-74.00, -73.45, 45.25, 45.80]; % Modified
    case 'Quebec'
        region_coordinates = [-71.55, -71.10, 46.63, 47.08];
    case 'LaMalbaie'
        region_coordinates = [-70.47, -70.02, 47.47, 47.92];
    otherwise
        % Default values if no match is found
        disp('Unknown site of interest. Please select a site.');
end

% Display the region coordinates for verification
disp(['Selected Region Coordinates: ',num2str(region_coordinates)]);
            
% Add resolution
if     NumberofSites == '121 (11*11)    '
    desired_resolution = 100;                     
elseif NumberofSites == '441 (21*11)    '
    desired_resolution = 400;
elseif NumberofSites == '961 (31*31)    '
    desired_resolution = 900;
elseif NumberofSites == '1681 (41*41)   '
    desired_resolution = 1600;
elseif NumberofSites == '2601 (51*51)   '
    desired_resolution = 2500;
elseif NumberofSites == '10201 (101*101)'
    desired_resolution = 10000;
end

desired_interval = (min((round(region_coordinates(4),2)-round(region_coordinates(3),2)),(round(region_coordinates(2),2)-round(region_coordinates(1),2))))/sqrt(desired_resolution);

[region_sites_of_interest_lat,region_sites_of_interest_lon] = meshgrid(round(region_coordinates(3),2):desired_interval:round(region_coordinates(4),2),round(region_coordinates(1),2):desired_interval:round(region_coordinates(2),2));
region_sites_of_interest = [region_sites_of_interest_lat(:) region_sites_of_interest_lon(:)];

num_region_sites_of_interest = size(region_sites_of_interest,1);

%% Exposure information

Exposure.id{1}         = {};
Exposure.number(1)     = 1;
Exposure.day(1)        = 0;
Exposure.night(1)      = 0;
Exposure.transit(1)    = 0;
Exposure.landusetyp{1} = {};
Exposure.GenOcc{1}     = {};
Exposure.GenOcc{1}     = {};
Exposure.hzOccType{1}  = {};
Exposure.eqOccType{1}  = {};
Exposure.BldgGen{1}    = {};
Exposure.hzTaxon{1}    = {};
Exposure.EqBldgType{1} = {};
Exposure.EqDesLev{1}   = {};
Exposure.sauid(1)      = 0;
Exposure.dauid(1)      = 0;
Exposure.adauid(1)     = 0;
Exposure.cdname{1}     = {};
Exposure.SAC(1)        = 0;
Exposure.ername{1}     = {};
Exposure.prname{1}     = {};
Exposure.csdname{1}    = {};

% Initialize exposure data
Exposure.structural    = StructuralValueCAD;
Exposure.nonstructural = NonStructuralValueCAD;
Exposure.contents      = ContentValueCAD;
Exposure.eqOccType     = {BuildingOccupancy};
Exposure.EqBldgType    = {BuildingType};
Exposure.EqDesLev      = {DesignLevel};
Exposure.taxonomy{1}   = char({[char(Exposure.eqOccType),'-',char(Exposure.EqBldgType),'-',char(Exposure.EqDesLev)]});
          
%% Seismic hazard and vulnerability information

F_GSC = 1 - PE(1,:)'; 

VulnerabilityHazard = zeros(num_region_sites_of_interest,22);

map_margin = 0.05;

map_range = [region_coordinates(1)-map_margin region_coordinates(2)+map_margin region_coordinates(3)-map_margin region_coordinates(4)+map_margin];

buffer = [-1 1 -1 1]; % Latitude west-east, longitude south-north

% Define the latitude and longitude limits for trimming
lat_min = region_coordinates(3) + buffer(1);
lat_max = region_coordinates(4) + buffer(2);
lon_min = region_coordinates(1) + buffer(3);
lon_max = region_coordinates(2) + buffer(4);

filtered_Site = Site(Site(:,1)>=lat_min & Site(:,1)<=lat_max & Site(:,2)>=lon_min & Site(:,2)<=lon_max,:);

[~,filtered_Indices] = ismember(filtered_Site,Site,'rows');  

filtered_Hazard_GSC2020_140  = Hazard_GSC2020_140(filtered_Indices,:,:);
filtered_Hazard_GSC2020_160  = Hazard_GSC2020_160(filtered_Indices,:,:);
filtered_Hazard_GSC2020_180  = Hazard_GSC2020_180(filtered_Indices,:,:);
filtered_Hazard_GSC2020_250  = Hazard_GSC2020_250(filtered_Indices,:,:);
filtered_Hazard_GSC2020_300  = Hazard_GSC2020_300(filtered_Indices,:,:);
filtered_Hazard_GSC2020_360  = Hazard_GSC2020_360(filtered_Indices,:,:);
filtered_Hazard_GSC2020_450  = Hazard_GSC2020_450(filtered_Indices,:,:);
filtered_Hazard_GSC2020_580  = Hazard_GSC2020_580(filtered_Indices,:,:);
filtered_Hazard_GSC2020_760  = Hazard_GSC2020_760(filtered_Indices,:,:);
filtered_Hazard_GSC2020_910  = Hazard_GSC2020_910(filtered_Indices,:,:);
filtered_Hazard_GSC2020_1100 = Hazard_GSC2020_1100(filtered_Indices,:,:);
filtered_Hazard_GSC2020_1500 = Hazard_GSC2020_1500(filtered_Indices,:,:);
filtered_Hazard_GSC2020_1600 = Hazard_GSC2020_1600(filtered_Indices,:,:);
filtered_Hazard_GSC2020_2000 = Hazard_GSC2020_2000(filtered_Indices,:,:);
filtered_Hazard_GSC2020_3000 = Hazard_GSC2020_3000(filtered_Indices,:,:);

clear Hazard_GSC2020_140 Hazard_GSC2020_160 Hazard_GSC2020_180 Hazard_GSC2020_250 Hazard_GSC2020_300 Hazard_GSC2020_360 Hazard_GSC2020_450 Hazard_GSC2020_580 Hazard_GSC2020_760 Hazard_GSC2020_910 Hazard_GSC2020_1100 Hazard_GSC2020_1500 Hazard_GSC2020_1600 Hazard_GSC2020_2000 Hazard_GSC2020_3000

% Identify the vulnerability function
if contains(string(Exposure.taxonomy(1)),'W3') == 1
    taxonomy_tmp = char(string(Exposure.taxonomy(1)));
    taxonomy_tmp(length(taxonomy_tmp)-3) = '2';
    vul_pick = find(contains(string(VulClassName(:,1)),taxonomy_tmp) == 1);
    clear taxonomy_tmp
elseif contains(string(Exposure.taxonomy(1)),'W4') == 1 
    taxonomy_tmp = char(string(Exposure.taxonomy(1)));
    taxonomy_tmp(length(taxonomy_tmp)-3) = '1';
    vul_pick = find(contains(string(VulClassName(:,1)),taxonomy_tmp) == 1);
    clear taxonomy_tmp
else
    vul_pick = find(contains(string(VulClassName(:,1)),string(Exposure.taxonomy(1))) == 1);
end
            
% Identify the intensity parameter for the vulnerability function
if     strcmp(string(VulIMName(vul_pick,1)),'SA(10.0)') == 1; t_pick = 1;
elseif strcmp(string(VulIMName(vul_pick,1)),'SA(5.0)')  == 1; t_pick = 2;
elseif strcmp(string(VulIMName(vul_pick,1)),'SA(2.0)')  == 1; t_pick = 3;
elseif strcmp(string(VulIMName(vul_pick,1)),'SA(1.0)')  == 1; t_pick = 4;
elseif strcmp(string(VulIMName(vul_pick,1)),'SA(0.5)')  == 1; t_pick = 5;
elseif strcmp(string(VulIMName(vul_pick,1)),'SA(0.3)')  == 1; t_pick = 6;
elseif strcmp(string(VulIMName(vul_pick,1)),'SA(0.2)')  == 1; t_pick = 7;
elseif strcmp(string(VulIMName(vul_pick,1)),'SA(0.1)')  == 1; t_pick = 8;
elseif strcmp(string(VulIMName(vul_pick,1)),'SA(0.05)') == 1; t_pick = 9;
elseif strcmp(string(VulIMName(vul_pick,1)),'PGA')      == 1; t_pick = 10;
elseif strcmp(string(VulIMName(vul_pick,1)),'PGV')      == 1; t_pick = 11;
end                

%% Creation of global and regional Vs30 maps in a region of interest

num_gscsite = size(filtered_Site,1);

filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,1 ) = filtered_Hazard_GSC2020_140(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,2 ) = filtered_Hazard_GSC2020_160(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,3 ) = filtered_Hazard_GSC2020_180(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,4 ) = filtered_Hazard_GSC2020_250(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,5 ) = filtered_Hazard_GSC2020_300(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,6 ) = filtered_Hazard_GSC2020_360(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,7 ) = filtered_Hazard_GSC2020_450(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,8 ) = filtered_Hazard_GSC2020_580(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,9 ) = filtered_Hazard_GSC2020_760(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,10) = filtered_Hazard_GSC2020_910(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,11) = filtered_Hazard_GSC2020_1100(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,12) = filtered_Hazard_GSC2020_1500(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,13) = filtered_Hazard_GSC2020_1600(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,14) = filtered_Hazard_GSC2020_2000(:,:,t_pick);
filtered_Hazard_GSC2020_combined(1:num_gscsite,1:10,15) = filtered_Hazard_GSC2020_3000(:,:,t_pick);

count_map = 0;
for kk = 1:length(Vs30_vector)
    for jj = 1:length(PE)

        count_map = count_map + 1;

        data = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_combined(:,jj,kk),region_sites_of_interest(:,2),region_sites_of_interest(:,1));

        filtered_Hazard_GSC2020_combined_sites(jj,kk,:) = data;

        % if kk == 6
        %     figure(1000+count_map) % Visual check
        %     scatter(region_sites_of_interest(:,2),region_sites_of_interest(:,1),9,data,'filled'); hold on;
        %     scatter(filtered_Site(:,2),filtered_Site(:,1),36,filtered_Hazard_GSC2020_combined(:,jj,kk),'filled'); hold on;
        %     scatter(filtered_Site(:,2),filtered_Site(:,1),'ko'); hold on;
        %     plot(PolyLon, PolyLat,'k-'); hold on;
        %     axis('square'); axis('equal'); axis(map_range);
        %     title(['RP = ',num2str(PE(2,jj)),' years & Vs30 = ',num2str(Vs30_vector(kk)),' m/s']);
        %     clim([min(data(:)) max(data(:))]); colorbar;
        % end

        clear data

    end
end

% Global Vs30 data
[X_global_vs30_trimmed,Y_global_vs30_trimmed] = meshgrid(x_global_vs30((x_global_vs30>=lon_min) & (x_global_vs30<=lon_max)),y_global_vs30((y_global_vs30>=lat_min) & (y_global_vs30<=lat_max)));
z_global_vs30_trimmed = double(z_global_vs30((y_global_vs30>=lat_min) & (y_global_vs30<=lat_max),(x_global_vs30>=lon_min) & (x_global_vs30<=lon_max)));
z_global_vs30_trimmed(z_global_vs30_trimmed==600) = NaN;
globalVs30Data_trimmed = [Y_global_vs30_trimmed(:), X_global_vs30_trimmed(:), z_global_vs30_trimmed(:)]; globalVs30Data_trimmed_tmp = globalVs30Data_trimmed;
globalVs30Data_trimmed = globalVs30Data_trimmed(~isnan(globalVs30Data_trimmed(:,3)),:);

% Regional Vs30 data
regionVs30Data_trimmed = [VS30_SLL(:,7), VS30_SLL(:,8), VS30_SLL(:,6)];
regionVs30Data_trimmed = [regionVs30Data_trimmed; MetroVancouver_Vs30Data];
regionVs30Data_trimmed = regionVs30Data_trimmed(regionVs30Data_trimmed(:,1)>=lat_min & regionVs30Data_trimmed(:,1)<=lat_max & regionVs30Data_trimmed(:,2)>=lon_min & regionVs30Data_trimmed(:,2)<=lon_max,:);
regionVs30Data_trimmed = regionVs30Data_trimmed(~isnan(regionVs30Data_trimmed(:,3)),:);

% Assign Global and Regional Vs30 data to grids
region_sites_of_interest(:,3) = griddata(globalVs30Data_trimmed(:,2),globalVs30Data_trimmed(:,1),globalVs30Data_trimmed(:,3),region_sites_of_interest(:,2),region_sites_of_interest(:,1));
if ~isempty(regionVs30Data_trimmed)
    region_sites_of_interest(:,4) = griddata(regionVs30Data_trimmed(:,2),regionVs30Data_trimmed(:,1),regionVs30Data_trimmed(:,3),region_sites_of_interest(:,2),region_sites_of_interest(:,1));
else
    region_sites_of_interest(:,4) = NaN;
end
region_sites_of_interest(:,5) = region_sites_of_interest(:,4)./region_sites_of_interest(:,3);
region_sites_of_interest(:,6) = griddata(globalVs30Data_trimmed_tmp(:,2),globalVs30Data_trimmed_tmp(:,1),globalVs30Data_trimmed_tmp(:,3),region_sites_of_interest(:,2),region_sites_of_interest(:,1),'nearest');

[Vs30_grid,PE_grid] = meshgrid(log10(Vs30_vector),PE(2,:));

for ii = 1:num_region_sites_of_interest

    data = filtered_Hazard_GSC2020_combined_sites(:,:,ii);

    if Vs30msPick == 'Global  '  
        im_gsc = griddata(Vs30_grid(:),PE_grid(:),data(:),log10(region_sites_of_interest(ii,3)),PE(2,:)');
    elseif Vs30msPick == 'Regional'   
        im_gsc = griddata(Vs30_grid(:),PE_grid(:),data(:),log10(region_sites_of_interest(ii,4)),PE(2,:)');
        if isnan(im_gsc)
            im_gsc = griddata(Vs30_grid(:),PE_grid(:),data(:),log10(region_sites_of_interest(ii,3)),PE(2,:)');
        end
    elseif Vs30msPick == 'Constant'   
        im_gsc = griddata(Vs30_grid(:),PE_grid(:),data(:),log10(Vs30msValue),PE(2,:)');
    end

    if t_pick == 11; im_gsc = 100*im_gsc; end % PGV (cm/s)

    % Tail approximation of seismic hazard curve: 1) Lognormal, 2) Gumbel, 3) Frechet, 4) Weibull
    RRtmp = zeros(4,4);
    for jj = 1:4
        [~,~,RRtmp(jj,:)] = ModelFit_GSCSeismicHazard_QLET(F_GSC,im_gsc,jj); 
    end
    [~,model_pick] = max(RRtmp(:,1));
    [para_im,coef_im,RR_im] = ModelFit_GSCSeismicHazard_QLET(F_GSC,im_gsc,model_pick);
    
    % Summary hazard-vulnerability information
    % 1    ) Vulnerability function ID (as in Vulnerability_NRC2019)
    % 2    ) Intensity measure ID (as in GSC2020_HazardMap)
    % 3    ) GSC hazard location ID (as in GSC2020_HazardMap)
    % 4 -13) GSC hazard values for 8 return periods (as in GSC2020_HazardMap)
    % 14   ) Probability distribution type for tail approximation 
    % 15-16) Parameters for tail approximation 
    % 17-18) Coefficients for tail approximation 
    % 19-22) Regression output for tail approximation (R^2, F, p-value, error variance)
    VulnerabilityHazard(ii,1:22) = [vul_pick t_pick 0 im_gsc(:,1)' model_pick para_im coef_im' RR_im];

    clear im_gsc RRtmp model_pick para_im coef_im RR_im

end

disp('          ');

%% Seismic vulnerability functions

im_risk_min = 0.05; % (g)

IM_vulfunc = 0.01:0.01:8.0; 

% Interpolation or lognormal approximation of vulnerability functions
if VulIM_ST(vul_pick,VulNumData(vul_pick,1)) >= 8
    VulMU_ST_approx(:,1) = interp1q(VulIM_ST(vul_pick,1:VulNumData(vul_pick,1))',VulMU_ST(vul_pick,1:VulNumData(vul_pick,1))',IM_vulfunc');    
    VulMU_ST_approx(IM_vulfunc<VulIM_ST(vul_pick,1),1) = 0;
else
    [datatmp1(:,1),datatmp2] = unique(VulMU_ST(vul_pick,1:VulNumData(vul_pick,1))); 
    datatmp3(:,1) = VulIM_ST(vul_pick,datatmp2);
    vul_para_approx(1,1) = interp1q(datatmp1,datatmp3,0.50);
    if isnan(vul_para_approx(1,1))
        point1 = [datatmp3(find(datatmp1>0.25,1,'first')) datatmp1(find(datatmp1>0.25,1,'first'))];
        point2 = [datatmp3(end) datatmp1(end)];
        vul_para_approx(1,1) = point1(1) + (0.5-point1(2))/(point2(2)-point1(2))*(point2(1)-point1(1));
        clear point1 point2
    end
    vul_para_approx(1,2) = interp1q(datatmp1,datatmp3,0.16);
    vul_para_approx(1,3) = interp1q(datatmp1,datatmp3,0.84);
    if isnan(vul_para_approx(1,3))
        vul_para_approx(1,4) = log(vul_para_approx(1,1)) - log(vul_para_approx(1,2));
    else
        vul_para_approx(1,4) = (log(vul_para_approx(1,3)) - log(vul_para_approx(1,2)))/2;
    end
    VulMU_ST_approx(:,1) = normcdf(log(IM_vulfunc/vul_para_approx(1,1))/vul_para_approx(1,4));    
    clear datatmp1 datatmp2 datatmp3
end
    
if VulIM_NS(vul_pick,VulNumData(vul_pick,2)) >= 8
    VulMU_NS_approx(:,1) = interp1q(VulIM_NS(vul_pick,1:VulNumData(vul_pick,2))',VulMU_NS(vul_pick,1:VulNumData(vul_pick,2))',IM_vulfunc');    
    VulMU_NS_approx(IM_vulfunc<VulIM_NS(vul_pick,1),1) = 0;
else
    [datatmp1(:,1),datatmp2] = unique(VulMU_NS(vul_pick,1:VulNumData(vul_pick,2))); 
    datatmp3(:,1) = VulIM_NS(vul_pick,datatmp2);
    vul_para_approx(1,5) = interp1q(datatmp1,datatmp3,0.50);
    if isnan(vul_para_approx(1,5))
        point1 = [datatmp3(find(datatmp1>0.25,1,'first')) datatmp1(find(datatmp1>0.25,1,'first'))];
        point2 = [datatmp3(end) datatmp1(end)];
        vul_para_approx(1,5) = point1(1) + (0.5-point1(2))/(point2(2)-point1(2))*(point2(1)-point1(1));
        clear point1 point2
    end
    vul_para_approx(1,6) = interp1q(datatmp1,datatmp3,0.16);
    vul_para_approx(1,7) = interp1q(datatmp1,datatmp3,0.84);
    if isnan(vul_para_approx(1,7))
        vul_para_approx(1,8) = log(vul_para_approx(1,5)) - log(vul_para_approx(1,6));
    else
        vul_para_approx(1,8) = (log(vul_para_approx(1,7)) - log(vul_para_approx(1,6)))/2;
    end
    VulMU_NS_approx(:,1) = normcdf(log(IM_vulfunc/vul_para_approx(1,5))/vul_para_approx(1,8));    
    clear datatmp1 datatmp2 datatmp3
end
    
if VulIM_CO(vul_pick,VulNumData(vul_pick,3)) >= 8
    VulMU_CO_approx(:,1) = interp1q(VulIM_CO(vul_pick,1:VulNumData(vul_pick,3))',VulMU_CO(vul_pick,1:VulNumData(vul_pick,3))',IM_vulfunc');    
    VulMU_CO_approx(IM_vulfunc<VulIM_CO(vul_pick,1),1) = 0;
else
    [datatmp1(:,1),datatmp2] = unique(VulMU_CO(vul_pick,1:VulNumData(vul_pick,3))); 
    datatmp3(:,1) = VulIM_CO(vul_pick,datatmp2);
    vul_para_approx(1,9) = interp1q(datatmp1,datatmp3,0.50);
    if isnan(vul_para_approx(1,9))
        point1 = [datatmp3(find(datatmp1>0.25,1,'first')) datatmp1(find(datatmp1>0.25,1,'first'))];
        point2 = [datatmp3(end) datatmp1(end)];
        vul_para_approx(1,9) = point1(1) + (0.5-point1(2))/(point2(2)-point1(2))*(point2(1)-point1(1));
        clear point1 point2
    end
    vul_para_approx(1,10) = interp1q(datatmp1,datatmp3,0.16);
    vul_para_approx(1,11) = interp1q(datatmp1,datatmp3,0.84);
    if isnan(vul_para_approx(1,11))
        vul_para_approx(1,12) = log(vul_para_approx(1,9)) - log(vul_para_approx(1,10));
    else
        vul_para_approx(1,12) = (log(vul_para_approx(1,11)) - log(vul_para_approx(1,10)))/2;
    end
    VulMU_CO_approx(:,1) = normcdf(log(IM_vulfunc/vul_para_approx(1,9))/vul_para_approx(1,12));    
    clear datatmp1 datatmp2 datatmp3
end

%% Seismic risk analysis

disp('--- Seismic Hazard & Risk Analysis ---');

Psimu = (1:NumberSample)/(NumberSample+1);

Loss_RP = 1./[100 200 500 1000 2500 5000 10000];

Loss_thresh = [1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000]; % dollars

% Uncertainty of vulnerability function - this is applied to structural, nonstructural, and contents individually.
% 1) Lognormal model with a CoV value
% 2) Beta model with a CoV value
option_vulnerability_uncertainty = 1;
if option_vulnerability_uncertainty == 1
    CoV_vul_logn = 0.3;
elseif option_vulnerability_uncertainty == 2
    CoV_vul_beta = 0.3;
    a_beta =  (1-(0.01:0.01:0.99))/CoV_vul_beta^2-(0.01:0.01:0.99);
    b_beta = ((1-(0.01:0.01:0.99))/CoV_vul_beta^2-(0.01:0.01:0.99)).*((1-(0.01:0.01:0.99))./(0.01:0.01:0.99));
    figure (99); plot(0.01:0.01:0.99,a_beta,'b-',0.01:0.01:0.99,b_beta,'r-');
    meanloss_upper = 1/(1+CoV_vul_beta^2);
end

LossSummary1 = zeros(num_region_sites_of_interest,length(Loss_RP));
LossSummary2 = zeros(num_region_sites_of_interest,length(Loss_thresh));

LOSS_all_risk_sites_of_interest = zeros(NumberSample,num_region_sites_of_interest);

for ii = 1:num_region_sites_of_interest
            
    if VulnerabilityHazard(ii,14) == 1
        ur_risk_min = logncdf(im_risk_min,VulnerabilityHazard(ii,15),VulnerabilityHazard(ii,16));
        ur_risk     = ur_risk_min + rand(NumberSample,1)*(1-ur_risk_min);
        im_risk     = exp(VulnerabilityHazard(ii,15)+VulnerabilityHazard(ii,16)*norminv(ur_risk));
    elseif VulnerabilityHazard(ii,14) == 2
        ur_risk_min = exp(-exp(-(im_risk_min-VulnerabilityHazard(ii,15))*VulnerabilityHazard(ii,16)));
        ur_risk     = ur_risk_min + rand(NumberSample,1)*(1-ur_risk_min);
        im_risk     = VulnerabilityHazard(ii,15)-log(-log(ur_risk))/VulnerabilityHazard(ii,16);
    elseif VulnerabilityHazard(ii,14) == 3
        ur_risk_min = exp(-(VulnerabilityHazard(ii,15)/im_risk_min)^VulnerabilityHazard(ii,16));
        ur_risk     = ur_risk_min + rand(NumberSample,1)*(1-ur_risk_min);
        im_risk     = exp(log(VulnerabilityHazard(ii,15))-log(-log(ur_risk))/VulnerabilityHazard(ii,16));
    elseif VulnerabilityHazard(ii,14) == 4
        ur_risk_min = 1-exp(-(im_risk_min/VulnerabilityHazard(ii,15))^VulnerabilityHazard(ii,16));
        ur_risk     = ur_risk_min + rand(NumberSample,1)*(1-ur_risk_min);
        im_risk     = exp(log(VulnerabilityHazard(ii,15))+log(-log(1-ur_risk))/VulnerabilityHazard(ii,16));
    end

    meanloss_st = interp1q(IM_vulfunc',VulMU_ST_approx(:,1),min(max(im_risk,im_risk_min),8));
    meanloss_ns = interp1q(IM_vulfunc',VulMU_NS_approx(:,1),min(max(im_risk,im_risk_min),8));
    meanloss_co = interp1q(IM_vulfunc',VulMU_CO_approx(:,1),min(max(im_risk,im_risk_min),8));
    
    if option_vulnerability_uncertainty == 0
        
        loss_st = meanloss_st;
        loss_ns = meanloss_ns;
        loss_co = meanloss_co;
        
    elseif option_vulnerability_uncertainty == 1
        
        loss_st = lognrnd(log(meanloss_st/sqrt(1+CoV_vul_logn^2)),sqrt(log(1+CoV_vul_logn^2)),NumberSample,1);
        loss_ns = lognrnd(log(meanloss_ns/sqrt(1+CoV_vul_logn^2)),sqrt(log(1+CoV_vul_logn^2)),NumberSample,1);
        loss_co = lognrnd(log(meanloss_co/sqrt(1+CoV_vul_logn^2)),sqrt(log(1+CoV_vul_logn^2)),NumberSample,1);
        
    elseif option_vulnerability_uncertainty == 2
        
        losscase1 = meanloss_st>0 & meanloss_st<meanloss_upper;
        losscase2 = meanloss_ns>0 & meanloss_ns<meanloss_upper;
        losscase3 = meanloss_co>0 & meanloss_co<meanloss_upper;
        
        loss_st(meanloss_st==0,1) = 0;
        loss_ns(meanloss_ns==0,1) = 0;
        loss_co(meanloss_co==0,1) = 0;

        loss_st(losscase1,1) = betarnd((1-meanloss_st(losscase1))/CoV_vul_beta^2-meanloss_st(losscase1),((1-meanloss_st(losscase1))/CoV_vul_beta^2-meanloss_st(losscase1)).*((1-meanloss_st(losscase1))./meanloss_st(losscase1)),length(find(losscase1==1)),1);
        loss_ns(losscase2,1) = betarnd((1-meanloss_ns(losscase2))/CoV_vul_beta^2-meanloss_ns(losscase2),((1-meanloss_ns(losscase2))/CoV_vul_beta^2-meanloss_ns(losscase2)).*((1-meanloss_ns(losscase2))./meanloss_ns(losscase2)),length(find(losscase2==1)),1);
        loss_co(losscase3,1) = betarnd((1-meanloss_co(losscase3))/CoV_vul_beta^2-meanloss_co(losscase3),((1-meanloss_co(losscase3))/CoV_vul_beta^2-meanloss_co(losscase3)).*((1-meanloss_co(losscase3))./meanloss_co(losscase3)),length(find(losscase3==1)),1);
        
        loss_st(meanloss_st>=meanloss_upper,1) = 1;
        loss_ns(meanloss_ns>=meanloss_upper,1) = 1;
        loss_co(meanloss_co>=meanloss_upper,1) = 1;
    
        clear losscase1 losscase2 losscase3
        
    end
    
    LOSS_all_risk = Exposure.structural*loss_st + Exposure.nonstructural*loss_ns + Exposure.contents*loss_co;
        
    LOSS_all_risk = sort(LOSS_all_risk);
    
    loss_ex_prob = (1-ur_risk_min) - (1-ur_risk_min)*Psimu;
    
    % Loss fractile (VaR) at Loss_RP
    for jj = 1:size(Loss_RP,2)
        
        [unique_ur_risk,unique_indices] = unique(ur_risk,'stable'); % Keep only unique values
        unique_im_risk = im_risk(unique_indices); % Apply the same filter to im_risk      
        
        Hazard_region_sites_of_interest(ii,jj) = interp1(1-sort(unique_ur_risk),sort(unique_im_risk),Loss_RP(1,jj));
     
        loss_ex_prob_pick = find(loss_ex_prob<Loss_RP(jj),1,'first');
        if loss_ex_prob_pick > 1
            loss_frac(jj,1) = LOSS_all_risk(loss_ex_prob_pick);
        else
            loss_frac(jj,1) = 0;
        end

        clear loss_ex_prob_pick

    end    

    % Loss frequency exceeding Loss_thresh 
    for jj = 1:length(Loss_thresh)
        loss_prob_thresh(jj,1) = (1-ur_risk_min)/NumberSample*length(find(LOSS_all_risk>=Loss_thresh(jj)));
    end    
        
    % LossSummary1: Loss fractiles (as specified in Loss_RP)
    LossSummary1(ii,:) = loss_frac(:);
    
    % LossSummary2: Loss frequency 
    LossSummary2(ii,:) = loss_prob_thresh(:);
    
    LOSS_all_risk_sites_of_interest(:,ii) = LOSS_all_risk;
    
    clear ur_risk_min ur_risk im_risk unique_ur_risk unique_indices unique_im_risk
    clear meanloss_st meanloss_ns meanloss_co loss_st loss_ns loss_co 
    clear LOSS_all_risk loss_ex_prob loss_frac loss_prob_thresh

    % Update progress every 1000 iterations or when the loop completes
    if mod(ii,10) == 0 || ii == num_region_sites_of_interest
        progress_percentage = ii/num_region_sites_of_interest*100;
        disp(['Calculating... ', num2str(round(progress_percentage,1)), '% completed.']);
    end    

end

%% Display: Seismic vulnerability functions
figure (1)
plot( ...
    VulIM_ST(vul_pick,1:VulNumData(vul_pick,1)),VulMU_ST(vul_pick,1:VulNumData(vul_pick,1)),'bo', ...
    VulIM_NS(vul_pick,1:VulNumData(vul_pick,2)),VulMU_NS(vul_pick,1:VulNumData(vul_pick,2)),'rs', ...
    VulIM_CO(vul_pick,1:VulNumData(vul_pick,3)),VulMU_CO(vul_pick,1:VulNumData(vul_pick,3)),'g^', ...
    IM_vulfunc,VulMU_ST_approx(:,1),'b-',IM_vulfunc,VulMU_NS_approx(:,1),'r-',IM_vulfunc,VulMU_CO_approx(:,1),'g-','LineWidth',1);
xlabel(string(VulIMName(vul_pick,1))); axis square; ax = gca; ax.XLim = [0 8]; ax.YLim = [0 1];
title(['Vulnerability function - ',char(string(VulClassName(vul_pick,1)))])
legend('Structural data','NonStructural data','Content data','Structural fit','NonStructural fit','Content fit','Location','southeast');    

%% Display: site information

figure (100)
plot(filtered_Site(:,2),filtered_Site(:,1),'b.'); hold on;
plot(PolyLon, PolyLat,'k-'); hold on;
axis('square'); axis('equal'); axis(map_range);
title('GSC hazard locations');

figure (101)
scatter(globalVs30Data_trimmed(:,2),globalVs30Data_trimmed(:,1),16,globalVs30Data_trimmed(:,3),'filled'); hold on;
plot(PolyLon, PolyLat,'k-'); hold on;
axis('square'); axis('equal'); axis(map_range);
title('Global Vs30 locations');
colormap('jet'); colorbar('Direction','reverse'); clim([0 1000]); 

figure (102)
scatter(regionVs30Data_trimmed(:,2),regionVs30Data_trimmed(:,1),16,regionVs30Data_trimmed(:,3),'filled'); hold on;
plot(PolyLon, PolyLat,'k-'); hold on;
axis('square'); axis('equal'); axis(map_range);
title('Regional Vs30 locations');
colormap('jet'); colorbar('Direction','reverse'); clim([0 1000]); 

figure (103)
scatter(region_sites_of_interest(:,2),region_sites_of_interest(:,1),16,region_sites_of_interest(:,3),'filled'); hold on;
plot(PolyLon, PolyLat,'k-'); hold on;
axis('square'); axis('equal'); axis(map_range);
title('Global Vs30 locations [interpolated]');
colormap('jet'); colorbar('Direction','reverse'); clim([0 1000]); 

figure (104)
scatter(region_sites_of_interest(:,2),region_sites_of_interest(:,1),16,region_sites_of_interest(:,4),'filled'); hold on;
plot(PolyLon, PolyLat,'k-'); hold on;
axis('square'); axis('equal'); axis(map_range);
title('Regional Vs30 locations [interpolated]');
colormap('jet'); colorbar('Direction','reverse'); clim([0 1000]); 

figure (105)
scatter(region_sites_of_interest(:,2),region_sites_of_interest(:,1),16,log10(region_sites_of_interest(:,5)),'filled'); hold on;
plot(PolyLon, PolyLat,'k-'); hold on;
axis('square'); axis('equal'); axis(map_range);
title('Regional/Global Vs30 ratios [log10, interpolated]');
colormap('jet'); colorbar; clim([-1 1]); 

%% Seismic hazard and risk maps
if Return_period == 100
    Hazard_level = 1;
elseif Return_period == 200
    Hazard_level = 2;
elseif Return_period == 500
    Hazard_level = 3;
elseif Return_period == 1000
    Hazard_level = 4;
elseif Return_period == 2500
    Hazard_level = 5;
elseif Return_period == 5000
    Hazard_level = 6;
elseif Return_period == 10000
    Hazard_level = 7;
end

if Vs30msPick == 'Regional'
    % Mask hazard and risk for sites without regional Vs30
    for i = 1:num_region_sites_of_interest 
        if isnan(region_sites_of_interest(i,4))
            Hazard_region_sites_of_interest(i,:) = NaN;
            LossSummary1(i,:) = NaN;
        end
    end
end
for i = 1:num_region_sites_of_interest 
    if isnan(region_sites_of_interest(i,6))
        Hazard_region_sites_of_interest(i,:) = NaN;
        LossSummary1(i,:) = NaN;
    end
end

[xq,yq] = meshgrid(linspace(min(region_sites_of_interest(:,2)),max(region_sites_of_interest(:,2))),linspace(min(region_sites_of_interest(:,1)),max(region_sites_of_interest(:,1))));

figure (301)
scatter(region_sites_of_interest(:,2),region_sites_of_interest(:,1),16,Hazard_region_sites_of_interest(:,Hazard_level),'filled'); hold on;
% F = scatteredInterpolant(region_sites_of_interest(:,2),region_sites_of_interest(:,1),Hazard_region_sites_of_interest(:,Hazard_level));
% contourf(xq,yq,F(xq,yq),1000,'LineColor','none'); hold('on'); 
plot(PolyLon,PolyLat,'k-'); hold('on');
colormap('jet'); colorbar(); 
axis('square'); axis('equal'); axis(map_range);
title(['Location of grids and ',num2str(1/Loss_RP(Hazard_level)),'-year SA(0.3) (g)']); 
xlabel('Longitude'); ylabel('Latitude');
hold('off'); 

figure (302)
scatter(region_sites_of_interest(:,2),region_sites_of_interest(:,1),16,(1/10^6).*LossSummary1(:,Hazard_level),'filled'); hold on;
% F = scatteredInterpolant(region_sites_of_interest(:,2),region_sites_of_interest(:,1),(1/10^6).*LossSummary1(:,Hazard_level));
% contourf(xq,yq,F(xq,yq),1000,'LineColor','none'); hold('on'); 
plot(PolyLon,PolyLat,'k-'); grid on;
colormap('jet'); colorbar(); 
axis('square'); axis('equal'); axis(map_range);
title([num2str(1/Loss_RP(Hazard_level)),'-year VaR loss (million CAD)']); 
xlabel('Longitude'); ylabel('Latitude');
hold('off');

%% Add NRCan building maps
switch char(sites_of_interest)

    case 'Victoria'
        Exposure_NRCan = readtable('BC.csv');                    
    case 'Vancouver'
        Exposure_NRCan = readtable('BC.csv');
    case 'Calgary'
        Exposure_NRCan = readtable('AB.csv');                   
    case 'Toronto'
        Exposure_NRCan = readtable('ON.csv');                   
    case 'Ottawa'
        Exposure_NRCan = readtable('ON.csv');                  
    case 'Montreal'
        Exposure_NRCan = readtable('QC.csv');                   
    case 'Quebec'
        Exposure_NRCan = readtable('QC.csv');                    
    case 'LaMalbaie'
        Exposure_NRCan = readtable('QC.csv');    

end

Exposure_NRCan_filtered = Exposure_NRCan(Exposure_NRCan.lon >= lon_min & Exposure_NRCan.lon <= lon_max & Exposure_NRCan.lat >= lat_min & Exposure_NRCan.lat <= lat_max,:);

Exposure_NRCan_filtered_wood     = Exposure_NRCan_filtered(strcmp(Exposure_NRCan_filtered.BldgGen,'Wood'),:);
Exposure_NRCan_filtered_concrete = Exposure_NRCan_filtered(strcmp(Exposure_NRCan_filtered.BldgGen,'Concrete'),:);

figure (303)
scatter(Exposure_NRCan_filtered_wood.lon,Exposure_NRCan_filtered_wood.lat,9,Exposure_NRCan_filtered_wood.number,'filled'); hold('on'); 
plot(PolyLon, PolyLat, 'k-'); grid on; box on;
colormap('jet'); colorbar(); 
axis('square'); axis('equal'); axis(map_range);
title(['Reference NRCan Map of Wood Buildings - Total: ',num2str(sum(Exposure_NRCan_filtered_wood.number)),' Buildings']); 
xlabel('Longitude'); ylabel('Latitude');
clim([0 10]); hold('off');

figure (304)
scatter(Exposure_NRCan_filtered_concrete.lon,Exposure_NRCan_filtered_concrete.lat,9,Exposure_NRCan_filtered_concrete.number,'filled'); hold('on'); 
plot(PolyLon,PolyLat,'k-'); grid on; box on;
colormap('jet'); colorbar(); 
axis('square'); axis('equal'); axis(map_range);
title(['Reference NRCan Map of Concrete Buildings - Total: ',num2str(sum(Exposure_NRCan_filtered_concrete.number)),' Buildings']); 
xlabel('Longitude'); ylabel('Latitude');
clim([0 10]); hold('off');
                 
disp('          ');

%% M-function for fitting probabilistic models to the GSC's hazard values

function [para,coef,corr] = ModelFit_GSCSeismicHazard_QLET(F,H,model_option)

    % Probability model fitting is based on the least squares method
    
    % Candidate models
    % 1) Lognormal
    % 2) Gumbel
    % 3) Frechet
    % 4) Weibull
    
    if model_option == 1 % Lognormal distribution
        
        % z = norminv(F) and x = log(y)
        % z = -para(1)/para(2) + 1/para(2)*x
        % coef(1) = -para(1)/para(2)
        % coef(2) = 1/para(2)
        [coef,~,~,~,corr] = regress(norminv(F),[ones(length(F),1) log(H)]); 
        para(1) = -coef(1)/coef(2);
        para(2) = 1/coef(2);
    
    elseif model_option == 2 % Gumbel distribution
        
        % z = -log(-log(F)) and x = x
        % z = -para(1)*para(2) + para(2)*x
        % coef(1) = -para(1)*para(2)
        % coef(2) = para(2)
        [coef,~,~,~,corr] = regress(-log(-log(F)),[ones(length(F),1) H]); 
        para(1) = -coef(1)/coef(2);
        para(2) = coef(2);
    
    elseif model_option == 3 % Frechet distribution
        
        % z = -log(-log(F)) and x = log(x)
        % z = -log(para(1))*para(2) + para(2)*x
        % coef(1) = -log(para(1))*para(2)
        % coef(2) = para(2)
        [coef,~,~,~,corr] = regress(-log(-log(F)),[ones(length(F),1) log(H)]); 
        para(1) = exp(-coef(1)/coef(2));
        para(2) = coef(2);
    
    elseif model_option == 4 % Weibull distribution
        
        % z = log(-log(1-F)) and x = log(x)
        % z = -log(para(1))*para(2) + para(2)*x
        % coef(1) = -log(para(1))*para(2)
        % coef(2) = para(2)
        [coef,~,~,~,corr] = regress(log(-log(1-F)),[ones(length(F),1) log(H)]); 
        para(1) = exp(-coef(1)/coef(2));
        para(2) = coef(2);
    
    end

end


