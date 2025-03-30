clear
close all

% Initialize values for user inputs (these can be modified as needed)
BuildingOccupancy       = 'RES1'; % Example: Selected building occupancy % Building Occupancy RES1,RES2,RES3A,RES3B,RES3C,RES3D,RES3E,RES3F,RES4,RES5,RES6,COM1,COM2,COM3,COM4A,COM4B,COM4C,COM5,COM6,COM7A,COM7B,COM7C,COM8,COM9,COM10,IND1,IND2,IND3,IND4,IND5,IND6,AGR1,REL1,GOV1,GOV2,EDU1,EDU2
BuildingType            = 'C2L'; % Example: Selected building type % Bilding Type W1,W2,W3,W4,C1H,C1M,C1L,C2H,C2M,C2L,C3H,C3M,C3L,MH,PC1,PC2H,PC2M,PC1L,RM1M,RM1L,RM2H,RM2M,RM2L,S1H,S1M,S1L,S2H,S2M,S2L,S3,S4H,S4M,S4L,URMM,URML
DesignLevel             = 'PC'; % Example: Selected design level % Design Level PC,LC,MC,HC 
ChooseArbitraryLocation = false; % Example: True (for arbitrary location), False (for selected site)
Longitude               = -123.3656; % Example: Arbitrary longitude
Latitude                = 48.4284; % Example: Arbitrary latitude
site_of_interest        = 'Montreal'; % Victoria, Vancouver, Calgary, Toronto, Ottawa, Montreal, Quebec, LaMalbaie
Vs30msPick              = 'Arbitrary'; % 'Global   ', 'Arbitrary', 'Regional ' % it should be 9 characters
Vs30msValue             = 360; % Arbitrary Vs30
StructuralValueCAD      = 120000; % Example: User input for structural value in CAD
NonStructuralValueCAD   = 360000; % Example: User input for non-structural value in CAD
ContentValueCAD         = 220000; % Example: User input for content value in CAD
NumAnnualMaxima         = 1000000; % Number of annual maxima

load('coastline_Canada_GEODAS.mat');
load('Vulnerability_NRC2019.mat');
load('fsauid_data.mat','Exposure'); 

disp('                                                                                             ');
disp('--- Seismic Risk Analysis for Canadian Buildings Using OpenQuake-Canada Model Information ---');
disp('---                                 Site-Specific Option                                  ---');
disp('                                                                                             ');

rng(0)

%% Exposure information

% Initialize exposure data
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

Exposure.structural    = StructuralValueCAD;
Exposure.nonstructural = NonStructuralValueCAD;
Exposure.contents      = ContentValueCAD;
Exposure.eqOccType     = {BuildingOccupancy};
Exposure.EqBldgType    = {BuildingType};
Exposure.EqDesLev      = {DesignLevel};
Exposure.taxonomy{1}   = char({[char(Exposure.eqOccType),'-',char(Exposure.EqBldgType),'-',char(Exposure.EqDesLev)]});

if ~ChooseArbitraryLocation 

    sites_of_interest = site_of_interest;

    load('fsauid_data.mat',['fsauid_',char(sites_of_interest)]);
    fsauid_data = eval(['fsauid_',char(sites_of_interest)]);
    
    Exposure.csdname{1} = char(sites_of_interest);
    Exposure.fsauid{1}  = fsauid_data.fsauid{1};
    Exposure.lon(1)     = fsauid_data.lon(1);
    Exposure.lat(1)     = fsauid_data.lat(1);

else

    sites_of_interest = {'Arbitrary Location'}

    load('fsauid_summary.mat');
    fsauid_data = eval('fsauid_summary');

    arbitrary_lcoation_lon = Longitude;
    arbitrary_location_lat = Latitude;

    [dist_check,closest_rows] = min(deg2km(distance(arbitrary_location_lat,arbitrary_lcoation_lon,fsauid_summary.AverageLat,fsauid_summary.AverageLon)));
    closest_fsauid = fsauid_summary.fsauid(closest_rows);

    Exposure.csdname{1} = 'Arbitrary Location';
    Exposure.fsauid{1}  = fsauid_summary.fsauid{closest_rows};
    Exposure.lon(1)     = arbitrary_lcoation_lon;
    Exposure.lat(1)     = arbitrary_location_lat;

end

%% Seismic hazard and vulnerability information

load GSC2020_HazardMap_Total_Classified_RegionalVsAdded 

F_GSC = 1 - PE(1,:)'; 

% Filter the Site variable based on these boundaries
filtered_Site = Site(Site(:,1) >= (min(Exposure.lat)-1) & Site(:,1) <= (max(Exposure.lat)+1) & Site(:,2) >= (min(Exposure.lon)-1) & Site(:,2) <= (max(Exposure.lon)+1), :);

% Find the indices of the filtered rows in the original Site matrix
[~,filtered_Indices] = ismember(filtered_Site,Site,'rows');   

% Use these indices to filter Hazard_GSC2020_Global and Hazard_GSC2020_Regional
filtered_Hazard_GSC2020_Global   = Hazard_GSC2020_Global(filtered_Indices,:,:);
filtered_Hazard_GSC2020_Regional = Hazard_GSC2020_Regional(filtered_Indices,:,:);
filtered_Hazard_GSC2020_140      = Hazard_GSC2020_140(filtered_Indices,:,:);
filtered_Hazard_GSC2020_160      = Hazard_GSC2020_160(filtered_Indices,:,:);
filtered_Hazard_GSC2020_180      = Hazard_GSC2020_180(filtered_Indices,:,:);
filtered_Hazard_GSC2020_250      = Hazard_GSC2020_250(filtered_Indices,:,:);
filtered_Hazard_GSC2020_300      = Hazard_GSC2020_300(filtered_Indices,:,:);
filtered_Hazard_GSC2020_360      = Hazard_GSC2020_360(filtered_Indices,:,:);
filtered_Hazard_GSC2020_450      = Hazard_GSC2020_450(filtered_Indices,:,:);
filtered_Hazard_GSC2020_580      = Hazard_GSC2020_580(filtered_Indices,:,:);
filtered_Hazard_GSC2020_760      = Hazard_GSC2020_760(filtered_Indices,:,:);
filtered_Hazard_GSC2020_910      = Hazard_GSC2020_910(filtered_Indices,:,:);
filtered_Hazard_GSC2020_1100     = Hazard_GSC2020_1100(filtered_Indices,:,:);
filtered_Hazard_GSC2020_1500     = Hazard_GSC2020_1500(filtered_Indices,:,:);
filtered_Hazard_GSC2020_1600     = Hazard_GSC2020_1600(filtered_Indices,:,:);
filtered_Hazard_GSC2020_2000     = Hazard_GSC2020_2000(filtered_Indices,:,:);
filtered_Hazard_GSC2020_3000     = Hazard_GSC2020_3000(filtered_Indices,:,:);

% Identify the vulnerability function
if contains(string(Exposure.taxonomy),'W3') == 1
    taxonomy_tmp = char(string(Exposure.taxonomy));
    taxonomy_tmp(length(taxonomy_tmp)-3) = '2';
    vul_pick = find(contains(string(VulClassName(:,1)),taxonomy_tmp) == 1);
    clear taxonomy_tmp
elseif contains(string(Exposure.taxonomy),'W4') == 1 
    taxonomy_tmp = char(string(Exposure.taxonomy));
    taxonomy_tmp(length(taxonomy_tmp)-3) = '1';
    vul_pick = find(contains(string(VulClassName(:,1)),taxonomy_tmp) == 1);
    clear taxonomy_tmp
else
    vul_pick = find(contains(string(VulClassName(:,1)),string(Exposure.taxonomy)) == 1);
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
    
% Extract seismic hazard information from the GSC2020 database
[~,gscsite_id] = min(sqrt((Site(:,2)-Exposure.lon).^2 + (Site(:,1)-Exposure.lat).^2));

% If one of cities chosen, Exposure coordinates turned to closest GSC site
if ~ChooseArbitraryLocation
    Exposure.lon = Site(gscsite_id,2);
    Exposure.lat = Site(gscsite_id,1);
end

disp(['Longitude: ',num2str(Exposure.lon)]);
disp(['Latitude : ',num2str(Exposure.lat)]);

% Acquire global or regional Vs30 by interpolation
if Vs30msPick == 'Global   '
    if ~ChooseArbitraryLocation
        im_gsc(:,1) = Hazard_GSC2020_Global(gscsite_id,:,t_pick);
    else
        for i = 1:size(filtered_Hazard_GSC2020_Global,2)
            im_gsc(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_Global(:,i,t_pick),Exposure.lon,Exposure.lat);
        end
    end
elseif Vs30msPick == 'Regional '
    if ~ChooseArbitraryLocation
        if ~isnan(Hazard_GSC2020_Regional(gscsite_id,1,t_pick))
            im_gsc(:,1) = Hazard_GSC2020_Regional(gscsite_id,:,t_pick);
        else
            disp('No regional Vs30 found for selected site');
            return
        end
    else
        if ~isnan(Hazard_GSC2020_Regional(gscsite_id,1,t_pick))
            for i = 1:size(filtered_Hazard_GSC2020_Regional,2)
                im_gsc(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_Regional(:,i,t_pick),Exposure.lon,Exposure.lat);
            end
        else
            disp('No regional Vs30 found for selected site');
            return
        end
    end                 
elseif Vs30msPick == 'Arbitrary'
    if Vs30msValue >= 140 && Vs30msValue < 160
        for i = 1:size(filtered_Hazard_GSC2020_140,2)            
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_140(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_160(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(140))/...
                (log10(160) - log10(140))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 160 && Vs30msValue < 180 
        for i = 1:size(filtered_Hazard_GSC2020_160,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_160(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_180(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(160))/...
                (log10(180) - log10(160))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 180 && Vs30msValue < 250 
        for i = 1:size(filtered_Hazard_GSC2020_180,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_180(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_250(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(180))/...
                (log10(250) - log10(180))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 250 && Vs30msValue < 300
        for i = 1:size(filtered_Hazard_GSC2020_250,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_250(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_300(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(250))/...
                (log10(300) - log10(250))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 300 && Vs30msValue < 360
        for i = 1:size(filtered_Hazard_GSC2020_300,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_300(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_360(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(300))/...
                (log10(360) - log10(300))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 360 && Vs30msValue < 450
        for i = 1:size(filtered_Hazard_GSC2020_360,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_360(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_450(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(360))/...
                (log10(450) - log10(360))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 450 && Vs30msValue < 580
        for i = 1:size(filtered_Hazard_GSC2020_450,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_450(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_580(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(450))/...
                (log10(580) - log10(450))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 580 && Vs30msValue < 760
        for i = 1:size(filtered_Hazard_GSC2020_580,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_580(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_760(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(580))/...
                (log10(760) - log10(580))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 760 && Vs30msValue < 910
        for i = 1:size(filtered_Hazard_GSC2020_760,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_760(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_910(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(760))/...
                (log10(910) - log10(760))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 910 && Vs30msValue < 1100
        for i = 1:size(filtered_Hazard_GSC2020_910,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_910(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_1100(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(910))/...
                (log10(1100) - log10(910))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 1100 && Vs30msValue < 1500
        for i = 1:size(filtered_Hazard_GSC2020_1100,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_1100(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_1500(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(1100))/...
                (log10(1500) - log10(1100))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 1500 && Vs30msValue < 1600
        for i = 1:size(filtered_Hazard_GSC2020_1500,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_1500(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_1600(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(1500))/...
                (log10(1600) - log10(1500))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 1600 && Vs30msValue < 2000
        for i = 1:size(filtered_Hazard_GSC2020_1600,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_1600(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_2000(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(1600))/...
                (log10(2000) - log10(1600))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    elseif Vs30msValue >= 2000 && Vs30msValue <= 3000
        for i = 1:size(filtered_Hazard_GSC2020_2000,2)
            temp_im_min(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_2000(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_max(i,1) = griddata(filtered_Site(:,2),filtered_Site(:,1),filtered_Hazard_GSC2020_3000(:,i,t_pick),Exposure.lon,Exposure.lat);
            temp_im_gsc(i,1) = log10(temp_im_min(i,1)) + ((log10(temp_im_max(i,1)) - log10(temp_im_min(i,1)))*...
                ((log10(Vs30msValue) - log10(2000))/...
                (log10(3000) - log10(2000))));
            temp_im_gsc(i,1) = 10^temp_im_gsc(i,1);
            im_gsc(i,1) = temp_im_gsc(i,1);
        end
    end
end
    
if t_pick == 11; im_gsc = 100*im_gsc; end % PGV (cm/s)

% Tail approximation of seismic hazard curve: 1) Lognormal, 2) Gumbel, 3) Frechet, 4) Weibull
RRtmp = zeros(4,4);
for jj = 1:4
    [~,~,RRtmp(jj,:)] = ModelFit_GSCSeismicHazard_QLE(F_GSC,im_gsc,jj); 
end
[~,model_pick] = max(RRtmp(:,1));
[para_im,coef_im,RR_im] = ModelFit_GSCSeismicHazard_QLE(F_GSC,im_gsc,model_pick);
    
% Summary hazard-vulnerability information
% 1    ) Vulnerability function ID (as in Vulnerability_NRC2019)
% 2    ) Intensity measure ID (as in GSC2020_HazardMap)
% 3    ) GSC hazard location ID (as in GSC2020_HazardMap)
% 4 -13) GSC hazard values for 10 return periods (as in GSC2020_HazardMap)
% 14   ) Probability distribution type for tail approximation 
% 15-16) Parameters for tail approximation 
% 17-18) Coefficients for tail approximation 
% 19-22) Regression output for tail approximation (R^2, F, p-value, error variance)
VulnerabilityHazard = [vul_pick t_pick gscsite_id im_gsc(:,1)' model_pick para_im coef_im' RR_im];

%% Seismic vulnerability functions

im_risk_min = 0.01; % (g)

IM_vulfunc = 0.01:0.01:8.0; 

% Interpolation or lognormal approximation of vulnerability functions
if VulIM_ST(vul_pick,VulNumData(vul_pick,1)) >= 8
    VulMU_ST_approx(:,1) = interp1q(VulIM_ST(vul_pick,1:VulNumData(vul_pick,1))',VulMU_ST(vul_pick,1:VulNumData(vul_pick,1))',IM_vulfunc');    
    VulMU_ST_approx(find(IM_vulfunc<VulIM_ST(vul_pick,1)),1) = 0;
else
    [datatmp1(:,1),datatmp2] = unique(VulMU_ST(vul_pick,1:VulNumData(vul_pick,1))); 
    datatmp3(:,1) = VulIM_ST(vul_pick,datatmp2);
    vul_para_approx(1,1) = interp1q(datatmp1,datatmp3,0.50);
    if isnan(vul_para_approx(1,1)) == 1
        point1 = [datatmp3(find(datatmp1>0.25,1,'first')) datatmp1(find(datatmp1>0.25,1,'first'))];
        point2 = [datatmp3(end) datatmp1(end)];
        vul_para_approx(1,1) = point1(1) + (0.5-point1(2))/(point2(2)-point1(2))*(point2(1)-point1(1));
        clear point1 point2
    end
    vul_para_approx(1,2) = interp1q(datatmp1,datatmp3,0.16);
    vul_para_approx(1,3) = interp1q(datatmp1,datatmp3,0.84);
    if isnan(vul_para_approx(1,3)) == 1
        vul_para_approx(1,4) = log(vul_para_approx(1,1)) - log(vul_para_approx(1,2));
    else
        vul_para_approx(1,4) = (log(vul_para_approx(1,3)) - log(vul_para_approx(1,2)))/2;
    end
    VulMU_ST_approx(:,1) = normcdf(log(IM_vulfunc/vul_para_approx(1,1))/vul_para_approx(1,4));    
    clear datatmp1 datatmp2 datatmp3
end
    
if VulIM_NS(vul_pick,VulNumData(vul_pick,2)) >= 8
    VulMU_NS_approx(:,1) = interp1q(VulIM_NS(vul_pick,1:VulNumData(vul_pick,2))',VulMU_NS(vul_pick,1:VulNumData(vul_pick,2))',IM_vulfunc');    
    VulMU_NS_approx(find(IM_vulfunc<VulIM_NS(vul_pick,1)),1) = 0;
else
    [datatmp1(:,1),datatmp2] = unique(VulMU_NS(vul_pick,1:VulNumData(vul_pick,2))); 
    datatmp3(:,1) = VulIM_NS(vul_pick,datatmp2);
    vul_para_approx(1,5) = interp1q(datatmp1,datatmp3,0.50);
    if isnan(vul_para_approx(1,5)) == 1
        point1 = [datatmp3(find(datatmp1>0.25,1,'first')) datatmp1(find(datatmp1>0.25,1,'first'))];
        point2 = [datatmp3(end) datatmp1(end)];
        vul_para_approx(1,5) = point1(1) + (0.5-point1(2))/(point2(2)-point1(2))*(point2(1)-point1(1));
        clear point1 point2
    end
    vul_para_approx(1,6) = interp1q(datatmp1,datatmp3,0.16);
    vul_para_approx(1,7) = interp1q(datatmp1,datatmp3,0.84);
    if isnan(vul_para_approx(1,7)) == 1
        vul_para_approx(1,8) = log(vul_para_approx(1,5)) - log(vul_para_approx(1,6));
    else
        vul_para_approx(1,8) = (log(vul_para_approx(1,7)) - log(vul_para_approx(1,6)))/2;
    end
    VulMU_NS_approx(:,1) = normcdf(log(IM_vulfunc/vul_para_approx(1,5))/vul_para_approx(1,8));    
    clear datatmp1 datatmp2 datatmp3
end
    
if VulIM_CO(vul_pick,VulNumData(vul_pick,3)) >= 8
    VulMU_CO_approx(:,1) = interp1q(VulIM_CO(vul_pick,1:VulNumData(vul_pick,3))',VulMU_CO(vul_pick,1:VulNumData(vul_pick,3))',IM_vulfunc');    
    VulMU_CO_approx(find(IM_vulfunc<VulIM_CO(vul_pick,1)),1) = 0;
else
    [datatmp1(:,1),datatmp2] = unique(VulMU_CO(vul_pick,1:VulNumData(vul_pick,3))); 
    datatmp3(:,1) = VulIM_CO(vul_pick,datatmp2);
    vul_para_approx(1,9)  = interp1q(datatmp1,datatmp3,0.50);
    if isnan(vul_para_approx(1,9)) == 1
        point1 = [datatmp3(find(datatmp1>0.25,1,'first')) datatmp1(find(datatmp1>0.25,1,'first'))];
        point2 = [datatmp3(end) datatmp1(end)];
        vul_para_approx(1,9) = point1(1) + (0.5-point1(2))/(point2(2)-point1(2))*(point2(1)-point1(1));
        clear point1 point2
    end
    vul_para_approx(1,10) = interp1q(datatmp1,datatmp3,0.16);
    vul_para_approx(1,11) = interp1q(datatmp1,datatmp3,0.84);
    if isnan(vul_para_approx(1,11)) == 1
        vul_para_approx(1,12) = log(vul_para_approx(1,9)) - log(vul_para_approx(1,10));
    else
        vul_para_approx(1,12) = (log(vul_para_approx(1,11)) - log(vul_para_approx(1,10)))/2;
    end
    VulMU_CO_approx(:,1) = normcdf(log(IM_vulfunc/vul_para_approx(1,9))/vul_para_approx(1,12));    
    clear datatmp1 datatmp2 datatmp3
end
    
%% Seismic risk analysis

Psimu = (1:NumAnnualMaxima)/(NumAnnualMaxima+1);

Frange = [0.9 0.95 0.98 0.99 0.995 0.998 0.999 0.9995 0.9998 0.9999];

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
    a_beta = (1-(0.01:0.01:0.99))/CoV_vul_beta^2-(0.01:0.01:0.99);
    b_beta = ((1-(0.01:0.01:0.99))/CoV_vul_beta^2-(0.01:0.01:0.99)).*((1-(0.01:0.01:0.99))./(0.01:0.01:0.99));
    figure (99); plot(0.01:0.01:0.99,a_beta,'b-',0.01:0.01:0.99,b_beta,'r-');
    meanloss_upper = 1/(1+CoV_vul_beta^2);
end

if VulnerabilityHazard(1,14) == 1
    ur_risk_min = logncdf(im_risk_min,VulnerabilityHazard(1,15),VulnerabilityHazard(1,16));
    ur_risk = ur_risk_min + rand(NumAnnualMaxima,1)*(1-ur_risk_min);
    im_risk = exp(VulnerabilityHazard(1,15)+VulnerabilityHazard(1,16)*norminv(ur_risk));
elseif VulnerabilityHazard(1,14) == 2
    ur_risk_min = exp(-exp(-(im_risk_min-VulnerabilityHazard(1,15))*VulnerabilityHazard(1,16)));
    ur_risk = ur_risk_min + rand(NumAnnualMaxima,1)*(1-ur_risk_min);
    im_risk = VulnerabilityHazard(1,15)-log(-log(ur_risk))/VulnerabilityHazard(1,16);
elseif VulnerabilityHazard(1,14) == 3
    ur_risk_min = exp(-(VulnerabilityHazard(1,15)/im_risk_min)^VulnerabilityHazard(1,16));
    ur_risk = ur_risk_min + rand(NumAnnualMaxima,1)*(1-ur_risk_min);
    im_risk = exp(log(VulnerabilityHazard(1,15))-log(-log(ur_risk))/VulnerabilityHazard(1,16));
elseif VulnerabilityHazard(1,14) == 4
    ur_risk_min = 1-exp(-(im_risk_min/VulnerabilityHazard(1,15))^VulnerabilityHazard(1,16));
    ur_risk = ur_risk_min + rand(NumAnnualMaxima,1)*(1-ur_risk_min);
    im_risk = exp(log(VulnerabilityHazard(1,15))+log(-log(1-ur_risk))/VulnerabilityHazard(1,16));
end
       
meanloss_st = interp1q(IM_vulfunc',VulMU_ST_approx(:,1),min(max(im_risk,im_risk_min),8));
meanloss_ns = interp1q(IM_vulfunc',VulMU_NS_approx(:,1),min(max(im_risk,im_risk_min),8));
meanloss_co = interp1q(IM_vulfunc',VulMU_CO_approx(:,1),min(max(im_risk,im_risk_min),8));
    
if option_vulnerability_uncertainty == 0
        
    loss_st = meanloss_st;
    loss_ns = meanloss_ns;
    loss_co = meanloss_co;
      
elseif option_vulnerability_uncertainty == 1
        
    loss_st = lognrnd(log(meanloss_st/sqrt(1+CoV_vul_logn^2)),sqrt(log(1+CoV_vul_logn^2)),NumAnnualMaxima,1);
    loss_ns = lognrnd(log(meanloss_ns/sqrt(1+CoV_vul_logn^2)),sqrt(log(1+CoV_vul_logn^2)),NumAnnualMaxima,1);
    loss_co = lognrnd(log(meanloss_co/sqrt(1+CoV_vul_logn^2)),sqrt(log(1+CoV_vul_logn^2)),NumAnnualMaxima,1);
        
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
for jj = 1:length(Loss_RP)
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
    loss_prob_thresh(jj,1) = (1-ur_risk_min)/NumAnnualMaxima*length(find(LOSS_all_risk>=Loss_thresh(jj)));
end

% Soil condition Global
for i = 1:size(Site_Vs30,1)
    if Site_Vs30(i,7) == 1
        Site_Vs30_class(i,1) = 'A';
    elseif Site_Vs30(i,7) == 2
        Site_Vs30_class(i,1) = 'B';
    elseif Site_Vs30(i,7) == 3
        Site_Vs30_class(i,1) = 'C';
    elseif Site_Vs30(i,7) == 4
        Site_Vs30_class(i,1) = 'D';
    elseif Site_Vs30(i,7) == 5
        Site_Vs30_class(i,1) = 'E';
    end
end
% Soil condition Regional
for i = 1:size(Site_Vs30,1)
    if Site_Vs30(i,8) == 1
        Site_Vs30_class(i,2) = 'A';
    elseif Site_Vs30(i,8) == 2
        Site_Vs30_class(i,2) = 'B';
    elseif Site_Vs30(i,8) == 3
        Site_Vs30_class(i,2) = 'C';
    elseif Site_Vs30(i,8) == 4
        Site_Vs30_class(i,2) = 'D';
    elseif Site_Vs30(i,8) == 5
        Site_Vs30_class(i,2) = 'E';
    end
end
           
% Define Vs30 input
if Vs30msPick == 'Global   '

    Vs30_site_temp = round(Site_Vs30(VulnerabilityHazard(1,3),3));
    Vs30_siteclass_temp = Site_Vs30_class(VulnerabilityHazard(1,3),1);

elseif Vs30msPick == 'Regional '

    Vs30_site_temp = round(Site_Vs30(VulnerabilityHazard(1,3),4));
    Vs30_siteclass_temp = Site_Vs30_class(VulnerabilityHazard(1,3),2);

elseif Vs30msPick == 'Arbitrary'

    Vs30_site_temp = round(Vs30msValue);
    % Site classification - 1:A 2:B 3:C 4:D 5:E 
    if Vs30_site_temp > 1500
        Vs30_siteclass_temp = 'A';
    elseif Vs30_site_temp > 760 && Vs30_site_temp <= 1500 
        Vs30_siteclass_temp = 'B';
    elseif Vs30_site_temp > 360 && Vs30_site_temp <= 760 
        Vs30_siteclass_temp = 'C';
    elseif Vs30_site_temp > 180 && Vs30_site_temp <= 360
        Vs30_siteclass_temp = 'D';
    elseif Vs30_site_temp > 140 && Vs30_site_temp <= 180
        Vs30_siteclass_temp = 'E';
    elseif Vs30_site_temp <= 140
        Vs30_siteclass_temp = 'F';
    end

end

disp(['Vs30      : ',num2str(Vs30_site_temp),'m/s']);
disp(['Site class: ',Vs30_siteclass_temp]);

figure
plot(VulnerabilityHazard(1,4:13),1-F_GSC,'bo',sort(im_risk),1-sort(ur_risk),'r.');
title([sites_of_interest,' - Number of simulations: ',num2str(NumAnnualMaxima)]);
xlabel(string(VulIMName(vul_pick,1))); ylabel('Annual exceedance frequency'); axis square;
ax = gca; ax.Box = 'on'; set(gca,'XScale','log','YScale','log');
ax.XLim = [0.001 10]; ax.XScale = 'log'; ax.XMinorTick = 'on'; ax.XGrid = 'on'; ax.XMinorGrid = 'on';
ax.YLim = [1e-5 0.1]; ax.YScale = 'log'; ax.YMinorTick = 'on'; ax.YGrid = 'on'; ax.YMinorGrid = 'on';

figure
plot( ...
    VulIM_ST(vul_pick,1:VulNumData(vul_pick,1)),VulMU_ST(vul_pick,1:VulNumData(vul_pick,1)),'bo',...
    VulIM_NS(vul_pick,1:VulNumData(vul_pick,2)),VulMU_NS(vul_pick,1:VulNumData(vul_pick,2)),'rs', ...
    VulIM_CO(vul_pick,1:VulNumData(vul_pick,3)),VulMU_CO(vul_pick,1:VulNumData(vul_pick,3)),'g^',...
    IM_vulfunc,VulMU_ST_approx,'b-',IM_vulfunc,VulMU_NS_approx,'r-',IM_vulfunc,VulMU_CO_approx,'g-','LineWidth',1);
xlabel(string(VulIMName(vul_pick,1))); ylabel('Loss ratio'); axis square;
title(['Vulnerability function - ',char(string(VulClassName(vul_pick,1)))])
legend('Structural data', 'NonStructural data', 'Content data', 'Structural fit', 'NonStructural fit', 'Content fit','Location','southeast'); 
ax = gca; ax.XLim = [0 8]; ax.YLim = [0 1];

figure
plot((1/10^6).*LOSS_all_risk,loss_ex_prob,'b-');
title([sites_of_interest,' - Number of simulations: ',num2str(NumAnnualMaxima)]);
xlabel('Loss (million CAD)'); ylabel('Annual exceedance frequency'); axis square;
ax = gca; ax.Box = 'on'; ax.XGrid = 'on';
ax.YLim = [1e-05 0.1]; ax.YScale = 'log'; ax.YMinorTick = 'on'; ax.YMinorGrid = 'on'; ax.YGrid = 'on';

loss_frac = round(loss_frac);
disp(['  100-year VaR loss: ',num2str(loss_frac(1))]);
disp(['  200-year VaR loss: ',num2str(loss_frac(2))]);
disp(['  500-year VaR loss: ',num2str(loss_frac(2))]);
disp([' 1000-year VaR loss: ',num2str(loss_frac(4))]);
disp([' 2500-year VaR loss: ',num2str(loss_frac(5))]);
disp([' 5000-year VaR loss: ',num2str(loss_frac(6))]);
disp(['10000-year VaR loss: ',num2str(loss_frac(7))]);

disp(['Selected Longitude : ',num2str(Exposure.lon)]);
disp(['Selected Latitude  : ',num2str(Exposure.lat)]);
disp(['Selected Vs30 (m/s): ',num2str(Vs30_site_temp)]);

disp(['Global Vs30 (m/s)  : ',num2str(round(Site_Vs30(VulnerabilityHazard(1,3),3)))]);
disp(['Regional Vs30 (m/s): ',num2str(round(Site_Vs30(VulnerabilityHazard(1,3),4)))]);

disp('          ');

%% M-function for fitting probabilistic models to the GSC's hazard values

function [para,coef,corr] = ModelFit_GSCSeismicHazard_QLE(F,H,model_option)

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

