function STORM(MODE,NUMSIMS,NUMSIMYRS,SEASONS,PTOT_SCENARIO,STORMINESS_SCENARIO,PTOT_SCENARIO2,STORMINESS_SCENARIO2,ET_SCENARIO)
%
%STORM - STOchastic Rainfall Model STORM produces realistic regional or watershed
%rainfall under various climate scenarios based on empirical-stochastic
%selection of historical rainfall characteristics.
%
%Based on Singer, M. B., and Michaelides, K. (2017), Deciphering the expression
%of climate change within the Lower Colorado River basin by stochastic simulation
%of convective rainfall.
%
%version name: STORM_v1
%
%
%Author: Michael Singer 2017 Date created: 2015-6 Last modified 03/23/2018
%----------------------------------------------------------------------------
%
%USAGE: STORM(MODE,NUMSIMS,NUMSIMYRS,SEASONS,PTOT_SCENARIO,STORMINESS_SCENARIO,PTOT_SCENARIO2,STORMINESS_SCENARIO2,ET_SCENARIO)
%
%MODE is a string in single quotations containing one of these two options:
%'Validation' for validating STORM's output against historical observations
%'Simulation' for simulating stochastic rainfall under various climate scenarios
%MODE affects where STORM generates output--either for a set of
%gauge locations for comparison with observational data (Validation) or on
%an aribitrary watershed grid (Simulation).
%
%NUMSIMS is the integer number of x-year (NUMSIMYRS) simulations to be run
%as a batch. Each simulation will produce its own output folder that is
%timestamped with relevant simulation information
%
%NUMSIMYRS is the integer number of years in each simulation. Recommended
%value >=30
%
%SEASONS is number specifying whether there are 'one' or 'two' seasons. If there are
%two seasons, pdfs will be required specifying rainfall characteristics for
%the 2nd season (Ptot, Duration, etc). Also, ptot and storminess scenarios will be required
%arguments in the function (see below). Currently, season one runs from May 1 to
%Sept 30 and season two runs from Oct 1 to Apr 30, but these dates can be
%changed.
%
%PTOT_SCENARIO is a string in single quotations specifying the type of Ptot
%simulation to be run. This variable characterizes the degree of overall
%wetness as measured by annual precipitation totals. There are currently 5
%options here:
%'ptotC' for control (stationary) conditions akin to observations of precipitation totals
%'ptot+' for step change of increased overall wetness from observed
%'ptot-' for step change of decreased overall wetness from observed
%'ptotT+' for progressive increasing trend in storminess from observed
%'ptotT-' for progressive decreasing trend in storminess from observed
%
%STORMINESS_SCENARIO is a string in single quotations specifying the type of
%Storminess simulation to be run. This variable characterizes the degree of
%overall storminess as measured by storm intensity for a particular
%duration. There are currently 5 options here:
%'stormsC' for no change in storminess -- use intensity-duration relationships from historical records
%'storms+' for step change of increased storminess from observed
%'storms-' for step change of decreased storminess from observed
%'stormsT+' for progressive increasing trend in storminess from observed
%'stormsT-' for progressive decreasing trend in storminess from observed.
%
%%PTOT_SCENARIO2 is a string in single quotations specifying the type of Ptot
%simulation to be run. This variable characterizes the degree of overall
%wetness as measured by annual precipitation totals. There are currently 5
%options here:
%'ptot2C' for control (stationary) conditions akin to observations of precipitation totals
%'ptot2+' for step change of increased overall wetness from observed
%'ptot2-' for step change of decreased overall wetness from observed
%'ptot2T+' for progressive increasing trend in storminess from observed
%'ptot2T-' for progressive decreasing trend in storminess from observed
%For SEASON = one, use [] for this argument.
%
%STORMINESS_SCENARIO2 is a string in single quotations specifying the type of
%Storminess simulation to be run. This variable characterizes the degree of
%overall storminess as measured by storm intensity for a particular
%duration. There are currently 5 options here:
%'storms2C' for no change in storminess -- use intensity-duration relationships from historical records
%'storms2+' for step change of increased storminess from observed
%'storms2-' for step change of decreased storminess from observed
%'storms2T+' for progressive increasing trend in storminess from observed
%'storms2T-' for progressive decreasing trend in storminess from observed.
%For SEASON = one, use [] for this argument.
%
%
%ET_SCENARIO is a string in single quotations specifying the type of
%Storminess simulation to be run. This variable characterizes the degree of
%evaporative demand and its reflection of climate change. There are currently 5 options here:
%'ETC' for no change in ET - use historical records
%'ET+' for step change of increased ET from observed
%'ET-' for step change of decreased ET from observed
%Note: TRENDS and SEASONAL CHANGES IN ET HAVE NOT YET BEEN CODED!!
%'ET_T+' for progressive increasing trend in ET from observed
%'ET_T-' for progressive decreasing trend in ET from observed.
%
%
%The code requires the following Matlab variables (.mat files) to be
%created and placed in a folder called 'model_input':
%
%%Pdfs listed below should be generated by the makedist function after exploring the appropriate distribution type in the dfittool
%'Ptot_pdf' - A pdf fitted to all available station precip totals data (mm).
%'Ptot_pdf2' - A pdf fitted to all available station precip totals data for season 2 (mm).
%'Ptot_pdf_cc_up' Historical pdf shifted toward wetter conditions
%'Ptot_pdf_cc_up2' Historical pdf for season 2 shifted toward wetter conditions
%'Ptot_pdf_cc_down' Historical pdf shifted toward drier conditions
%'Ptot_pdf_cc_down2' Historical pdf for season 2 shifted toward drier conditions
%'Duration_pdf' A pdf fitted to all available storm duration data (min)
%'Duration_pdf2' A pdf fitted to all available storm duration data of season 2 (min)
%'Area_pdf' A pdf fitted to all available storm area data (m^2)
%'Area_pdf2' A pdf fitted to all available storm area data for season 2 (m^2)
%'Int_arr_pdf' A pdf fitted to all available interarrival times data (hrs)
%'Int_arr_pdf2' A pdf fitted to all available interarrival times data for season 2 (hrs)
%'Recess_pdf' A pdf of recession coefficients of intensity with distance from storm center (mm/hr/km).
%
%'Easting' Actual rain gauge location. Used when MODE = 'Validation'
%'Northing' Actual rain gauge location. Used when MODE = 'Validation'
%'gauges' Numbers of actual rain gauges. Used when MODE = 'Validation'
%'gauge_elev' Elevations of actual rain gauges. Elevation value is extracted from field called 'RASTERVALU'. Used when  MODE = 'Validation'
%'boundary/boundary.shp' Shapefile of the watershed boundary
%'point_elevations.shp' Shapefile of the elevations at each 'gauging' point on the grid within the watershed boundary
%'X' Longitudinal data for each grid point. Used when MODE = 'Simulation'
%'Y' Latitudinal data for each grid point. Used when MODE = 'Simulation'. X & Y will be sampled below to determine storm center location.
%'fuzz' A vector of fuzzy tolerance values for intensity selection, assumed based on differences between PI-PD curves(mm/hr).
%'bndry_buffer.shp' A shapefile of the watershed boundary buffer (5000 m). Used for storm center determination allowing for partial storm coverage.
%'ET_monthly_day' A matrix of daytime ET values organized with one month per column, drawn from data or calculations (mm).
%'ET_monthly_night' A matrix of nighttime ET values organized with one month per column, drawn from data or calculations (mm).
%
%----------------------------------------------------------------------------

%Restatement of arguments in function
XXX = ('Thank you for using STORM. Your variable settings are: ');
disp(XXX)
XXX = ['Type of Model Run = ',MODE];
disp(XXX)
XXX = ['Number of Simulations = ',num2str(NUMSIMS)];
disp(XXX)
XXX = ['Number of Simulation Years = ',num2str(NUMSIMYRS)];
disp(XXX)
if  strcmpi('two',SEASONS) == 1
    XXX = ('Number of Seasons = Two');
elseif  strcmpi('one',SEASONS) == 1
    XXX = ('Number of Seasons = One');
end
disp(XXX)
if  strcmpi('ptotC',PTOT_SCENARIO) == 1
    XXX = ('Scenario of Total Precipitation (Season 1) = Control climate');
elseif  strcmpi('ptot+',PTOT_SCENARIO) == 1
    XXX = ('Scenario of Total Precipitation (Season 1) = Step change increase in seasonal precipitation');
elseif  strcmpi('ptot-',PTOT_SCENARIO) == 1
    XXX = ('Scenario of Total Precipitation (Season 1) = Step change decrease in seasonal precipitation');
elseif  strcmpi('ptotT+',PTOT_SCENARIO) == 1
    XXX = ('Scenario of Total Precipitation (Season 1) = Postive trend in seasonal precipitation');
elseif  strcmpi('ptotT-',PTOT_SCENARIO) == 1
    XXX = ('Scenario of Total Precipitation (Season 1) = Negative trend in seasonal precipitation');
end
disp(XXX)
if  strcmpi('stormsC',STORMINESS_SCENARIO) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 1) = Control climate');
elseif  strcmpi('storms+',STORMINESS_SCENARIO) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 1) = Step change increase in seasonal storminess');
elseif  strcmpi('storms-',STORMINESS_SCENARIO) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 1) = Step change decrease in seasonal storminess');
elseif  strcmpi('stormsT+',STORMINESS_SCENARIO) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 1) = Postive trend in seasonal storminess');
elseif  strcmpi('stormsT-',STORMINESS_SCENARIO) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 1) = Negative trend in seasonal storminess');
end
disp(XXX)
if  strcmpi('ptot2C',PTOT_SCENARIO2) == 1
    XXX = ('Scenario of Total Precipitation (Season 2) = Control climate');
elseif  strcmpi('ptot2+',PTOT_SCENARIO2) == 1
    XXX = ('Scenario of Total Precipitation (Season 2) = Step change increase in seasonal precipitation');
elseif  strcmpi('ptot2-',PTOT_SCENARIO2) == 1
    XXX = ('Scenario of Total Precipitation (Season 2) = Step change decrease in seasonal precipitation');
elseif  strcmpi('ptot2T+',PTOT_SCENARIO2) == 1
    XXX = ('Scenario of Total Precipitation (Season 2) = Postive trend in seasonal precipitation');
elseif  strcmpi('ptot2T-',PTOT_SCENARIO2) == 1
    XXX = ('Scenario of Total Precipitation (Season 2) = Negative trend in seasonal precipitation');
end
disp(XXX)
if  strcmpi('storms2C',STORMINESS_SCENARIO2) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 2) = Control climate');
elseif  strcmpi('storms2+',STORMINESS_SCENARIO2) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 2) = Step change increase in seasonal storminess');
elseif  strcmpi('storms2-',STORMINESS_SCENARIO2) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 2) = Step change decrease in seasonal storminess');
elseif  strcmpi('storms2T+',STORMINESS_SCENARIO2) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 2) = Postive trend in seasonal storminess');
elseif  strcmpi('storms2T-',STORMINESS_SCENARIO2) == 1
    XXX = ('Scenario of Precipitation Intensity (Season 2) = Negative trend in seasonal storminess');
end
disp(XXX)
if  strcmpi('ETC',ET_SCENARIO) == 1
    XXX = ('Scenario of ET =  Control climate');
elseif  strcmpi('ET+',ET_SCENARIO) == 1
    XXX = ('Scenario of ET = Step change increase in PET');
elseif  strcmpi('ET-',ET_SCENARIO) == 1
    XXX = ('Scenario of ET = Step change decrease in PET');
end
disp(XXX)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: There are several hard-coded parameters below, all of which can be
%modified. They should indeed be modified for extension to any particular basin.


storminess_scaling_factor = 0.05; %This scalar specifies the fractional change in intensity per year when storm_trend is applied in STORMINESS_SCENARIO
storminess_scaling_factor2 = 0.05; %This scalar specifies the fractional change in intensity per year when storm_trend is applied in STORMINESS_SCENARIO2
storm_stepchange = 0.25; %This scalar specifies the value of fractional step change in intensity when storms+ or storms- are applied in STORMINESS_SCENARIO
storm_stepchange2 = 0.25; %This scalar specifies the value of fractional step change in intensity when storms+ or storms- are applied in STORMINESS_SCENARIO2
ptot_scaling_factor = 0.05; %This scalar specifies the fractional change in wetness per year when storm_trend is applied in PTOT_SCENARIO
ptot_scaling_factor2 = 0.05; %This scalar specifies the fractional change in wetness per year when storm_trend is applied in PTOT_SCENARIO2
ET_scaling_factor = 0.25; %This scalar specifies the value of fractional step change in ET when ET+ or ET- are applied in ET_SCENARIO
%Note: These scalars can all be changed to reflect observed or projected trends in climate

numcurves = 11; %number of intensity-duration curves (see below for curve equations)


%IINITIALIZE VARIABLES

%Add variable for number of simulations of simyears
simyears = NUMSIMYRS; %number of years to simulate

fs = filesep;
%load pdfs for season one
if  strcmpi('ptotC',PTOT_SCENARIO) == 1
    load(strcat('model_input',fs,'Ptot_pdf.mat'), 'Ptot_pdf')% This is the pdf fitted to all available station precip data (normal dist). It will be sampled below.
    mu = mean(Ptot_pdf);
    sigma = std(Ptot_pdf);
elseif  strcmpi('ptot+',PTOT_SCENARIO) == 1
    load(strcat('model_input',fs,'Ptot_pdf_cc_up.mat'), 'Ptot_pdf_cc_up')% This is the pdf fitted to all available station precip data (normal dist). It will be sampled below.
    mu = mean(Ptot_pdf_cc_up);
    sigma = std(Ptot_pdf_cc_up);
elseif strcmpi('ptot-',PTOT_SCENARIO) == 1
    load(strcat('model_input',fs,'Ptot_pdf_cc_down.mat'), 'Ptot_pdf_cc_down')% This is the pdf fitted to all available station precip data (normal dist). It will be sampled below.
    mu = mean(Ptot_pdf_cc_down);
    sigma = std(Ptot_pdf_cc_down);
elseif strcmpi('ptotT+',PTOT_SCENARIO) == 1
    load(strcat('model_input',fs,'Ptot_pdf.mat'), 'Ptot_pdf')% This is the pdf fitted to all available station precip data (normal dist). It will be modified as a trend and sampled below.
    mu = mean(Ptot_pdf);
    sigma = std(Ptot_pdf);
elseif strcmpi('ptotT-',PTOT_SCENARIO) == 1
    load(strcat('model_input',fs,'Ptot_pdf.mat'), 'Ptot_pdf')% This is the pdf fitted to all available station precip data (normal dist). It will be modified as a trend and sampled below.
    mu = mean(Ptot_pdf);
    sigma = std(Ptot_pdf);
end
load(strcat('model_input',fs,'Duration_pdf.mat'), 'Duration_pdf')% This is the pdf fitted to all available station duration data (GEV dist). It will be sampled below.
load(strcat('model_input',fs,'Area_pdf.mat'), 'Area_pdf')% This is the pdf fitted to all available station area data (EV dist). It will be sampled below.
load(strcat('model_input',fs,'Int_arr_pdf.mat'), 'Int_arr_pdf')% This is the pdf fitted to all available station interarrival time data (GEV dist). It will be sampled below.

%Load pdfs for season two
if  strcmpi('two',SEASONS) == 1
    load(strcat('model_input',fs,'Duration_pdf2.mat'), 'Duration_pdf2')% This is the pdf for season2 fitted to all available station duration data (log-logistic dist). It will be sampled below.
    load(strcat('model_input',fs,'Area_pdf2.mat'), 'Area_pdf2')% This is the pdf fitted for season2 to all available station area data (EV dist). It will be sampled below.
    load(strcat('model_input',fs,'Int_arr_pdf2.mat'), 'Int_arr_pdf2')% This is the pdf for season2 fitted to all available station interarrival time data (GEV dist). It will be sampled below.
    if  strcmpi('ptot2C',PTOT_SCENARIO2) == 1
        load(strcat('model_input',fs,'Ptot_pdf2.mat'), 'Ptot_pdf2')% This is the pdf fitted to all available station precip data (gamma dist). It will be sampled below.
        %mu2 = mean(Ptot_pdf2);
        %sigma2 = std(Ptot_pdf2);
    elseif  strcmpi('ptot2+',PTOT_SCENARIO2) == 1
        load(strcat('model_input',fs,'Ptot_pdf_cc_up2.mat'), 'Ptot_pdf_cc_up2')% This is the pdf for season2 fitted to all available station precip data (gamma dist). It will be sampled below.
        %mu2 = mean(Ptot_pdf_cc_up2);
        %sigma2 = std(Ptot_pdf_cc_up2);
    elseif strcmpi('ptot2-',PTOT_SCENARIO2) == 1
        load(strcat('model_input',fs,'Ptot_pdf_cc_down2.mat'), 'Ptot_pdf_cc_down2')% This is the pdf for season2 fitted to all available station precip data (gamma dist). It will be sampled below.
        %mu2 = mean(Ptot_pdf_cc_down2);
        %sigma2 = std(Ptot_pdf_cc_down2);
    elseif strcmpi('ptot2T+',PTOT_SCENARIO2) == 1
        load(strcat('model_input',fs,'Ptot_pdf2.mat'), 'Ptot_pdf2')% This is the pdf fitted to all available station precip data (gamma dist). It will be modified as a trend and sampled below.
        %mu2 = mean(Ptot_pdf2);
        gamma_a2 = 1.65466; %this is actually the shape parameter 'a' for a gamma distribution
        %sigma2 = std(Ptot_pdf2);
        gamma_b2 = 52.4768; %this is actually the shape parameter 'b' for a gamma distribution
    elseif strcmpi('ptot2T-',PTOT_SCENARIO2) == 1
        load(strcat('model_input',fs,'Ptot_pdf2.mat'), 'Ptot_pdf2')% This is the pdf fitted to all available station precip data (gamma dist). It will be modified as a trend and sampled below.
        %mu2 = mean(Ptot_pdf2);
        gamma_a2 = 1.65466;
        %sigma2 = std(Ptot_pdf2);
        gamma_b2 = 52.4768;
    end
end

load(strcat('model_input',fs,'Recess_pdf.mat'), 'Recess_pdf')% This is the pdf of storm gradient recession coefficiencts from Morin et al, 2005 (normal dist). It will be sampled below.
if strcmpi('Validation',MODE) == 1
    load(strcat('model_input',fs,'Easting.mat'), 'Easting')% This is the Longitudinal data for each gauge.
    load(strcat('model_input',fs,'Northing.mat'), 'Northing')% This is the Latitudinal data for each gauge. It will be sampled below.
    load(strcat('model_input',fs,'gauges.mat'), 'gauges')% This is the list of gauge numbers. It will be sampled below.
    load(strcat('model_input',fs,'gauge_elev.mat'), 'gauge_elev')% This is the list of gauge numbers. It will be sampled below.
    Zz = gauge_elev;
    numgauges = length(gauges);
elseif strcmpi('Simulation',MODE) == 1
    a = shaperead(strcat('boundary',fs,'boundary.shp')); %This is shapefile of the watershed boundary
    b = shaperead(strcat('point_elevations',fs,'point_elevations.shp')); %This is shapefile of the elevations at each 'gauging' point on the grid within the watershed boundary
    Yy = extractfield(a,'Y');
    Xx = extractfield(a,'X');
    Zz = extractfield(b,'RASTERVALU');
    [X1,Y1]=meshgrid(linspace(min(Xx),max(Xx),round((max(Xx)-min(Xx))/1000)),linspace(min(Yy),max(Yy),round((max(Yy)-min(Yy))/1000))); %creates a mesh with 1000m spacings
    isin=inpolygon(X1(:),Y1(:),Xx,Yy);
    Yin = Y1(isin);
    Xin = X1(isin); %creates a list of points, Xin, Yin on the grid that are within the watershed boundary and which will be used as 'gauging' output
    save Xin.txt Xin -ASCII
    save Yin.txt Yin -ASCII
    numgauges = length(Xin); %number of rain gauges in the basin. NOTE: In this version this produces output on a grid, rather than at real gauge locations.
end
load(strcat('model_input',fs,'X.mat'), 'X')% This is the Longitudinal data for each grid point. It will be sampled below to determine storm center location.
load(strcat('model_input',fs,'Y.mat'), 'Y')% This is the Latitudinal data for each grid point. It will be sampled below to determine storm center location.
load(strcat('model_input',fs,'fuzz.mat'), 'fuzz')% This is a vector of fuzzy tolerace values for intensity selection.
load(strcat('model_input',fs,'ET_monthly_day.mat'), 'ET_monthly_day')% This is a matrix of averaged daytime values of ET (mm) grouped as one column per month.
load(strcat('model_input',fs,'ET_monthly_night.mat'), 'ET_monthly_night')% This is a matrix of averaged nighttime values of ET (mm) grouped as one column per month.

%replaced above to compute these metrics based on the pdf used.
% mu = 207; %Mean of historical Ptotal pdf used as initial value.
% sigma = 64; %Std of historical Ptotal pdf used as initial value.
% mu2 = 207; %Mean of historical Ptotal pdf for season 2 used as initial value.
% sigma2 = 64; %Std of historical Ptotal pdf for season 2 used as initial value.

Gauges = 1:numgauges;
c = shaperead(strcat('bndry_buffer',fs,'bndry_buffer.shp')); %This is a shapefile of the watershed boundary buffer (5000 m)
Yyy = extractfield(c,'Y');
Xxx = extractfield(c,'X');
[Xx1,Yy1]=meshgrid(linspace(min(Xxx),max(Xxx),round((max(Xxx)-min(Xxx))/500)),linspace(min(Yyy),max(Yyy),round((max(Yyy)-min(Yyy))/500))); %creates a mesh with 500m spacings
isin=inpolygon(Xx1(:),Yy1(:),Xxx,Yyy);
Yyin = Yy1(isin);
Xxin = Xx1(isin); %creates a list of points, Xxin, Yyin on the grid that are within the watershed boundary buffer and which will be used as storm center locations

%These are elevation ranges for the 3 orographic groups in Walnut Gulch,
%AZ. They will need to be adapted for other basins (e.g., via hypsometric
%curve).
OroGrp1 = round(min(Zz)):1:1350;
OroGrp2 = 1350:1:1500;
OroGrp3 = 1501:1:round(max(Zz));

%lambda, kappa, and C are parameters of the intensity-duration curves from
%WGEW of the form:
%intensity1 = lambda*exp(delta*duration)+kappa*exp(gamma*duration)+C
delta = -0.508;
gamma = -0.008;
lambda = [642.2 578.0 513.8 449.5 385.3 321.1 256.9 192.7 128.4 64.1 21];
kappa = [93.1 83.8 74.5 65.2 55.9 46.6 37.2 27.9 18.6 9.3 0.9];
C = [4.5 4 3.5 3 2.5 2 1.5 1 0.5 0.25 0.05];%constant value of intensity for each recession curve at high duration values
%season2
%note: C=0 for season2
delta2 = -0.04;
gamma2 = -0.001038;
lambda2 = [72.71 72.71*.9 72.71*.8 72.71*.7 72.71*.6 72.71*.5 72.71*.4 72.71*.3 72.71*.2 72.71*.1 72.71*.05];
kappa2 = [5.432 5.432*.9 5.432*.8 5.432*.7 5.432*.6 5.432*.5 5.432*.4 5.432*.3 5.432*.2 5.432*.1 5.432*.05];
%Note: These curves were obtained by finding the maximum intensity value
%for each duration, followed by LOWESS smoothing with a span of 0.1. Next we fit a double
%exponential form to this LOWESS fit and obtained the curve
%parameters,lambda, kappa, and C (C=0 for season 2). Finally, we obtain the other curves as
%percentiles of the maximum fitted curve.

max_numstorms = 10000; %This is for initializing matrices. Trailing zeros are deleted from the matrix at the end of the code.

%Everything above this line is model input. Model operation begins below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create named and timestamped model output folder for the simulation set
t1 = [datestr(floor(now)) '_' datestr(now,13)];
t2=regexprep(t1,'[^\w'']',''); %for naming output directories and files by current date/time
t2 = t2(:,1:14);
tx = strcat(MODE(1:3),num2str(NUMSIMS),'_',num2str(NUMSIMYRS),'Y_',SEASONS, 'S_',PTOT_SCENARIO,'_',STORMINESS_SCENARIO,'_',PTOT_SCENARIO2,'_',STORMINESS_SCENARIO2,'_',ET_SCENARIO,'_');
tx0 = strcat('model_output',fs,tx,t2);
mkdir(tx0)

for simulations = 1:NUMSIMS %number of simulations to run. A separate timestamped output folder is created for each
    t0 = clock;
    t1a = [datestr(floor(now)) '_' datestr(now,13)];
    t2a=regexprep(t1a,'[^\w'']',''); %for naming output directories and files by current date/time
    t2a = t2a(:,1:14);
    tx2 = strcat(tx0,fs,t2a);
    mkdir(tx2)
    
    %next section of code simulates ET on a 12-hour basis (average day and
    %night values for each month of each simulation year) based on sampling
    %a pdf composed from regionally measured data
    
    jan_index = ~isnan(ET_monthly_day(:,1));
    feb_index = ~isnan(ET_monthly_day(:,2));
    mar_index = ~isnan(ET_monthly_day(:,3));
    apr_index = ~isnan(ET_monthly_day(:,4));
    may_index = ~isnan(ET_monthly_day(:,5));
    jun_index = ~isnan(ET_monthly_day(:,6));
    jul_index = ~isnan(ET_monthly_day(:,7));
    aug_index = ~isnan(ET_monthly_day(:,8));
    sep_index = ~isnan(ET_monthly_day(:,9));
    oct_index =  ~isnan(ET_monthly_day(:,10));
    nov_index = ~isnan(ET_monthly_day(:,11));
    dec_index = ~isnan(ET_monthly_day(:,12));
    ET_final = NaN;
    for yooo = 1:simyears
        yoo = num2str(yooo);
        ET_day_jan = datasample(ET_monthly_day(jan_index,1),31);
        ET_night_jan = datasample(ET_monthly_night(jan_index,1),31);
        ET_day_feb = datasample(ET_monthly_day(feb_index,2),28);
        ET_night_feb = datasample(ET_monthly_night(feb_index,2),28);
        ET_day_mar = datasample(ET_monthly_day(mar_index,3),31);
        ET_night_mar = datasample(ET_monthly_night(mar_index,3),31);
        ET_day_apr = datasample(ET_monthly_day(apr_index,4),30);
        ET_night_apr = datasample(ET_monthly_night(apr_index,4),30);
        ET_day_may = datasample(ET_monthly_day(may_index,5),31);
        ET_night_may = datasample(ET_monthly_night(may_index,5),31);
        ET_day_jun = datasample(ET_monthly_day(jun_index,6),30);
        ET_night_jun = datasample(ET_monthly_night(jun_index,6),30);
        ET_day_jul = datasample(ET_monthly_day(jul_index,7),31);
        ET_night_jul = datasample(ET_monthly_night(jul_index,7),31);
        ET_day_aug = datasample(ET_monthly_day(aug_index,8),31);
        ET_night_aug = datasample(ET_monthly_night(aug_index,8),31);
        ET_day_sep = datasample(ET_monthly_day(sep_index,9),30);
        ET_night_sep = datasample(ET_monthly_night(sep_index,9),30);
        ET_day_oct = datasample(ET_monthly_day(oct_index,10),31);
        ET_night_oct = datasample(ET_monthly_night(oct_index,10),31);
        ET_day_nov = datasample(ET_monthly_day(nov_index,11),30);
        ET_night_nov = datasample(ET_monthly_night(nov_index,11),30);
        ET_day_dec = datasample(ET_monthly_day(dec_index,12),31);
        ET_night_dec = datasample(ET_monthly_night(dec_index,12),31);
        jan_interleave = reshape([ET_day_jan(:) ET_night_jan(:)]',2*length(ET_day_jan), []); %#ok
        feb_interleave = reshape([ET_day_feb(:) ET_night_feb(:)]',2*length(ET_day_feb), []); %#ok
        mar_interleave = reshape([ET_day_mar(:) ET_night_mar(:)]',2*length(ET_day_mar), []); %#ok
        apr_interleave = reshape([ET_day_apr(:) ET_night_apr(:)]',2*length(ET_day_apr), []); %#ok
        may_interleave = reshape([ET_day_may(:) ET_night_may(:)]',2*length(ET_day_may), []); %#ok
        jun_interleave = reshape([ET_day_jun(:) ET_night_jun(:)]',2*length(ET_day_jun), []); %#ok
        jul_interleave = reshape([ET_day_jul(:) ET_night_jul(:)]',2*length(ET_day_jul), []); %#ok
        aug_interleave = reshape([ET_day_aug(:) ET_night_aug(:)]',2*length(ET_day_aug), []); %#ok
        sep_interleave = reshape([ET_day_sep(:) ET_night_sep(:)]',2*length(ET_day_sep), []); %#ok
        oct_interleave = reshape([ET_day_oct(:) ET_night_oct(:)]',2*length(ET_day_oct), []); %#ok
        nov_interleave = reshape([ET_day_nov(:) ET_night_nov(:)]',2*length(ET_day_nov), []); %#ok
        dec_interleave = reshape([ET_day_dec(:) ET_night_dec(:)]',2*length(ET_day_dec), []); %#ok
        eval(['ET_Y',yoo,' = [may_interleave; jun_interleave; jul_interleave; aug_interleave; sep_interleave; oct_interleave; nov_interleave; dec_interleave; jan_interleave; feb_interleave; mar_interleave; apr_interleave];'])
        eval(['ET_final = [ET_final; ET_Y',yoo,'];'])
        
    end
    ET_matrix = ET_final(2:length(ET_final));
    %This matrix contains a complete time series of day/night ET values for the length of the simulation.
    if  strcmpi('ET+',ET_SCENARIO) == 1
        ET_matrix = ET_matrix+ET_matrix.*ET_scaling_factor;
    elseif  strcmpi('ET-',ET_SCENARIO) == 1
        ET_matrix = ET_matrix-ET_matrix.*ET_scaling_factor;
        %Trends and seasonality of changes in ET not yet included
        %     elseif  strcmpi('ET_T+',ET_SCENARIO) == 1
        %         ET_matrix = ET_matrix+ET_matrix.*ET_scaling_factor;
        %     elseif  strcmpi('ET_T-',ET_SCENARIO) == 1
        %         ET_matrix = ET_matrix-ET_matrix.*ET_scaling_factor;
    end
    %figure
    
    %Initialize output matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Storm matrix = [Storm# StormArea StormDuration Int_DurCurve# StormIntensity #GaugesHit RecessionVal StormTotal UTM_Longitude UTM_Latitude Year SimTime]
    %The values in Storm_matrix have the following corresponding units: [# km^2 min # mm/hr mm/hr/km mm m m y hr]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Storm_matrix = zeros(max_numstorms*simyears,12); %based on presumed maximum number of storms per year (see above)
    Ptot_ann_global = zeros(1,simyears);
    Ptot_ann_global2 = zeros(1,simyears);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gauge matrix = [Year Storm# StormIntensity StormDuration StormTotal Ann_Cum_PTotal Int_Arr_Time SimTime]
    %The values in Gauge_matrix have the following corresponding units: [y # mm/hr min mm mm hr hr]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ii = 1:numgauges %initialize matrices of storm output for each gauge
        eval(['Gauge_matrix_',num2str(ii),' = zeros(max_numstorms*simyears,8);']) %based on presumed maximum number of storms per year = 300
    end
    
    Intensity_local_all = 0; %initialize all variables (concatenated matrices of generated output)
    Storm_totals_all = 0;
    Duration_local_all = 0;
    Gauges_hit_all = 0;
    storm_count = 0;
    master_storm_count = 0;
    storm_trend = 0;
    storm_trend2 = 0;
    calendar_time = 0; %tracks simulation time per year in hours
    for syear = 1:simyears
        calendar_year_time = 0; %tracks simulation time per year in hours
        clear leftovers*
        storm_trend = storm_trend+storminess_scaling_factor;
        Ptotal = 0; %#ok
        if  strcmpi('ptotC',PTOT_SCENARIO) == 1
            Ptot_ann_global(syear) = random(Ptot_pdf,1,true); %samples from the control normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
        elseif  strcmpi('ptot+',PTOT_SCENARIO) == 1
            Ptot_ann_global(syear) = random(Ptot_pdf_cc_up,1,true); %samples from the increased wetness normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
        elseif strcmpi('ptot-',PTOT_SCENARIO) == 1
            Ptot_ann_global(syear) = random(Ptot_pdf_cc_down,1,true); %samples from the decreased wetness normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
        elseif  strcmpi('ptotT+',PTOT_SCENARIO) == 1
            mu = mu+mu*ptot_scaling_factor; %scale the mean of the PTotal distribution each year by the scaling factor defined above and add it to the prior mean value.
            Ptot_pdf_cc_trendup = makedist('Normal',mu,sigma); %create new normal pdf with the updated mean scaled by the climate change trend
            Ptot_ann_global(syear) = random(Ptot_pdf_cc_trendup,1,true); %samples from normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
        elseif strcmpi('ptotT-',PTOT_SCENARIO) == 1
            mu = mu-mu*ptot_scaling_factor; %scale the mean of the PTotal distribution each year by the scaling factor defined above and subtract it from the prior mean value.
            Ptot_pdf_cc_trenddown = makedist('Normal',mu,sigma); %create new normal pdf with the updated mean scaled by the climate change trend
            Ptot_ann_global(syear) = random(Ptot_pdf_cc_trenddown,1,true); %samples from normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
        end
        Storm_total_local_year = zeros(1,numgauges);
        ann_cum_Ptot_gauge = zeros(1,numgauges);
        for storm = 1:max_numstorms
            int_arr_val = random(Int_arr_pdf,1,true); %Samples from distribution of interarrival times (hr). This can be used to develop STORM output for use in rainfall-runoff models or any water balance application.
            clear cx cy r mask_name aa bb cc dd* North_hit North East East_hit rain_int_gauge intensity_val duration_val recess_val Ptotal* gdist center_val* int_dur_curve_num*
            rain_int_gauge = zeros(1,numgauges);
            center_val_X = datasample(Xxin,1); %sample uniformly from storm center matrix from grid with 500 m spacings within a 5000 m buffer around the basin.
            center_val_Y = datasample(Yyin,1);
            North = center_val_Y;
            East = center_val_X;
            area_val = random(Area_pdf,1,true); %Samples from distribution of storm areas in m2
            cx = East; %value of coord should be set to storm center selected (same below)
            cy = North;
            r = sqrt(area_val/pi); %value here should be selected based on area above in meters to match the UTM values in North and East vectors.
            if strcmpi('Validation',MODE) == 1
                mask_name = (((Easting(:)-cx).^2 + (Northing(:)-cy).^2) <= r^2); % %determine which gauges are hit by Euclidean geometry
            elseif strcmpi('Simulation',MODE) == 1
                mask_name = (((Xin(:)-cx).^2 + (Yin(:)-cy).^2) <= r^2); %determine which grid locations are hit by Euclidean geometry
            else
            end
            aa = find(mask_name == 1);
            if isempty(aa) %this short circuits the storm loop in the case that the storm does not affect any 'gauging' location
                continue
            end
            storm_count = storm_count+1;
            master_storm_count = master_storm_count+1;
            gauges_hit = Gauges(aa);
            num_gauges_hit = length(gauges_hit);
            Storm_matrix(master_storm_count,1) = master_storm_count; %storm #
            Storm_matrix(master_storm_count,2) = area_val;
            Storm_matrix(master_storm_count,9) = cx;
            Storm_matrix(master_storm_count,10) = cy;
            Storm_matrix(master_storm_count,11) = syear;
            %this routine below allows for orography in precip by first
            %determining the closest gauge and then determining its
            %orographic grouping
            if strcmpi('Validation',MODE) == 1
                gdist = (Easting(:)-cx).^2 + (Northing(:)-cy).^2;
            elseif strcmpi('Simulation',MODE) == 1
                gdist = (Xin(:)-cx).^2 + (Yin(:)-cy).^2;
            end
            bb = min(gdist);
            closest_gauge = round(Zz(gdist == bb)); %this will be compared against orographic gauge groupings to determine the appropriate set of intensity-duration curves
            Storm_matrix(master_storm_count,6) = num_gauges_hit;
            Gauges_hit_all = [Gauges_hit_all, gauges_hit]; %#ok<AGROW>
            
            %this routine below determines to which orographic group the
            %closest gauge to the storm center belongs to, and censors the
            %number of curves accordingly missing top curve in GR1, top and
            %bottom curves for GR2, and bottom curve for GR3
            dd = ismember(OroGrp1,closest_gauge); %new version of orography compares local 'gauge' elevation to elevation bands called OroGrp, defined above
            if ~isempty(dd)
                baa = 'a';
            end
            dd = ismember(OroGrp2,closest_gauge);
            if ~isempty(dd)
                %int_dur_curve_num = (2:numcurves-1)'; % these were empirically determined based on data from WG (monsoon rainfall only)-middle orographic group.
                baa = 'b';
            end
            dd = ismember(OroGrp3,closest_gauge);
            if ~isempty(dd)
                baa = 'c';
            end
            int_dur_curve_num = (1:numcurves)';
            duration_val = random(Duration_pdf,1,true);
            duration_val = round(duration_val); %round to nearest minute for consistency with measured data
            Storm_matrix(master_storm_count,3) = duration_val;
            calendar_year_time = calendar_year_time+int_arr_val+duration_val/60; %hours
            calendar_time = calendar_time+int_arr_val+duration_val/60; %hours
            Storm_matrix(master_storm_count,12) = calendar_time; %stored cumulative time of simulation
            
            
            % original curve# probs for 30%-20%-10%: [0.0636 0.0727 0.0819 0.0909 0.0909 0.0909 0.0909 0.0909 0.1001 0.1090 0.1182]
            % original curve# probs are modified as below
            if baa == 'a'
                int_dur_curve_val = datasample(int_dur_curve_num,1,'weights',[0.0318 0.0759 0.0851 0.0941 0.0941 0.0941 0.0941 0.0941 0.1033 0.1121 0.1213]); %weights to reflect probabilities that favor selection of lower curves.
            elseif baa == 'b'
                int_dur_curve_val = datasample(int_dur_curve_num,1,'weights',[0.0478 0.0778 0.0869 0.0959 0.0959 0.0959 0.0959 0.0959 0.1051 0.1141 0.0888]); %each of these is adapted based on its orographic grouping
            elseif baa == 'c'
                int_dur_curve_val = datasample(int_dur_curve_num,1,'weights',[0.0696 0.0786 0.0878 0.0968 0.0968 0.0968 0.0968 0.0968 0.1060 0.1149 0.0591]);
            else
            end
            Storm_matrix(master_storm_count,4) = int_dur_curve_val;
            intensity_val = lambda(int_dur_curve_val)*exp(delta*duration_val)+kappa(int_dur_curve_val)*exp(gamma*duration_val)+ C(int_dur_curve_val);  %these curves are based on empirical data from WG.
            fuzz_int_val = datasample(fuzz,1);
            intensity_val2 = intensity_val+fuzz_int_val; %allowing for +/-5 mm/hr fuzzy tolerance around selected intensity
            if intensity_val2 < 0.5 %cannot have zero or negative intensity
                intensity_val = 0;
            elseif intensity_val2 < 1 && intensity_val >=0.5
                intensity_val = 1;
            else
                intensity_val = intensity_val2;
            end
            if strcmpi('stormsT+',STORMINESS_SCENARIO) == 1
                intensity_val = round(intensity_val + intensity_val*storm_trend); %storminess trend is applied and its effect rises each year of simulation
            elseif strcmpi('stormsT-',STORMINESS_SCENARIO) == 1
                intensity_val = round(intensity_val - intensity_val*storm_trend);
            elseif strcmpi('storms+',STORMINESS_SCENARIO) == 1
                intensity_val = round(intensity_val + intensity_val*storm_stepchange); %storminess change is applied as a step change uniformly over the simulation
            elseif strcmpi('storms-',STORMINESS_SCENARIO) == 1
                intensity_val = round(intensity_val - intensity_val*storm_stepchange);
            end
            Storm_matrix(master_storm_count,5) = intensity_val;
            
            %area to determine which gauges are hit
            recess_val = random(Recess_pdf,1,true); %this pdf of recession coefficients determines how intensity declines with distance from storm center (see below)
            Storm_matrix(master_storm_count,7) = recess_val;
            for j = 1:numgauges %determine cartesian distances to all hit gauges and associated intensity values at each gauge hit by the storm
                if ismember(j,gauges_hit)
                    if strcmpi('Validation',MODE) == 1
                        gauge_dist = sqrt(((Easting(j)-cx).^2 + (Northing(j)-cy).^2));
                    elseif strcmpi('Simulation',MODE) == 1
                        gauge_dist = sqrt(((Xin(j)-cx).^2 + (Yin(j)-cy).^2));
                    end
                    gauge_dist_km = gauge_dist./1000; %change to km
                    rain_int_gauge(j) = intensity_val*exp(-2*(recess_val^2)*gauge_dist_km.^2); %Rodriguez-Iturbe et al., 1986; Morin et al., 2005 but sampled from a distribution
                elseif ~ismember(j,gauges_hit)
                    rain_int_gauge(j) = 0; %give zero intensity to gauges not hit
                end
                ann_cum_Ptot_gauge(j) = ann_cum_Ptot_gauge(j)+rain_int_gauge(j)*duration_val/60;
            end
            
            for jj = 1:numgauges
                eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,1) = syear;']) %year
                eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,2) = master_storm_count;']) %storm #
                eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,3) = rain_int_gauge(jj);'])
                eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,4) = duration_val;'])
                eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,5) = rain_int_gauge(jj)*duration_val/60;']) %storm total
                eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,6) = ann_cum_Ptot_gauge(jj);']) %ann cum total (Ptot)
                eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,7) = int_arr_val;']) %interarrival time in hours
                eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,8) = calendar_time;']) %simulation time in hours
            end
            Intensity_local_all = [Intensity_local_all,rain_int_gauge]; %#ok<AGROW> %collect into vector of all simulated intensities (at all gauges)
            dur_step(1:numgauges) = duration_val;
            Duration_local_all = [Duration_local_all,dur_step]; %#ok<AGROW>
            Storm_total_local_year(storm,1:numgauges) = rain_int_gauge.*duration_val/60; % %collect storm total data for all gauges into rows by storm
            Storm_totals_all = [Storm_totals_all,rain_int_gauge.*duration_val/60]; %#ok<AGROW>
            Storm_matrix(master_storm_count,8) = intensity_val*duration_val/60;
            Ptotal = zeros(1,numgauges);
            for k = 1:numgauges
                Ptotal(k) = nansum(Storm_total_local_year(:,k)); % %sum the annual storm total at each gauge
            end
            
            Ptotal_test = find(nanmedian(Ptotal) > Ptot_ann_global(syear), 1); %once the median of all gauges exceeds the selected annual storm total, a new simulation year beginsPtotal_test = find(Ptotal > Ptot_ann_global(syear));    %once the first gauge exceeds the selected annual storm total, a new simulation year begins
            if ~isempty(Ptotal_test)
                for l = 1:numgauges
                    Ptotal_local(syear,l) = Ptotal(l); %#ok
                end
                break %end storm loop and start a new simulation year
            end
            eval(['Storm_total_local_year_',num2str(syear),'(storm,1:numgauges) = Storm_total_local_year(storm,1:numgauges);']) %collect all local annual storm totals for each gauge.
            
        end %storm loop
        stormcount = storm+1;
        %         leftstuff = find(Gauge_matrix_1(:,1) == syear);
        %         leftstuff2 = length(leftstuff);
        %         leftovers1 = 153*24*60-sum(Gauge_matrix_1(leftstuff,7).*60)-sum(Gauge_matrix_1(leftstuff,4)); %#ok %remaining monsoon time in minutes (not occupied by storms or interstorm periods)
        %         leftovers11 = leftovers1/leftstuff2/60;
        %redistribute leftover time in the season equally as additions to interstorm periods
        %Gauge_matrix_1(leftstuff,7) = Gauge_matrix_1(leftstuff,7) + leftovers11;
        %         for j = 1:numgauges
        %             eval(['Gauge_matrix_',j,'(leftstuff,7) = Gauge_matrix_',j,'(leftstuff,7) + leftovers11'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SEASON 2
        %calendar_time = 0; %tracks simulation time per year in hours
        storm_trend2 = storm_trend2+storminess_scaling_factor2;
        if  strcmpi('two',SEASONS) == 1
            Ptotal = 0; %#ok
            if  strcmpi('ptot2C',PTOT_SCENARIO2) == 1
                Ptot_ann_global2(syear) = random(Ptot_pdf2,1,true); %samples from the control normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
            elseif  strcmpi('ptot2+',PTOT_SCENARIO2) == 1
                Ptot_ann_global2(syear) = random(Ptot_pdf_cc_up2,1,true); %samples from the increased wetness normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
            elseif strcmpi('ptot2-',PTOT_SCENARIO2) == 1
                Ptot_ann_global2(syear) = random(Ptot_pdf_cc_down2,1,true); %samples from the decreased wetness normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
            elseif  strcmpi('ptot2T+',PTOT_SCENARIO2) == 1
                gamma_a2 = gamma_a2+ptot_scaling_factor2;
                Ptot_pdf_cc_trendup2 = makedist('Gamma','a',gamma_a2,'b',gamma_b2); %create new normal pdf with the updated mean scaled by the climate change trend
                Ptot_ann_global2(syear) = random(Ptot_pdf_cc_trendup2,1,true); %samples from normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
            elseif strcmpi('ptot2T-',PTOT_SCENARIO2) == 1
                gamma_a2 = gamma_a2-ptot_scaling_factor2;
                Ptot_pdf_cc_trenddown2 = makedist('Gamma','a',gamma_a2,'b',gamma_b2); %create new normal pdf with the updated mean scaled by the climate change trend
                Ptot_ann_global2(syear) = random(Ptot_pdf_cc_trenddown2,1,true); %samples from normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
            end
            
            for storm = stormcount:max_numstorms  %continue counting storms in season 2
                int_arr_val = random(Int_arr_pdf2,1,true); %Samples from distribution of interarrival times (hr). This can be used to develop STORM output for use in rainfall-runoff models or any water balance application.
                clear cx cy r mask_name aa bb cc dd* North_hit North East East_hit rain_int_gauge intensity_val duration_val recess_val Ptotal* gdist center_val* int_dur_curve_num*
                rain_int_gauge = zeros(1,numgauges);
                center_val_X = datasample(Xxin,1); %sample uniformly from storm center matrix from grid with 500 m spacings within a 5000 m buffer around the basin.
                center_val_Y = datasample(Yyin,1);
                North = center_val_Y;
                East = center_val_X;
                area_val = random(Area_pdf2,1,true); %Samples from distribution of storm areas in m2
                cx = East; %value of coord should be set to storm center selected (same below)
                cy = North;
                r = sqrt(area_val/pi); %value here should be selected based on area above in meters to match the UTM values in North and East vectors.
                if strcmpi('Validation',MODE) == 1
                    mask_name = (((Easting(:)-cx).^2 + (Northing(:)-cy).^2) <= r^2); % %determine which gauges are hit by Euclidean geometry
                elseif strcmpi('Simulation',MODE) == 1
                    mask_name = (((Xin(:)-cx).^2 + (Yin(:)-cy).^2) <= r^2); %determine which grid locations are hit by Euclidean geometry
                else
                end
                aa = find(mask_name == 1);
                if isempty(aa) %this short circuits the storm loop in the case that the storm does not affect any 'gauging' location
                    continue
                end
                storm_count = storm_count+1;
                master_storm_count = master_storm_count+1;
                gauges_hit = Gauges(aa);
                num_gauges_hit = length(gauges_hit);
                Storm_matrix(master_storm_count,1) = master_storm_count; %storm #
                Storm_matrix(master_storm_count,2) = area_val;
                Storm_matrix(master_storm_count,9) = cx;
                Storm_matrix(master_storm_count,10) = cy;
                Storm_matrix(master_storm_count,11) = syear;
                %this routine below allows for orography in precip by first determining the closest gauge and then determining its orographic grouping
                if strcmpi('Validation',MODE) == 1
                    gdist = (Easting(:)-cx).^2 + (Northing(:)-cy).^2;
                elseif strcmpi('Simulation',MODE) == 1
                    gdist = (Xin(:)-cx).^2 + (Yin(:)-cy).^2;
                end
                bb = min(gdist);
                closest_gauge = round(Zz(gdist == bb)); %this will be compared against orographic gauge groupings to determine the appropriate set of intensity-duration curves
                Storm_matrix(master_storm_count,6) = num_gauges_hit;
                Gauges_hit_all = [Gauges_hit_all, gauges_hit]; %#ok<AGROW>
                
                %this routine below determines to which orographic group the closest gauge to the storm center belongs to, and censors the
                %number of curves accordingly missing top curve in GR1, top and bottom curves for GR2, and bottom curve for GR3
                dd = ismember(OroGrp1,closest_gauge); %new version of orography compares local 'gauge' elevation to elevation bands called OroGrp, defined above
                if ~isempty(dd)
                    %if ~isempty(dd1,1)
                    baa = 'a';
                end
                dd = ismember(OroGrp2,closest_gauge);
                if ~isempty(dd)
                    baa = 'b';
                end
                dd = ismember(OroGrp3,closest_gauge);
                if ~isempty(dd)
                    baa = 'c';
                end
                int_dur_curve_num = (1:numcurves)';
                duration_val = random(Duration_pdf2,1,true);
                duration_val = round(duration_val); %round to nearest minute for consistency with measured data
                Storm_matrix(master_storm_count,3) = duration_val;
                calendar_time = calendar_time+int_arr_val+duration_val/60;
                Storm_matrix(master_storm_count,12) = calendar_time; %stored cumulative time of simulation
                
                if baa == 'a'
                    int_dur_curve_val = datasample(int_dur_curve_num,1,'weights',[0.0318 0.0759 0.0851 0.0941 0.0941 0.0941 0.0941 0.0941 0.1033 0.1121 0.1213]); %weights to reflect probabilities that favor selection of lower curves.
                elseif baa == 'b'
                    int_dur_curve_val = datasample(int_dur_curve_num,1,'weights',[0.0478 0.0778 0.0869 0.0959 0.0959 0.0959 0.0959 0.0959 0.1051 0.1141 0.0888]); %each of these is adapted based on its orographic grouping
                elseif baa == 'c'
                    int_dur_curve_val = datasample(int_dur_curve_num,1,'weights',[0.0696 0.0786 0.0878 0.0968 0.0968 0.0968 0.0968 0.0968 0.1060 0.1149 0.0591]);
                else
                end
                Storm_matrix(master_storm_count,4) = int_dur_curve_val;
                intensity_val = lambda2(int_dur_curve_val)*exp(delta2*duration_val)+kappa2(int_dur_curve_val)*exp(gamma2*duration_val) + C(int_dur_curve_val);  %these curves are based on empirical data from WG.
                fuzz_int_val = datasample(fuzz,1);
                intensity_val2 = intensity_val+fuzz_int_val; %allowing for +/-5 mm/hr fuzzy tolerance around selected intensity
                if intensity_val2 < 0.5 %cannot have zero or negative intensity
                    intensity_val = 0;
                elseif intensity_val2 < 1 && intensity_val >=0.5
                    intensity_val = 1;
                else
                    intensity_val = intensity_val2;
                end
                if strcmpi('storms2T+',STORMINESS_SCENARIO2) == 1
                    intensity_val = round(intensity_val + intensity_val*storm_trend2); %storminess trend is applied and its effect rises each year of simulation
                elseif strcmpi('storms2T-',STORMINESS_SCENARIO2) == 1
                    intensity_val = round(intensity_val - intensity_val*storm_trend2);
                elseif strcmpi('storms2+',STORMINESS_SCENARIO2) == 1
                    intensity_val = round(intensity_val + intensity_val*storm_stepchange2); %storminess change is applied as a step change uniformly over the simulation
                elseif strcmpi('storms2-',STORMINESS_SCENARIO2) == 1
                    intensity_val = round(intensity_val - intensity_val*storm_stepchange2);
                end
                Storm_matrix(master_storm_count,5) = intensity_val;
                
                %area to determine which gauges are hit
                recess_val = random(Recess_pdf,1,true); %this pdf of recession coefficients determines how intensity declines with distance from storm center (see below)
                Storm_matrix(master_storm_count,7) = recess_val;
                for j = 1:numgauges %determine cartesian distances to all hit gauges and associated intensity values at each gauge hit by the storm
                    if ismember(j,gauges_hit)
                        if strcmpi('Validation',MODE) == 1
                            gauge_dist = sqrt(((Easting(j)-cx).^2 + (Northing(j)-cy).^2));
                        elseif strcmpi('Simulation',MODE) == 1
                            gauge_dist = sqrt(((Xin(j)-cx).^2 + (Yin(j)-cy).^2));
                        end
                        gauge_dist_km = gauge_dist./1000; %change to km
                        rain_int_gauge(j) = intensity_val*exp(-2*(recess_val^2)*gauge_dist_km.^2); %Rodriguez-Iturbe et al., 1986; Morin et al., 2005 but sampled from a distribution
                    elseif ~ismember(j,gauges_hit)
                        rain_int_gauge(j) = 0; %give zero intensity to gauges not hit
                    end
                    ann_cum_Ptot_gauge(j) = ann_cum_Ptot_gauge(j)+rain_int_gauge(j)*duration_val/60;
                end
                
                for jj = 1:numgauges
                    eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,1) = syear;']) %year
                    eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,2) = master_storm_count;']) %storm #
                    eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,3) = rain_int_gauge(jj);'])
                    eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,4) = duration_val;'])
                    eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,5) = rain_int_gauge(jj)*duration_val/60;']) %storm total
                    eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,6) = ann_cum_Ptot_gauge(jj);']) %ann cum total (Ptot)
                    eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,7) = int_arr_val;']) %interarrival time in hours
                    eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,8) = calendar_time;']) %total simulation time in hours
                end
                Intensity_local_all = [Intensity_local_all,rain_int_gauge]; %#ok<AGROW> %collect into vector of all simulated intensities (at all gauges)
                dur_step(1:numgauges) = duration_val;
                Duration_local_all = [Duration_local_all,dur_step]; %#ok<AGROW>
                Storm_total_local_year(storm,1:numgauges) = rain_int_gauge.*duration_val/60; % %collect storm total data for all gauges into rows by storm
                Storm_totals_all = [Storm_totals_all,rain_int_gauge.*duration_val/60]; %#ok<AGROW>
                Storm_matrix(master_storm_count,8) = intensity_val*duration_val/60;
                Ptotal = zeros(1,numgauges);
                for k = 1:numgauges
                    Ptotal(k) = nansum(Storm_total_local_year(:,k)); % %sum the annual storm total at each gauge
                end
                Ptotal_test = find(nanmedian(Ptotal) > Ptot_ann_global2(syear), 1); %once the median of all gauges exceeds the selected annual storm total, a new simulation year beginsPtotal_test = find(Ptotal > Ptot_ann_global(syear));    %once the first gauge exceeds the selected annual storm total, a new simulation year begins
                if ~isempty(Ptotal_test)
                    for l = 1:numgauges
                        Ptotal_local(syear,l) = Ptotal(l);
                    end
                    break %end storm lopp and start a new simulation year
                end
                eval(['Storm_total_local_year_',num2str(syear),'(storm,1:numgauges) = Storm_total_local_year(storm,1:numgauges);']) %collect all local annual storm totals for each gauge.
            end %storm loop
        end %end of 2nd season loop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end %year loop
    %hold on plot(Easting,Northing,'o')
    %grid on
    Storm_matrix = Storm_matrix(Storm_matrix(:,9)>0,:); % %gets rid of trailing zeros from initialized matrix
    AB = find(Gauge_matrix_1(:,2)>0); %#ok
    for CC = 1:numgauges
        eval(['Gauge_matrix_',num2str(CC),' = Gauge_matrix_',num2str(CC),'(AB,:);'])
    end
    leftstuff2 = max(Gauge_matrix_1(:,8)); %simulated time in hours
    nn = length(Gauge_matrix_1(:,8));
    leftovers = 365.25*NUMSIMYRS*24-leftstuff2; % %remaining annual time in hours (not occupied by storms or interstorm periods)
    leftovers22 = leftovers/nn; %adjustment to add to each interarrival time in hours
    Gauge_matrix_1back = Gauge_matrix_1;
    Gauge_matrix_1(:,7) = Gauge_matrix_1(:,7) + leftovers22;
    cal_time2 = 0;
    for jjj = 1:nn
        Gauge_matrix_1(jjj,8) = cal_time2 + Gauge_matrix_1(jjj,7) + Gauge_matrix_1(jjj,4)./60;
        cal_time2 = Gauge_matrix_1(jjj,8);
    end
    for jj = 2:numgauges
        eval(['Gauge_matrix_',num2str(jj),'(:,7) = Gauge_matrix_1(:,7);'])
        eval(['Gauge_matrix_',num2str(jj),'(:,8) = Gauge_matrix_1(:,8);'])
    end
    Storm_matrix(:,12) = Gauge_matrix_1(:,8);
    Gauges_hit_all(Gauges_hit_all == 0) = NaN; %#ok
    GZ = find(Intensity_local_all>0);
    Intensity_all = Intensity_local_all(GZ); %#ok
    Duration_all = Duration_local_all(GZ); %#ok
    Storm_totals_all = Storm_totals_all(GZ); %#ok
    Ptot_ann_global = Ptot_ann_global(2:length(Ptot_ann_global)); %#ok
    if  strcmpi('two',SEASONS) == 1
        Ptot_ann_global2 = Ptot_ann_global2(2:length(Ptot_ann_global2)); %#ok
    end
    Gauges_hit_all = Gauges_hit_all(2:length(Gauges_hit_all)); %#ok
    
    eval(['save ',tx2,fs,'Ptot_ann_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_global_',t2,' Ptot_ann_global'])
    if  strcmpi('two',SEASONS) == 1
        eval(['save ',tx2,fs,'Ptot_ann2_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_global2_',t2,' Ptot_ann_global2']);
    end
    eval(['save ',tx2,fs,'Storm_matrix_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_',t2,' Storm_matrix'])
    eval(['save ',tx2,fs,'Gauges_hit_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_all_',t2,' Gauges_hit_all'])
    eval(['save ',tx2,fs,'Storm_totals_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_all_',t2,' Storm_totals_all'])
    eval(['save ',tx2,fs,'Intensity_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_selected_',t2,' Intensity_all'])
    eval(['save ',tx2,fs,'Duration_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_selected_',t2,' Duration_all'])
    
    eval(['save ',tx2,fs,'ET_matrix_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_',t2,' ET_matrix'])
    
    for kk = 1:numgauges
        eval(['save ',tx2,fs,'Gauge_matrix',num2str(kk),'_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_',t2,' Gauge_matrix_',num2str(kk);])
    end
    
    boo = etime(clock,t0);
    runtime_seconds = boo; %#ok
    runtime_minutes = boo/60 %#ok
    
end



