clc;
clear all;
close all;

%% Inputs

% Export global variables from script
global numPass avgW payload dist cruiseSpeed Vfwd cruiseTime Vmaxfwd altitude h hoverTime climbDist rateClimb climbTime reserve
global atmData alt temp pressure density vSound kinVisc rho nu BLdata advRatio bladeLoading flatPlateData flatPlateWeightData flatPlateAreaData
global Ptotal_fwd Ptotal_hover

% Import special global variables
global Wg_new We_new W_battery Ec_tot Ptotal_hover Ptotal_fwd R Omega energies weights radii hoverpowers

% Mission Profile
numPass = 2; % number of passengers (including pilot)
avgW = 200; % average weight of person [lbs]
payload = avgW * numPass; % total payload weight [lbs]

dist = 25*1609; % distance [m] (25 miles)
cruiseSpeed = 100 * .5144; % [m/s] 100 knots
Vfwd = cruiseSpeed; % forward velocity [m/s]
cruiseTime = dist/Vfwd; % cruise time[s]

Vmaxfwd = 120 * .5144; % maximum forward velocity [m/s]
altitude = 1000; % altitude [ft]

h = 10 * 0.3048; % hover height [m]
hoverTime = 240 ; % [s]
 
climbDist = altitude * 0.3048 - h; % vertical climb/landing distance [m]
rateClimb = 2.54; % rate of climb [m/s] (500 fpm)
climbTime = climbDist/rateClimb; % climbing time [s]

reserve = 20 * 60; % reserve requirement [s] (20 min) (FAA requirements)

%% Data Loading
% Atmospheric Data for Interpolation based on Altitude
atmData = xlsread('atmospheredata.xlsx'); % load atmospheric data
alt = atmData(:,1); % altitude [ft]
temp = atmData(:,2); % temperature [R]
pressure = atmData(:,3); % pressure [lb/ft^2]
density = atmData(:,4); % density [slugs/ft^3]
vSound = atmData(:,5); % speed of sound [ft/s]
kinVisc = atmData(:,6); % kinematic viscosity [ft^2/s]

rho = interp1(alt, density, altitude) * 515.379; % density of air [kg/m^3]
nu = interp1(alt, kinVisc, altitude) * 0.092903; % kinematic viscosity of air [m^2/s]

BLdata = xlsread('bladeloadingdata.xlsx'); % load data for blade loading vs. advance ratio
advRatio = BLdata(:,1); % advance ratio data
bladeLoading = BLdata(:,2); % blade loading data

flatPlateData = xlsread('flatplateareadata.xlsx'); % load data for flat plate area vs. gross weight
flatPlateWeightData = flatPlateData(:,1); % gross weight [lbs]
flatPlateAreaData = flatPlateData(:,2); % equivalent flat plate area [ft^2]

%% Analysis

energies = [];
weights = [];
radii = [];
hoverpowers = [];

% Run analysis code and output plots
% First parameter - sweep test case
    % 'passenger'
    % 'speed'
    % 'distance'
    % 'hover'
% Second parameter - helicopter type
    % 'compound'
    % 'electric'

analysis('speed', 'electric');

%% Outputs
GrossWeightLbs = Wg_new*0.2247
EmptyWeightLbs = We_new*0.2247
BatteryWeightLbs = W_battery*0.2247
TotalEnergyCapacity_kWh = Ec_tot/1000
TotalHoverPower_kW = Ptotal_hover/1000
TotalCruisePower_kW = Ptotal_fwd/1000
Radius_ft = R*3.28
RPM = Omega*9.549