clc;
clear all;
%close all;

%% Inputs

% % Mission Profile
% 
numPass = 2; % number of passengers (including pilot)
avgW = 200; % average weight of person [lbs]
payload = avgW * numPass; % total payload weight [lbs]
% 
% payload = 200 ; % [lbs]
dist = 25*1609; % distance [m] (25 miles)
cruiseSpeed = 100 * .5144; % [m/s] 100 knots
Vfwd = cruiseSpeed; % forward velocity [m/s]
cruiseTime = dist/Vfwd; % cruise time[s]
% cruiseTime = 30*60;
% dist = cruiseSpeed * cruiseTime;
% 
Vmaxfwd = 120 * .5144; % maximum forward velocity [m/s]
altitude = 1000; % altitude [ft]
% 
h = 10 * 0.3048; % hover height [m]
hoverTime = 240 ; % [s]
 
climbDist = altitude * 0.3048 - h; % vertical climb/landing distance [m]
rateClimb = 2.54; % rate of climb [m/s] (500 fpm)
climbTime = climbDist/rateClimb; % climbing time [s]
% 
reserve = 20 * 60; % reserve requirement [s] (20 min) (FAA requirements)
% 
Ed = 140.6; % Energy Density [W*h/kg] (Kokam LiPo battery energy density - used in Freidrich and Robertson 2015)

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

Ed_sweep = [144 250 400];
% pass_sweep = [2 4 6 8];

for j = 1:length(Ed_sweep)
    Ed = Ed_sweep(j);

    %numPass = pass_sweep(j);
    %payload = avgW * numPass;
    
    if Ed == 144
%         passengers = 1:10;
%         speeds = [60:120]*.5144;
         distances = [5:69]*1609;
%        hovers = 10:10:1220;
    elseif Ed == 250
        %passengers = 1:10;
                %speeds = [35:120]*.5144;
        distances = [5:183]*1609; 
%        hovers = 10:10:3810;
    elseif Ed ==400
        %passengers = 2:10;
                %speeds = [25:120]*.5144;
         distances = [5:329]*1609;
%       hovers = 10:10:6740;
    end
    
    for i = 1:length(distances)
        
        %numPass = passengers(i);
        %payload = avgW * numPass;
        
%         cruiseSpeed = speeds(i);
%         Vfwd = cruiseSpeed;
%         cruiseTime = dist/Vfwd;
%         
          dist = distances(i);
          cruiseTime = dist/Vfwd;
%          
%        hoverTime = hovers(i);

        
        % HELICOPTER
        
        % Following process from: Guide for Conceptual Helicopter Design (Kee)
        
        % Main Rotor Design
        
        % Determine initial weight estimate
        % Relation created from electric helicopter statistics gathered
        m = 1.175; % slope of the line
        b = 104.93; % y-intercept of line
        Wg_lbs = (payload - b/m)/(1 - 1/m); % Initial calculation of gross weight [lbs]
        Wg_init = Wg_lbs * 4.45; % convert weight to SI units [N]
        We_init = (Wg_lbs - payload) * 4.45; % empty weight [N]
        
        % Calculate the maximum tip velocity
        Mmaxtip = 0.65; % main rotor tip Mach number (assumption - value given in Kee)
        a = interp1(alt, vSound, altitude) * 0.3048; % speed of sound [m/s]
        Vmaxtip = Mmaxtip * a; % maximum tip velocity [m/s]
        
        % Determine number of main rotor blades
        numBlades = 2; % number of blades (assumption - typical for lighter weight helicopters - Leishman)
%         numBlades = 4; % run with 4 blades instead of 2
        
        % Estimate Reynolds number
        R_init = (0.0011 * (Wg_init * 0.2247) + 11.496) * 0.3048; % rotor radius [m] (relation determined from data in Leishman)
        T_init = Wg_init; % set thrust equal to gross weight [N]
        A_init = pi * R_init^2; % disk area [m^2]
        Omega_init = Vmaxtip / R_init; % maximum rotational velocity [rad/s]
        Vtip_init = Omega_init * R_init; % tip velocity [m/s]
        Ct_init = T_init/(A_init * rho * Vtip_init^2); % thrust coefficient
        muMax_init = Vmaxfwd / Vtip_init; % maximum advance ratio
        BL_init = interp1(advRatio, bladeLoading, muMax_init); % maximum blade loading
        sigma_init = Ct_init / BL_init; % blade solidity
        c_init = (sigma_init * pi * R_init)/numBlades; % chord [m]
        V_blade = sqrt(Vmaxfwd^2 + (Omega_init * 0.75 * R_init)^2); % velocity of blade at 75% of radius [m/s]
        Re = V_blade * c_init/nu; % Reynolds number
        
        % Choose an airfoil section
        % Assumption - NACA 0012 historically used most often (Kee)
        
        % Determine Average Lift Curve Slope and Average Profile Drag Coefficient
        Cl_alpha = 2 * pi; % lift curve slope (thin airfoil theory)
        Cdo = 0.011; % 0.00512; % profile drag coefficient (at 0 angle of attack because symmetric airfoil)
        
       
        
        %% Analysis
        
        % energies = [5:1:110] * 1000;
        % grossweights = [];
        % batteries = [];
        % radii = [];
        % energycap = [];
        %
        % for i = 1:length(energies)
        %     Ec = energies(i);
        %
        % START OUTER LOOP
        cond1 = 1; % conditional for when to exit loop
        Wg = Wg_init; % initialize gross weight [N]
        We = We_init; % initialize empty weight [N]
        P_hover = 0; % initialize power [W]
        Mbatt = 100; %kg 
        Ec = Ed * Mbatt; % initial guess for energy capacity [W*hr]
        FM = 0; % used for changing variables for figure of merit adjustments
        Rflag = 0; % used for changing radius for aspect ratio adjustments
        powerCalc = 0; % used for when to move to final power calculations
        Omega = Omega_init;
        skipInitial = 0;
        
         
        mainLoop_counter = 0;
        innerLoop_counter = 0;
        Ec_counter = 0;
        
        % Ecs = [Ec];
        
        while cond1 == 1
  
            mainLoop_counter = mainLoop_counter+1;
            
            if skipInitial == 0
                
                % Determine rotor radius
                R = (0.0011 * (Wg * 0.2247) + 11.496) * 0.3048; % rotor radius [m] (relation determined from data in Leishman)
                
                % Determine maximum rotational velocity
                Omega_max = Vmaxtip / R; % maximum rotational velocity [rad/s]
                
                % Determine thrust coefficient
                T = Wg; % set thrust equal to gross weight [N]
                A = pi * R^2; % disk area [m^2]
                
                if Rflag == 1
                    R = R_new; % use new radius
                    Omega_max = Vmaxtip / R; % recalculate maximum angular velocity
                    A = pi * R^2;  % recalculate area
                end
                
                % Begin inner loop
                cond2 = 1; % conditional for when to exit inner loop
                
                while cond2 == 1
                    innerLoop_counter = innerLoop_counter + 1;
                    Ec = Ed * Mbatt; % initial guess for energy capacity [W*hr]
                    % If figure of merit is not in the right range, use updated value
                    % of rotational velocity instead
                    if FM == 1
                        Omega = Omega_new;
                    end
                    
                    if R < 3 % if R is too low, increase R
                        Rflag = 1; % allows radius change at the beginning of the loop
                        R_new = R + 0.01*R; % increase radius
                        R = R_new; % use new radius
                        Omega_max = Vmaxtip / R; % recalculate maximum angular velocity
                        A = pi * R^2;  % recalculate area
                    elseif R >10 % if R is too high, decrease radius
                        Rflag = 1; % allows radius change at the beginning of the loop
                        R_new = R - 0.01*R; % decrease radius
                        R = R_new; % use new radius
                        Omega_max = Vmaxtip / R; % recalculate maximum angular velocity
                        A = pi * R^2;  % recalculate area
                    end
                    
                    
                    if Omega < 20 % If omega is too low, increase omega
                        Omega = Omega + 0.01 * Omega;
                    elseif Omega > Omega_max % if Omega is too high (constrained by tip mach #), decrease omega
                        Omega = Omega - 0.01*Omega;
                    end
                    
% % OLD - tip mach # too high - had to eliminate arbitrary range of 20-50
%                     if Omega < 20 % If omega is too low, increase omega
%                         Omega = Omega + 0.01 * Omega;
%                     elseif Omega > 50 % if Omega is too high, decrease omega
%                         Omega = Omega - 0.01*Omega;
%                     end
                    
                    Vtip = Omega * R; % tip velocity [m/s]
                    Ct = T/(A * rho * Vtip^2); % thrust coefficient
                    
                    % Determine Blade Solidity
                    muMax = Vmaxfwd / Vtip; % maximum advance ratio
                    BL = interp1(advRatio, bladeLoading, muMax); % maximum blade loading
                    sigma = Ct / BL; % blade solidity
                    
                    % Determine the chord and the aspect ratio
                    c = (sigma * pi * R)/numBlades; % chord [m]
                    
                    % If figure of merit is not in the right range, use updated value
                    % of chord length instead
                    if FM == 2
                        c = c_new;
                        sigma = c * numBlades/ (pi * R); % recalculate solidity
                        BL = Ct/sigma; % recalculate blade loading
                    end
                    
                    AR = R/c; % aspect ratio
                    
                    % Adjust rotational velocity/radius until aspect ratio is between 15 and 20
                    if AR >= 15 && AR <= 20 && R >=3 && R <= 10 && Omega >= 20 && Omega <=Omega_max
                        cond2 = 0; % exit loop
                    elseif AR > 20
                        % decrease radius and omega
                        Rflag = 1; % allows radius change at the beginning of the loop
                        R_new = R - 0.01*R; % decrease radius
                        R = R_new; % use new radius
                        Omega_max = Vmaxtip / R; % recalculate maximum angular velocity
                        A = pi * R^2;  % recalculate area
                        
                        Omega = Omega - 0.01 * Omega;
                    elseif AR < 15
                        Rflag = 1; % allows radius change at the beginning of the loop
                        R_new = R + 0.01*R; % increase radius
                        R = R_new; % use new radius
                        Omega_max = Vmaxtip / R; % recalculate maximum angular velocity
                        A = pi * R^2;  % recalculate area
                    end
                end
                % End of inner loop
            end
            
            
            % Determine the average lift coefficient
            Cl = 6 * Ct/sigma; % average lift coefficient
            
            % Preliminary Power Calculations
            
            % Estimate of power required to hover
            B = 1 - sqrt(2 * Ct)/numBlades; % main rotor tip loss factor
            Pi_hover = (1/B) * T^(1.5)/sqrt(2 * rho * A); % induced power of main rotor with tip loss [W]
            Po_hover = 0.125 * sigma * Cdo * rho * A * Vtip^3; % profile power of main rotor [W]
            P_hover_new = Pi_hover + Po_hover; % calculated total power of main rotor in hover [W]
            
            % Estimate battery weight
                mb = Ec/Ed; % battery mass [kg]
%             W_battery = mb*9.81 % [N] 
%             
             % Power Density Calculations - Battery needs enough Power to cover highest component of Power 
             
                if mainLoop_counter == 1 %initialize through the first iteration, until Pow_max is calculated later
                    Pow_max = P_hover_new;
                    Pow_max_main = P_hover_new;
                    Pow_max_tail = P_hover_new;
                else 
                    Pow_max = Pow_max;
                end 
             % Estimate Motor and Invertor Weight more advanced

                % Use current technology levels for the test cases (2 for inv. 2.2 for
                % motor SP - similar to how we use 144 
                % if ed = 144, etc. etc. 

            if Ed <= 144 % this constraint will include the test cases 
                SPinv = 2200; %W/kg 
                SPmot = 2000; %W/kg 
            elseif Ed == 250
                SPinv = 9000;
                SPmot = 9000; 
            else 
                SPinv = 19000;
                SPmot = 16000; 
            end 
            % the power that sizes the invertor, is the power flowing into
            % the invertor, this needs to be higher than the mission power
            % due to the losses associated with imperfect efficiencies of
            % the motor and invertor. The power flowing into the motor, is
            % the power through the invertor multiplied by the efficiency
            % of the invertor. The power that sizes the battery needs to be
            % the power flowing into the invertor.
            
            inv_eff = 0.99; % Uranda
            mot_eff = 1 ;
            if mainLoop_counter == 1
                Pmot = Pow_max/mot_eff; % mult by motor efficiency when found 
                Mm = Pmot/SPmot;
                Ptotal_hover = P_hover_new;
                Ptotal_fwd = P_hover_new;
                PtTail_hover = 0;
                PtTail_fwd = 0;
                Pt_hover = 0;
                Pt_fwd = 0;
            else
                if max([Pt_hover Pt_fwd]) == Pt_hover
                    Mm_main = Pt_hover/SPmot;
                elseif max([Pt_hover Pt_fwd]) == Pt_fwd
                    Mm_main = Pt_fwd/SPmot;
                end
                
                if max([PtTail_hover PtTail_fwd]) == PtTail_hover
                    Mm_tail = PtTail_hover/SPmot;
                elseif max([PtTail_hover PtTail_fwd]) == PtTail_fwd
                    Mm_tail = PtTail_fwd/SPmot;
                end
                
                Mm = Mm_tail + Mm_main;
            end
                    
%             Pmot = Pow_max/mot_eff; % mult by motor efficiency when found 
%             Mm = Pmot/SPmot; 
            
            Pinv = (Pow_max_main + Pow_max_tail)/(inv_eff*mot_eff);% determines the power of the invertor as the maximum power needed during the mission
            Minv = Pinv/SPinv;
      
            W_propulsion = (Mm + Minv) * 9.81 ; 
            
            Pow_dens = 3 * Ed; % W/kg
            bm_Ed = mb;
            
            bm_pw = Pinv/Pow_dens;
                
%                 % If energy capacity used in battery weight calculations has
%                 % converged, exit the outer loop
                           if bm_pw > bm_Ed
                               bm = bm_pw;
                           else 
                               bm = bm_Ed;
                           end 
            
               W_battery = bm * 9.81; % Weight of the battery [N]


      
            % Estimate motor weight
%             Mm = P_hover_new/5200; % motor mass [kg] (assumption from Siemans motor - 260 kW/50 kg)
%             W_propulsion = Mm * 9.81; % motor weight [N]
            


            % Make second gross weight estimate
            
            % Attempt 1
            %     W_blades1 = .2247*(0.06 * (We - W_battery) * R^0.4 * sigma^0.33); % [N]
            %     W_hub1 = .2247*(0.0135 * (We - W-battery) * R^0.42); % [N]
            %     W_fuselage = 0.21 * (We - W_battery); % [N]
            %     W_controls = 0.06 * (We - W_battery); % [N]
            %     W_electrical = 0.06 * (We - W_battery); % [N]
            %     W_fixedequip = 0.28 * (We - W_battery); % [N]
            
            % Attempt 2
            %     W_mainrotor = 1.7 * (Wg * 0.22481)^(0.342) * (R * 3.28)^1.58 * sigma^0.63; % [lbs]
            %     W_tailrotor = 7.12 * (Wg * 0.22481/1000)^0.446 * (0.2*R*3.28)^1.62 * (2*sigma)^0.66; % [lbs]
            %     W_flightcontrol = 0.0226 * (Wg * 0.22481)^0.712 * (Vfwd/0.5144)^0.653; % [lbs]
            %     W_landinggear = 0.0475 * (Wg * 0.22481/1000)^0.975; % [lbs]
            %     W_fuselage = 0.37 * (Wg * 0.22481)^0.598 * (R *3.28)^.942; % [lbs]
            
            % Attempt 3 - NDARC AFDD Weight Correlations
            x = 1; % technology factor
            
            numRotor = 1; % number of rotors
            vblade = 1.25; % flap natural frequency [per rev] (Johnson)
            W_blades = 0.02606 * numRotor * numBlades^0.6592 * (R*3.28)^1.3371 * (c*3.28)^0.9959 * (Vtip*3.28)^0.6682 * vblade^2.5279; % [lbs]
            
            vhub = 1; % flap natural freqency [per rev] (assumption)
            W_hub = 0.003722 * numRotor * numBlades^0.2807 * (R*3.28)^1.5377 * (Vtip*3.28)^0.4290 * vhub^2.1414 * (W_blades/numRotor)^0.5505; % [lbs]
            
            Rtr = 1.3 * sqrt(Wg * 0.2247/1000); % [ft] radius of tail rotor (from Kee)
            Pds = P_hover_new * 0.00134102/0.9; % drive system rated power [hp] (used total hover power - put in efficiency?)
            W_tailrotor = 1.3778 * Rtr^0.0897 * (Pds * R *3.28/(Vtip*3.28))^0.8951; % [lbs]
            
            f_ramp = 1; % assume 1 for no cargo ramp
            nz = 1.5; % design ultimate load factor at structural design gross weight
            Sbody = interp1(flatPlateWeightData, flatPlateAreaData, (Wg * 0.2247)); % wetted area of the body [ft^2] (Leishman equivalent flat plate area)
            len = 1.25 * (R * 3.28 + Rtr + 0.5); % length of the fuselage [ft] - use calculated one from Kee plus 25% (assumption)
            W_fuselage = 5.896 * f_ramp * (Wg * 0.22481/1000)^0.4908 * nz^0.1323 * Sbody^0.2544 * len^0.6100; % [lbs]
            
            % HT and VT sizing based on tail rotor radius (page 227 Weisner 1974)
            cHT = (0.8 + 0.48)/2 * Rtr; % mean chord of HT [ft]
            bHT = 1.25 * 2 * Rtr; % span of HT [ft]
            sHT = cHT * bHT; % planform area of HT [ft^2]
            AR_HT = bHT^2/sHT; % aspect ratio of HT
            W_HT = 0.7176 * sHT^1.1881 * AR_HT^0.3173; % weight of horizontal tail [lbs]
            
            cVT = (0.47 + 0.94)/2 * Rtr; % mean chord of VT [ft]
            bVT = 1.05 * Rtr; % span of VT [ft]
            sVT = cVT * bVT; % planform area of VT
            AR_VT = bVT^2/sVT; % aspect ratio of VT
            ftr = 1.6311; % for tail rotor located on vertical tail (=1 if not)
            W_VT = 1.0460 * ftr * sVT^0.9441 * AR_VT^0.5332; % weight of vertical tail [lbs]
            
            fLG = 0.0325; % landing gear weight fraction
            W_landinggear = fLG * Wg * 0.2247; % landing gear weight [lbs]
            
            Omega_eng = (Omega * 9.549) * 6; % engine output speed [rpm] (ratio ranges from 6-1 to 9-1)
            frs = 0.13; % rotor shaft weight fraction
            W_gearbox = (1 - frs) * 95.7634 * numRotor^0.38553 * Pds^0.78137 * Omega_eng^0.09899/(Omega * 9.549)^0.80686; % weight of gearbox [lbs]
            W_rotorshaft = frs * 95.7634 * numRotor^0.38553 * Pds^0.78137 * Omega_eng^0.09899/(Omega * 9.549)^0.80686; % weight of rotor shaft [lbs]
            
            Qds = Pds / (Omega * 9.549); % [hp/rpm]
            xhub = R * 3.28 + Rtr; % length of drive shaft between rotors [ft] !!guess - added rotor radii!!
            Nds = 1; % number of intermediate drive shafts !!guess!!
            fP = 0.15; % rotor rated power fraction (for single main rotor and tail rotor)
            W_driveshaft = 1.166 * Qds^0.3828 * xhub^1.0455 * Nds^0.3909 * (0.01*fP)^0.2693; % driveshaft weight [lbs]
            
            W_rotorbrakes = 0.000871 * W_blades * (0.01 * Vtip * 3.28)^2; % rotor brake weight [lbs]
            
            fnbsv = 1; % (=1.8984 for ballistically survivable, 1 otherwise)
            W_controls = 2.1785 * fnbsv * (Wg * 0.2247)^0.3999 * numRotor^1.3855; % flight controls weight [lbs]
            
            We_new =  x * (W_blades + W_hub + W_tailrotor + W_fuselage + W_HT + W_VT + W_landinggear + W_rotorshaft + W_rotorbrakes + W_controls) * 4.45 + 1.1*(W_battery + W_propulsion); % [N]

            Wg_new = We_new + (payload * 4.45); % new estimate of gross weight [N]
            
            Cp_ideal = Ct^(3/2)/sqrt(2);
            Cp = P_hover_new/(rho * A *Vtip^3);
            
            figureMerit = Cp_ideal/Cp;
            %     figureMerit = Pi/Pt_new; % Figure of merit
            
            % If total power and gross weight have converged within 10%, adjust
            % geometry to get figure of merit between 0.7 and 0.8
            if abs(P_hover_new - P_hover)/P_hover < 0.01 && abs(Wg_new - Wg)/Wg < 0.01
                if figureMerit < 0.7
                    FM = 1; % first case for figure of merit adjustments
                    % reduce rotational velocity
                    skipInitial = 0;
                    Omega_new = Omega - 0.01 * Omega;
                    
                    % Recalculate hover power/weight using adjustments
                    P_hover = P_hover_new;
                    Wg = Wg_new;
                    We = We_new;
                    powerCalc = 0;
                    %         elseif figureMerit > 0.8
                    %             skipInitial = 0;
                    %             FM = 2; % second case for figure of merit adjustments
                    %             % increase the chord length
                    %             c_new = c + 0.01*c;
                    %
                    %             % Recalculate hover power/weight using adjustments
                    %             Pt = Pt_new;
                    %             Wg = Wg_new;
                    %             We = We_new;
                    %             powerCalc = 0;
                else
                    powerCalc = 1; % move to final power calculations if power/weight have converged and figure of merit is within specifications
                    T = Wg_new; % set thrust equal to new gross weight
                end
            else
                FM = 0;
                skipInitial = 0;
                % Recalculate hover power using new gross weight estimate
                P_hover = P_hover_new;
                Wg = Wg_new;
                We = We_new;
                powerCalc = 0;
            end
            
            if powerCalc == 1
                
                % Determine the power required to hover
                PP_OGE = 0.708 * (h/(R*2))^3 - 1.4569 * (h/(R*2))^2 + 1.3432 * (h/(R*2)) + 0.5147; % Power/Power out of ground effect
                Pi_GE = PP_OGE * Pi_hover; % induced power of rotor with tiploss in ground effect [W]
                Pt_hover = Pi_GE + Po_hover; % total hover power [W]
                
                % Determine the parasite power required in forward flight
                EFPA = interp1(flatPlateWeightData, flatPlateAreaData, (Wg_new * 0.2247)) * 0.0929; % equivalent flat plate area [m^2] (interpolated from relation in Leishman)
                Pp_fwd = 0.5 * rho * Vfwd^3 * EFPA; % parasite power in forward flight[W]
                
                D = 0.5 * rho * Vfwd^2 * EFPA; % [N]
                T = sqrt(Wg_new^2 + D^2);
                
                % Determine main rotor power required and mach number of advancing blade
                % tip for forward flight
                mu = Vfwd/Vtip; % advance ratio
                Po_fwd = (1 + 4.3 * mu^2) * Po_hover; % profile power in forward flight [W]
                Vi_fwd = T/(2 * rho * A * Vfwd); % induced velocity in forward flight [m/s] (equation from Leishman)
                
                Pi_fwd = (1.25) * T * Vi_fwd; % induced power in forward flight [W]
                Pt_fwd = Pp_fwd + Po_fwd + Pi_fwd; % total power in forward flight [W]
                Mtip = (Vfwd + Vtip)/a; % Tip Mach number
                
                % Tail Rotor Design
                RTail = 1.3 * sqrt(Wg * 0.2247/1000) * 0.3048; % tail rotor radius [m]
                OmegaTail = 4.5 * Omega; % tail rotor rotational velocity [rad/s]
                CdoTail = 0.0138 * Cdo; % tail rotor drag coefficient
                numBladesTail = 4; % number of blades on tail rotor (assumption based on Leishman)
                L = R + RTail + 0.5 * 0.3048; % length of the fuselage from the center of gravity to the tail rotor hub [m]
                ARTail = 6.25; % aspect ratio of the tail (average of range 4.5-8 given in Kee)
                cTail = RTail/ARTail; % chord of the tail [m]
                sigmaTail = numBladesTail * cTail/(pi * RTail); % tail rotor aspect ratio
                
                TTail = Pt_hover/(Omega * L); % thrust of tail rotor [N]
                ATail = pi * RTail^2; % area of tail rotor [m^2]
                VtipTail = OmegaTail * RTail; % tip velocity of tail rotor
                CtTail = TTail/(ATail * rho * VtipTail^2); % thrust coefficient of tail rotor
                BTail = 1 - sqrt(2 * CtTail)/numBladesTail; % tail rotor tiploss factor
                PiTail_hover = (1/BTail) * TTail^1.5/sqrt(2 * rho * ATail); % induced power of tail rotor in hover [W]
                PoTail_hover = 0.125 * sigmaTail * CdoTail * rho * ATail * VtipTail^3; % profile power of tail rotor in hover [W]
                PtTail_hover = PiTail_hover + PoTail_hover; % total power of tail rotor in hover [W]
                
                muTail = Vfwd/VtipTail; % advance ratio of the tail rotor
                PoTail_fwd = PoTail_hover * (1 + 4.3 *muTail^2); % profile power of tail rotor in forward flight [W]
                TTail_fwd = Pt_fwd / (Omega * L); % thrust of tail rotor in forward flight [N]
                CtTail_fwd = TTail_fwd/(ATail * rho * VtipTail^2); % thrust coefficient of tail in forward flight
                BTail_fwd = 1 - sqrt(2 * CtTail_fwd)/numBladesTail; % tail rotor tiploss factor in forward flight
                ViTail_fwd = TTail_fwd/(2 * rho * ATail * Vfwd); % induced velocity of tail rotor in forward flight [m/s]
                PiTail_fwd = (1/BTail_fwd) * TTail_fwd * ViTail_fwd; % induced power of tail rotor in forward flight [W]
                PtTail_fwd = PoTail_fwd + PiTail_fwd; % total power of tail rotor in forward flight [W]
                MtipTail = (Vfwd + VtipTail)/a; % tip Mach number of tail rotor
                
                % Calculate Climb Power
                Pc = T * rateClimb; % climb power [W]
                
                % Recalculate total powers with tail rotor included
                Ptotal_hover = (Pt_hover + PtTail_hover);
                Ptotal_fwd = (Pt_fwd + PtTail_fwd);
                
                % Calculate energy capacity needed
                Ec_hover = (hoverTime)/3600 * Ptotal_hover; % energy capacity required for hover [W*h]
                Ec_cruise = (cruiseTime+reserve)/3600 * Ptotal_fwd; % energy capacity required for cruise [W*h]
                Ec_climb = climbTime/3600 * Pc; % energy capacity required to climb [W*h]
                Ec_tot = (Ec_hover + Ec_cruise + Ec_climb)/0.9; % total energy capacity required [W*h] (10% unusable energy - McDonald and German)
                
                Ec_counter = Ec_counter+1;
                % Ecs = [Ecs, Ec_tot];
                
                Pow_max = max([Pc Ptotal_hover Ptotal_fwd]); % resets the value of Pow_max 
                Pow_max_tail = max([PtTail_hover PtTail_fwd]);
                Pow_max_main = max([Pt_hover Pt_fwd]);
                
                % Recalculate BM based on new energy and mass requirements
                
                bm_Ed = Ec_tot/Ed; % battery mass [kg] 
%             
             % Power Density Calculations - Battery needs enough Power to cover highest component of Power 
                bm_pw = Pinv/Pow_dens;
                
%                 % If energy capacity used in battery weight calculations has
%                 % converged, exit the outer loop
                           if bm_pw > bm_Ed
                               bm = bm_pw;
%                                powcount(i,j) = 1 ;
                           else 
                               bm = bm_Ed;
%                                powcount(i,j) = 0 ;
                           end 
%                            bm_ratio(i,j) = bm_Ed/bm_pw;

         
               % If energy capacity used in battery weight calculations has

                % converged, exit the outer loop

% P_max = max([Pc, Ptotal_hover, Ptotal_fwd])

%            energyCSBatt = bm * Ed; %energy of the constraint sized (CS) battery based on either the power constraint or energy constraint (whichever is driving)
%   
%             if abs(energyCSBatt - Ec_tot)/Ec_tot > 0.01
%                     del = Ec_tot - Ec;
%                     Ec = Ec + 0.1 * del;
%                     % Ec = Ec_tot;
%                     % skipInitial = 1;
%                 else
%                     cond1 = 0;
%                 end   
               

                    
                     if abs(bm - Mbatt)/Mbatt > 0.01
                        del = bm - Mbatt;
                        Mbatt = Mbatt + 0.1 * del;
                        % Ec = Ec_tot;
                        % skipInitial = 1;
                    else
                        cond1 = 0;
                    end   

            end
            
        end


        energies(i,j) = Ec_tot/1000;
        weights(i,j) = Wg_new * 0.2247;
        radii(i,j) = R*3.28;
        hoverpowers(i,j) = Ptotal_hover/1000;
tailHover(i) = PtTail_hover;
tailFwd(i) = PtTail_fwd;
mainHover(i,j) = Pt_hover;
mainFwd(i,j) = Pt_fwd;
     end
 end





% grossweights(i) = Wg_new * 0.2247;
% batteries(i) = W_battery * 0.2247;
% radii(i) = R*3.28;
% energycap(i) = Ec_hover/1000;
%
%
% end



%% Outputs
GrossWeightLbs = Wg_new * 0.2247
EmptyWeightLbs = We_new * 0.2247
BatteryWeightLbs = W_battery * 0.2247
TotalEnergyCapacity_kWh = Ec_tot/1000
TotalHoverPower_kW = Ptotal_hover/1000
TotalCruisePower_kW = Ptotal_fwd/1000
Radius_ft = R*3.28
RPM = Omega * 9.549

% iterations = [0:Ec_counter];
% plot(iterations, Ecs)

%% Plots

% PASSENGERS
% figure(1)
% pass1 = 1:10;
% pass2 = 1:10;
% pass3 = 2:10;
% plot(pass1, energies(1:length(pass1), 1), ':k', 'LineWidth', 3);
% hold on
% plot(pass2, energies(1:length(pass2), 2), '--k', 'LineWidth', 3);
% hold on
% plot(pass3, energies(1:length(pass3), 3), 'k', 'LineWidth', 3);
% box off
% set(gcf,'color','w');
% xlabel('Number of Passengers', 'FontSize', 17, 'FontWeight', 'bold')
% ylabel('Total Energy (kWh)', 'FontSize', 17, 'FontWeight', 'bold')
% ylim([0 250])
% set(gca, 'linewidth', 2, 'FontSize', 15)
% 
% % 
% % 
% grossweights = [1000 2000 3000 4000 5000];
% for k = 1:length(grossweights)
%     findWeight = grossweights(k);
%     numAtWeight = interp1(weights(1:length(pass1)), pass1, findWeight);
%     energyAtNum = interp1(pass1, energies(1:length(pass1)), numAtWeight);
% 
%     numAtWeight2 = interp1(weights(1:length(pass2),2), pass2, findWeight);
%     energyAtNum2 = interp1(pass2, energies(1:length(pass2),2), numAtWeight2);
% 
%     numAtWeight3 = interp1(weights(1:length(pass3),3), pass3, findWeight);
%     energyAtNum3 = interp1(pass3, energies(1:length(pass3),3), numAtWeight3);
% 
%     numbers = [numAtWeight numAtWeight2 numAtWeight3];
%     energies2 = [energyAtNum energyAtNum2 energyAtNum3];
% 
%     hold on
%     plot(numbers, energies2, 'k', 'LineWidth', 1.5)
%     text(numbers(3)+0.1, energies2(3)-3, strcat(num2str(findWeight), ' lbs'), 'FontSize', 13);
% end
% text(numbers(2)+0.1, energies2(2)-3, '5000 lbs', 'FontSize', 13);
% leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
% title(leg, 'Battery Energy Density', 'FontSize', 12)
% leg.FontSize = 13;

% % CRUISE SPEED
% figure(2)
% speed1 = 60:120;
% speed2 = 35:120;
% speed3 = 25:120;
% plot(speed1, energies(1:length(speed1), 1), ':k', 'LineWidth', 2)
% hold on
% plot(speed2, energies(1:length(speed2), 2), '--k', 'LineWidth', 2)
% hold on
% plot(speed3, energies(1:length(speed3), 3), 'k', 'LineWidth', 2)
% box off
% set(gcf,'color','w');
% xlabel('Cruise Speed (kt)', 'Color', 'k')
% ylabel('Total Energy (kWh)')
% leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
% title(leg, 'Battery Energy Density')


% % DISTANCE
figure(3)
dist1 = 5:69;
dist2 = 5:183;
dist3 = 5:329;


plot(dist1, energies(1:length(dist1), 1), ':r', 'LineWidth', 2)
hold on
plot(dist2, energies(1:length(dist2), 2), '--r', 'LineWidth', 2)
hold on
plot(dist3, energies(1:length(dist3), 3), 'r', 'LineWidth', 2)
box off
set(gcf,'color','w');

xlabel('Distance (miles)', 'FontSize', 17)
ylabel('Total Energy (kWh)', 'FontSize', 17)
set(gca, 'linewidth', 2, 'FontSize', 15)

xlabel('Distance (miles)')
ylabel('Total Energy (kWh)')

% grossweights = [2000 3000 4000 5000];
% for k = 1:length(grossweights)
%     findWeight = grossweights(k);
%     numAtWeight = interp1(weights(1:length(dist1)), dist1, findWeight);
%     energyAtNum = interp1(dist1, energies(1:length(dist1)), numAtWeight);
% 
%     numAtWeight2 = interp1(weights(16:length(dist2),2), dist2(16:end), findWeight);
%     energyAtNum2 = interp1(dist2, energies(1:length(dist2),2), numAtWeight2);
% 
%     numAtWeight3 = interp1(weights(111:length(dist3),3), dist3(111:end), findWeight);
%     energyAtNum3 = interp1(dist3, energies(1:length(dist3),3), numAtWeight3);
% 
%     numbers = [numAtWeight numAtWeight2 numAtWeight3];
%     energies2 = [energyAtNum energyAtNum2 energyAtNum3];
% 
%     hold on
%     plot(numbers, energies2, 'k')
% 
%     text(numbers(3)+3, energies2(3)-5, strcat(' ', num2str(findWeight), ' lbs'), 'FontSize', 13);
%    
% end 


 leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
 title(leg, 'Battery Energy Density')
 leg.FontSize = 13;

% 
% figure;
% plot(bm_ratio(1:33,1),':k','LineWidth',2)
% hold on
% plot(bm_ratio(1:100,2),'--k','LineWidth',2)
% plot(bm_ratio(:,3),'k','LineWidth',2)
% leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
% xlabel('Distance (miles)', 'FontSize', 17, 'FontWeight', 'bold')
% ylabel('Battery Mass Ratio', 'FontSize', 17, 'FontWeight', 'bold')
% box off
% set(gcf,'color','w');
% set(gca, 'linewidth', 2, 'FontSize', 12)
% leg.FontSize = 10;
% title(leg, 'Battery Energy Density')



% HOVER TIME
% figure(4)
% time1 = 10:10:1220;
% time2 = 10:10:3810;
% time3 = 10:10:6740;
% plot(time1, energies(1:length(time1), 1), ':r', 'LineWidth', 2)
% hold on
% plot(time2, energies(1:length(time2), 2), '--r', 'LineWidth', 2)
% hold on
% plot(time3, energies(1:length(time3), 3), 'r', 'LineWidth', 2)
% box off
% set(gcf,'color','w');
% xlabel('Hover Time (sec)', 'FontSize', 17)
% ylabel('Total Energy (kWh)', 'FontSize', 17)
% set(gca, 'linewidth', 2, 'FontSize', 15)
% 
% % grossweights = [2000 3000 4000 5000];
% % grossweights = 6000;
% % for k = 1:length(grossweights)
% %     findWeight = grossweights;
% %     numAtWeight = interp1(weights(1:length(time1)), time1, findWeight);
% %     energyAtNum = interp1(time1, energies(1:length(time1)), numAtWeight);
% % 
% %     numAtWeight2 = interp1(weights(1:length(time2),2), time2, findWeight);
% %     energyAtNum2 = interp1(time2, energies(1:length(time2),2), numAtWeight2);
% % 
% %     numAtWeight3 = interp1(weights(111:length(time3),3), time3(111:end), findWeight);
% %     energyAtNum3 = interp1(time3, energies(1:length(time3),3), numAtWeight3);
% % 
% %     numbers = [numAtWeight numAtWeight2 numAtWeight3];
% %     energies2 = [energyAtNum energyAtNum2 energyAtNum3];
% % 
% %     hold on
% %     plot(numbers, energies2, 'k', 'LineWidth', 1.5)
% %     text(numbers(3)+3, energies2(3)-5, strcat(' ', num2str(findWeight), ' lbs'), 'FontSize', 13);
% 
% % end
% 
% leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
% title(leg, 'Battery Energy Density')
% leg.FontSize = 13;


% figure(5)
% plot(distances/1609, energies(:,1), 'k', 'LineWidth', 2)
% hold on
% plot(distances/1609, energies(:,2), '-.k', 'LineWidth', 2)
% hold on
% plot(distances/1609, energies(:,3), '--k', 'LineWidth', 2)
% hold on
% plot(distances/1609, energies(:,4), ':k', 'LineWidth', 2)
% box off
% set(gcf,'color','w');
% xlabel('Distance (miles)', 'FontSize', 14)
% ylabel('Total Energy (kWh)', 'FontSize', 14)
% set(gca, 'linewidth', 1.5, 'FontSize', 14)
% leg = legend('2', '4', '6', '8','Location', 'NW');
% title(leg, 'Number of Passengers')
% 
% figure(6)
% plot(distances/1609, weights(:,1),'k', 'LineWidth', 2)
% hold on
% plot(distances/1609, weights(:,2), '-.k', 'LineWidth', 2)
% hold on
% plot(distances/1609, weights(:,3), '--k', 'LineWidth', 2)
% hold on
% plot(distances/1609, weights(:,4), ':k', 'LineWidth', 2)
% box off
% set(gcf,'color','w');
% xlabel('Distance (miles)', 'FontSize', 14)
% ylabel('Gross Weight (lbs)', 'FontSize', 14)
% set(gca, 'linewidth', 1.5, 'FontSize', 14)
% leg = legend('2', '4', '6', '8','Location', 'NW');
% title(leg, 'Number of Passengers')