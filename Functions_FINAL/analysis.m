function [] = analysis(test_case, heli_type)

% Import global variables from script
global numPass avgW payload dist Vfwd cruiseTime Vmaxfwd altitude h hoverTime rateClimb climbTime reserve Ed
global alt vSound rho nu advRatio bladeLoading flatPlateWeightData flatPlateAreaData
global Ec_hover Ec_climb Ec_cruise Ec_counter Pc Ptotal_fwd Ptotal_hover A Wg

% Export special global variables
global AR c c_wing Wg_new We_new W_battery Ec_tot Ptotal_hover Ptotal_fwd R Omega energies weights radii hoverpowers

Ed_sweep = [144 250 400];

for j = 1:length(Ed_sweep)
    Ed = Ed_sweep(j);
    
    % Check helicopter type and test case to set variables for each Ed
    if strcmp(heli_type, 'compound')
        % Check test case to set variables for each Ed
        if Ed == 144
            if strcmp(test_case, 'passenger')
                passengers = 1:14;
                sweep = passengers;
            elseif strcmp(test_case, 'speed')
                speeds = [28:120]*.5144;
                sweep = speeds;
            elseif strcmp(test_case, 'distance')
                distances = [5:58]*1609;
                sweep = distances;
            elseif strcmp(test_case, 'hover')
                hovers = 10:10:1080;
                sweep = hovers;
            end
        elseif Ed == 250
            if strcmp(test_case, 'passenger')
                passengers = 1:14;
                sweep = passengers;
            elseif strcmp(test_case, 'speed')
                speeds = [25:120]*.5144;
                sweep = speeds;
            elseif strcmp(test_case, 'distance')
                distances = [5:141]*1609;
                sweep = distances;
            elseif strcmp(test_case, 'hover')
                hovers = 10:10:3330;
                sweep = hovers;
            end
        elseif Ed == 400
            if strcmp(test_case, 'passenger')
                passengers = 1:14;
                sweep = passengers;
            elseif strcmp(test_case, 'speed')
                speeds = [25:120]*.5144;
                sweep = speeds;
            elseif strcmp(test_case, 'distance')
                distances = [5:257]*1609;
                sweep = distances;
            elseif strcmp(test_case, 'hover')
                hovers = 10:10:6150;
                sweep = hovers;
            end
        end
    elseif strcmp(heli_type, 'electric')
        % Check test case to set variables for each Ed
        if Ed == 144
            if strcmp(test_case, 'passenger')
                passengers = 1:10;
                sweep = passengers;
            elseif strcmp(test_case, 'speed')
                speeds = [50:120]*.5144;
                sweep = speeds;
            elseif strcmp(test_case, 'distance')
                distances = [5:61]*1609;
                sweep = distances;
            elseif strcmp(test_case, 'hover')
                hovers = 10:10:1060;
                sweep = hovers;
            end
        elseif Ed == 250
            if strcmp(test_case, 'passenger')
                passengers = 1:10;
                sweep = passengers;
            elseif strcmp(test_case, 'speed')
                speeds = [35:120]*.5144;
                sweep = speeds;
            elseif strcmp(test_case, 'distance')
                distances = [5:170]*1609;
                sweep = distances;
            elseif strcmp(test_case, 'hover')
                hovers = 10:10:3520;
                sweep = hovers;
            end
        elseif Ed == 400
            if strcmp(test_case, 'passenger')
                passengers = 2:10;
                sweep = passengers;
            elseif strcmp(test_case, 'speed')
                speeds = [25:120]*.5144;
                sweep = speeds;
            elseif strcmp(test_case, 'distance')
                distances = [5:309]*1609;
                sweep = distances;
            elseif strcmp(test_case, 'hover')
                hovers = 10:10:6630;
                sweep = hovers;
            end
        end
    end
    
    for i = 1:length(sweep)
        % Check test case again to set other sweep variables
        if strcmp(test_case, 'passenger')
            numPass = passengers(i);
            payload = avgW * numPass;
        elseif strcmp(test_case, 'speed')
            cruiseSpeed = speeds(i);
            Vfwd = cruiseSpeed;
            cruiseTime = dist/Vfwd;
        elseif strcmp(test_case, 'distance')
            dist = distances(i);
            cruiseTime = dist/Vfwd;
        elseif strcmp(test_case, 'hover')
            hoverTime = hovers(i);
        end
        
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
        
        % Calculate wing parameters if helicopter type is compound
        if strcmp(heli_type, 'compound')
            % WING PARAMETERS (added for compounds)
            x_wing = 0.6; % percent of weight the wing lifts in cruise
            Cdo_aircraft = 0.03; % zero-lift drag coefficient of aircraft (dirty fixed gear prop aircraft - Raymer)
            Cdo_wing2 = 0.005; % zero-lift drag coeff of the wing (Assuming NACA 0009 and Re = 3m - Schlicting)
        end
        
        %% Analysis
        
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
        
        
        while cond1 == 1
            mainLoop_counter = mainLoop_counter+1;
            
            if skipInitial == 0
                if01 = 1;
                
                % Determine rotor radius
                R = (0.0011 * (Wg * 0.2247) + 11.496) * 0.3048; % rotor radius [m] (relation determined from data in Leishman)
                
                % Determine maximum rotational velocity
                Omega_max = Vmaxtip / R; % maximum rotational velocity [rad/s]
                
                % Determine thrust coefficient
                T = Wg; % set thrust equal to gross weight [N]
                A = pi * R^2; % disk area [m^2]
                
                if Rflag == 1
                    if02 = 1;
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
                        if03 = 1;
                        Omega = Omega_new;
                    end
                    
                    if R < 3 % if R is too low, increase R
                        if04 = 1;
                        Rflag = 1; % allows radius change at the beginning of the loop
                        R_new = R + 0.01*R; % increase radius
                        R = R_new; % use new radius
                        Omega_max = Vmaxtip / R; % recalculate maximum angular velocity
                        A = pi * R^2;  % recalculate area
                    elseif R >10 % if R is too high, decrease radius
                        if05 = 1;
                        Rflag = 1; % allows radius change at the beginning of the loop
                        R_new = R - 0.01*R; % decrease radius
                        R = R_new; % use new radius
                        Omega_max = Vmaxtip / R; % recalculate maximum angular velocity
                        A = pi * R^2;  % recalculate area
                    end
                    
                    
                    
                    if Omega < 20 % If omega is too low, increase omega
                        if06 = 1;
                        Omega = Omega + 0.01 * Omega;
                    elseif Omega > Omega_max % if Omega is too high, decrease omega
                        if07 = 1;
                        Omega = Omega - 0.01*Omega;
                    end
                    
                    
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
                        if08 = 1;
                        c = c_new;
                        sigma = c * numBlades/ (pi * R); % recalculate solidity
                        BL = Ct/sigma; % recalculate blade loading
                    end
                    
                    AR = R/c; % aspect ratio
                    
                    % Adjust rotational velocity/radius until aspect ratio is between 15 and 20
                    if AR >= 15 && AR <= 20 && R >=3 && R <= 10 && Omega >= 20 && Omega <=Omega_max
                        if09 = 1;
                        cond2 = 0; % exit loop
                    elseif AR > 20
                        if10 = 1;
                        % decrease radius and omega
                        Rflag = 1; % allows radius change at the beginning of the loop
                        R_new = R - 0.01*R; % decrease radius
                        R = R_new; % use new radius
                        Omega_max = Vmaxtip / R; % recalculate maximum angular velocity
                        A = pi * R^2;  % recalculate area
                        
                        Omega = Omega - 0.01 * Omega;
                    elseif AR < 15
                        if11 = 1;
                        Rflag = 1; % allows radius change at the beginning of the loop
                        R_new = R + 0.01*R; % increase radius
                        R = R_new; % use new radius
                        Omega_max = Vmaxtip / R; % recalculate maximum angular velocity
                        A = pi * R^2;  % recalculate area
                    end
                end
                % End of inner loop
            end
            
            % Calculate lift parameters depending on helicopter type
            if strcmp(heli_type, 'compound')
                % WING DESIGN (added for compounds)
                AR_wing = 10; % aspect ratio of the wing (Russell and Johnson)
                lambda = 0.8; % taper ratio of wing (Russell, Silva, Johnson, Yeo)
                e_wing = 0.8; % oswald's efficiency factor of wing (typical for propeller powered aircraft - Raymer)
                s_wing = (x_wing * Wg)/(0.5 * rho * Vfwd^2 * sqrt(pi * AR_wing * e_wing * Cdo_aircraft)); % wing area [m^2](equation from page 136 raymer)
                b_wing = sqrt(AR_wing * s_wing); % wing span [m]
                Cl_wing = x_wing*Wg/(0.5*rho*Vfwd^2*s_wing); % lift coefficient of the wing
                c_wing = s_wing/b_wing; % mean wing chord [m]
                rc_wing = 2 * c_wing/(1 + lambda); % root chord of wing [m]
                tc_wing = rc_wing * lambda; % tip chord of wing [m]
                
                % Keys Paper
                % hover download calculation for flaps deflected 80 deg down
                
                %     step 1: distances
                f_width = 6 * 0.3048; % estimate for width of fuselage (6 ft) [m]
                
                if R*2 <= b_wing
                    if12 = 1;
                    zeta1 = 0;
                else
                    if13 = 1;
                    zeta1 = R - b_wing/2;
                end
                
                zeta2 = R - f_width;
                zeta3 = R;
                
                %step 2: ratio of distance to radius
                zetaR1 = zeta1/R *100;
                zetaR2 = zeta2/R *100;
                zetaR3 = zeta3/R *100;
                
                %step 3/4: kv
                kv_plot = [0 11.85185185 63.7037037 97.77777778 130.3703704 162.962963 197.037037 202.962963 204.4444444 202.962963 208.8888889 214.8148148 225.1851852 241.4814815 262.2222222 280 312.5925926 343.7037037 376.2962963 392.5925926 410.3703704 411.8518519];
                zetaR_plot = [0 17.1957672 32.07010582 40.43650794 49.15343915 60.01322751 72.66534392 83.76322751 97.6984127 109.8346561 120.5753968 128.8161376 134.5767196 141.4351852 148.3134921 154.4642857 161.7526455 170.462963 177.7513228 183.1812169 192.1891534 201.8386243];
                
                kv1 = interp1(zetaR_plot, kv_plot, zetaR1);
                kv2 = interp1(zetaR_plot, kv_plot, zetaR2);
                kv3 = interp1(zetaR_plot, kv_plot, zetaR3);
                
                %step 5: delta kv
                del_kv1 = kv2 - kv1;
                del_kv2 = kv3 - kv2;
                
                %step 6: Cdv
                Cdv = 1.25;
                
                %step 7: section width
                w1 = tc_wing;
                w2 = rc_wing;
                
                %step 8: calc
                DvT1 = 2* (Cdv * w1 /(4*pi*R) * del_kv1);
                DvT2 = 2* (Cdv * w2 /(4*pi*R) * del_kv2);
                
                wingDvT = DvT1 + DvT2; % percentage increase in thrust due to download
            elseif strcmp(heli_type, 'electric')
                % Determine the average lift coefficient
                Cl = 6 * Ct/sigma; % average lift coefficient
            end
            
            % Preliminary Power Calculations depending on helicopter type
            if strcmp(heli_type, 'compound')
                % Estimate of power required to hover (changed for compounds)
                B = 1 - sqrt(2 * Ct)/numBlades; % main rotor tip loss factor
                Pi_hover = (1/B) * (T*(1 + wingDvT/100))^(1.5)/sqrt(2 * rho * A); % changed to add download increase to thrust
                Po_hover = 0.125 * sigma * Cdo * rho * A * Vtip^3; % profile power of main rotor [W]
                P_hover_new = Pi_hover + Po_hover; % calculated total power of main rotor in hover
            elseif strcmp(heli_type, 'electric')
                % Estimate of power required to hover
                B = 1 - sqrt(2 * Ct)/numBlades; % main rotor tip loss factor
                Pi_hover = (1/B) * T^(1.5)/sqrt(2 * rho * A); % induced power of main rotor with tip loss [W]
                Po_hover = 0.125 * sigma * Cdo * rho * A * Vtip^3; % profile power of main rotor [W]
                P_hover_new = Pi_hover + Po_hover; % calculated total power of main rotor in hover [W]
            end
            
            % Estimate battery weight
            mb = Ec/Ed; % battery mass [kg]
            
            % Power Density Calculations - Battery needs enough Power to cover highest component of Power
            
            if mainLoop_counter == 1 %initialize through the first iteration, until Pow_max is calculated later
                if14 = 1;
                Pow_max = P_hover_new;
            else
                if15 = 1;
                Pow_max = Pow_max;
            end
            % Estimate Motor and Invertor Weight more advanced
            
            % Use current technology levels for the test cases (2 for inv. 2.2 for
            % motor SP - similar to how we use 144
            % if ed = 144, etc. etc.
            
            if Ed <= 144 % this constraint will include the test cases
                if16 = 1;
                SPinv = 2200; %W/kg
                SPmot = 2000; %W/kg
            elseif Ed == 250
                if17 = 1;
                SPinv = 9000;
                SPmot = 9000;
            else
                if18 = 1;
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
            else
                if max([Ptotal_hover Ptotal_fwd]) == Ptotal_hover
                    Mm_tail = PtTail_hover/SPmot;
                    Mm_main = Pt_hover/SPmot;
                elseif max([Ptotal_hover Ptotal_fwd]) == Ptotal_fwd
                    Mm_tail = PtTail_fwd/SPmot;
                    Mm_main = Pt_fwd/SPmot;
                end
                Mm = Mm_tail + Mm_main;
            end
            
            Pinv = Pow_max/(inv_eff*mot_eff);% determines the power of the invertor as the maximum power needed during the mission
            Minv = Pinv/SPinv;
            
            W_propulsion = (Mm + Minv) * 9.81 ;
            
            Pow_dens = 3 * Ed; % W/kg
            bm_Ed = mb;
            
            bm_pw = Pinv/Pow_dens;
            
            %                 % If energy capacity used in battery weight calculations has
            %                 % converged, exit the outer loop
            if bm_pw > bm_Ed
                if19 = 1;
                bm = bm_pw;
            else
                if20 = 1;
                bm = bm_Ed;
            end
            
            W_battery = bm * 9.81; % Weight of the battery [N]
            
            
            % NDARC AFDD Weight Correlations
            
            
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
            
            
            % Attempt 3
            We_new =  x * (W_blades + W_hub + W_tailrotor + W_fuselage + W_HT + W_VT + W_landinggear +  W_gearbox + W_rotorshaft + W_driveshaft + W_rotorbrakes + W_controls) * 4.45 + 1.1*(W_battery + W_propulsion); % [N]
            
            % Calculate helicopter weights depending on helicopter type
            if strcmp(heli_type, 'compound')
                % WING WEIGHT (added for compounds)
                fLGloc = 1; % 1 if landing gear is not on wing
                sweepAngle = 0;
                tau = 0.12; % wing airfoil thickness to chord ratio
                bfold = 0; % fraction wing span that folds
                W_wing = 5.66411 * fLGloc * (Wg * 0.2247/(1000 * cos(sweepAngle)))^0.847 * nz^0.39579 * (s_wing * 10.7639)^0.21754 * AR_wing^0.50016 * ((1+lambda)/tau)^0.09359 * (1-bfold)^-0.14356; % wing weight [lbs]
                
                % New empty weight
                We_new =  x * (W_blades + W_hub + W_tailrotor + W_fuselage + W_HT + W_VT + W_landinggear + W_rotorshaft + W_rotorbrakes + W_controls + W_wing) * 4.45 + 1.1*(W_battery + W_propulsion); % [N]
            end
            
            Wg_new = We_new + (payload * 4.45); % new estimate of gross weight [N]
            
            % If compound helicopter, calculate Ct
            if strcmp(heli_type, 'compound')
                Ct = (T*(1 + wingDvT/100))/(A * rho * Vtip^2); % changed to add download increase to thrust
            end
            
            Cp_ideal = Ct^(3/2)/sqrt(2);
            Cp = P_hover_new/(rho * A *Vtip^3);
            
            figureMerit = Cp_ideal/Cp;
            
            
            % If total power and gross weight have converged within 10%, adjust
            % geometry to get figure of merit above 0.7
            if abs(P_hover_new - P_hover)/P_hover < 0.01 && abs(Wg_new - Wg)/Wg < 0.01
                if21 = 1;
                if figureMerit < 0.7
                    if22 = 1;
                    FM = 1; % first case for figure of merit adjustments
                    % reduce rotational velocity
                    skipInitial = 0;
                    Omega_new = Omega - 0.01 * Omega;
                    
                    % Recalculate hover power/weight using adjustments
                    P_hover = P_hover_new;
                    Wg = Wg_new;
                    We = We_new;
                    powerCalc = 0;
                    
                    
                else
                    if23 = 1;
                    powerCalc = 1; % move to final power calculations if power/weight have converged and figure of merit is within specifications
                    T = Wg_new; % set thrust equal to new gross weight
                end
            else
                if24 = 1;
                FM = 0;
                skipInitial = 0;
                % Recalculate hover power using new gross weight estimate
                P_hover = P_hover_new;
                Wg = Wg_new;
                We = We_new;
                powerCalc = 0;
            end
            
            if powerCalc == 1
                if25 = 1;
                
                % Run specific power calculations depending on helicopter
                % type
                if strcmp(heli_type, 'compound')
                    % Determine the power required to hover (changed for compounds - based on Pi_hover)
                    PP_OGE = 0.708 * (h/(R*2))^3 - 1.4569 * (h/(R*2))^2 + 1.3432 * (h/(R*2)) + 0.5147; % Power/Power out of ground effect
                    Pi_GE = PP_OGE * Pi_hover; % induced power of rotor with tiploss in ground effect [W]
                    Pt_hover = Pi_GE + Po_hover; % total hover power [W]
                    
                    % Determine the parasite power required in forward flight (CHANGED for compounds)
                    EFPA = interp1(flatPlateWeightData, flatPlateAreaData, (Wg_new * 0.2247)) * 0.0929; % equivalent flat plate area [m^2] (interpolated from relation in Leishman)
                    e_span = 0.95; % span efficiency (assuming taper ratio of 0.8 and aspect ratio of 10 - Pope)
                    
                    Pp_fwd = (0.5 * rho * Vfwd^3 * (EFPA + s_wing*(Cdo_wing2 + Cl_wing ^ 2 / (pi * e_span * AR_wing) ) )); % parasite power in forward flight [W] (wing drag component added)
                    D = ((0.5 * rho * Vfwd^2 * EFPA) + (0.5 * rho * Vfwd^2 * s_wing * (Cdo_wing2+Cl_wing^2/(pi*e_span*AR_wing)))); % total drag: rotor + wing [N]
                    
                    T = sqrt((Wg_new - x_wing*Wg_new)^2 + D^2); % thrust updated with percent weight offload
                elseif strcmp(heli_type, 'electric')
                    % Determine the power required to hover
                    PP_OGE = 0.708 * (h/(R*2))^3 - 1.4569 * (h/(R*2))^2 + 1.3432 * (h/(R*2)) + 0.5147; % Power/Power out of ground effect
                    Pi_GE = PP_OGE * Pi_hover; % induced power of rotor with tiploss in ground effect [W]
                    Pt_hover = Pi_GE + Po_hover; % total hover power [W]
                    
                    % Determine the parasite power required in forward flight
                    EFPA = interp1(flatPlateWeightData, flatPlateAreaData, (Wg_new * 0.2247)) * 0.0929; % equivalent flat plate area [m^2] (interpolated from relation in Leishman)
                    Pp_fwd = 0.5 * rho * Vfwd^3 * EFPA; % parasite power in forward flight[W]
                    
                    D = 0.5 * rho * Vfwd^2 * EFPA; % [N]
                    T = sqrt(Wg_new^2 + D^2);
                end
                
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
                PiTail_hover = (1/BTail) * TTail^1.5/sqrt(2 * rho * ATail); % induced power of tail rotor in hover [W] (CHANGE?)
                PoTail_hover = 0.125 * sigmaTail * CdoTail * rho * ATail * VtipTail^3; % profile power of tail rotor in hover [W] (CHANGE?)
                PtTail_hover = PiTail_hover + PoTail_hover; % total power of tail rotor in hover [W]
                
                muTail = Vfwd/VtipTail; % advance ratio of the tail rotor
                PoTail_fwd = PoTail_hover * (1 + 4.3 *muTail^2); % profile power of tail rotor in forward flight [W] (CHANGE?)
                TTail_fwd = Pt_fwd / (Omega * L); % thrust of tail rotor in forward flight [N]
                CtTail_fwd = TTail_fwd/(ATail * rho * VtipTail^2); % thrust coefficient of tail in forward flight
                BTail_fwd = 1 - sqrt(2 * CtTail_fwd)/numBladesTail; % tail rotor tiploss factor in forward flight
                ViTail_fwd = TTail_fwd/(2 * rho * ATail * Vfwd); % induced velocity of tail rotor in forward flight [m/s]
                PiTail_fwd = (1/BTail_fwd) * TTail_fwd * ViTail_fwd; % induced power of tail rotor in forward flight [W] (CHANGE?)
                PtTail_fwd = PoTail_fwd + PiTail_fwd; % total power of tail rotor in forward flight [W]
                MtipTail = (Vfwd + VtipTail)/a; % tip Mach number of tail rotor
                
                % Calculate climb power depending on helicopter type
                if strcmp(heli_type, 'compound')
                    % Calculate Climb Power
                    Pc = (T - .0142*T) * rateClimb; % climb power [W] (CHANGED - 1.42% less gross weight - Keys p156)
                elseif strcmp(heli_type, 'electric')
                    % Calculate Climb Power
                    Pc = T * rateClimb; % climb power [W]
                end
                
                % Recalculate total powers with tail rotor included
                Ptotal_hover = (Pt_hover + PtTail_hover);
                Ptotal_fwd = (Pt_fwd + PtTail_fwd);
                
                % Calculate energy capacity needed
                Ec_hover = (hoverTime)/3600 * Ptotal_hover; % energy capacity required for hover [W*h]
                Ec_cruise = (cruiseTime+reserve)/3600 * Ptotal_fwd; % energy capacity required for cruise [W*h]
                Ec_climb = climbTime/3600 * Pc; % energy capacity required to climb [W*h]
                Ec_tot = (Ec_hover + Ec_cruise + Ec_climb)/0.9; % total energy capacity required [W*h] (10% unusable energy - McDonald and German)
                
                Ec_counter = Ec_counter+1;
                
                Pow_max = max([Pc Ptotal_hover Ptotal_fwd]); % resets the value of Pow_max
                
                % Recalculate BM based on new energy and mass requirements
                
                bm_Ed = Ec_tot/Ed; % battery mass [kg]
                %
                % Power Density Calculations - Battery needs enough Power to cover highest component of Power
                bm_pw = Pinv/Pow_dens;
                
                %                 % If energy capacity used in battery weight calculations has
                %                 % converged, exit the outer loop
                if bm_pw > bm_Ed
                    if26 = 1;
                    bm = bm_pw;
                    %powcal(i,j)=1;
                else
                    if27 = 1;
                    bm = bm_Ed;
                    %powcal(i,j) = 0;
                end
                
                if abs(bm - Mbatt)/Mbatt > 0.01
                    if28 = 1;
                    del = bm - Mbatt;
                    Mbatt = Mbatt + 0.1 * del;
                    % Ec = Ec_tot;
                    % skipInitial = 1;
                else
                    if29 = 1;
                    cond1 = 0;
                end
            end
        end
        
        energies(i,j) = Ec_tot/1000;
        weights(i,j) = Wg_new * 0.2247;
        radii(i,j) = R*3.28;
        hoverpowers(i,j) = Ptotal_hover/1000;
        
    end
end
plot_sweep(test_case, heli_type)
end