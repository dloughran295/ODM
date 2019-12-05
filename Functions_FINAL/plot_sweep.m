function [] = plot_sweep(test_case, heli_type)
% Export special global variables
global energies weights

% Check helicopter type
if strcmp(heli_type, 'compound')
    % Check sweep test case and output plot
    if strcmp(test_case, 'passenger')
        
        % PASSENGERS
        figure(1)
        pass1 = 1:20;
        pass2 = 1:20;
        pass3 = 1:20;
        plot(pass1(1:length(energies(:, 1))), energies(:, 1), ':k', 'LineWidth', 2);
        hold on
        plot(pass2(1:length(energies(:, 2))), energies(:, 2), '--k', 'LineWidth', 2);
        hold on
        plot(pass3(1:length(energies(:, 3))), energies(:, 3), 'k', 'LineWidth', 2);
        box off
        set(gcf,'color','w');
        xlabel('Number of Passengers', 'FontSize', 14)
        ylabel('Total Energy (kWh)', 'FontSize', 14)
        set(gca, 'linewidth', 2, 'FontSize', 12)
        
        grossweights = [1000 2000 3000 4000 5000];
        for k = 1:length(grossweights)
            findWeight = grossweights(k);
            numAtWeight = interp1(weights(:, 1), pass1(1:length(weights(:, 1))), findWeight);
            energyAtNum = interp1(pass1(1:length(energies(:, 1))), energies(:, 1), numAtWeight);
            
            numAtWeight2 = interp1(weights(:, 2), pass2(1:length(weights(:, 2))), findWeight);
            energyAtNum2 = interp1(pass2(1:length(energies(:, 2))), energies(:, 2), numAtWeight2);
            
            numAtWeight3 = interp1(weights(:, 3), pass3(1:length(weights(:, 3))), findWeight);
            energyAtNum3 = interp1(pass3(1:length(energies(:, 3))), energies(:, 3), numAtWeight3);
            
            numbers = [numAtWeight numAtWeight2 numAtWeight3];
            energies2 = [energyAtNum energyAtNum2 energyAtNum3];
            
            hold on
            plot(numbers, energies2, 'k', 'LineWidth', 1.5)
            text(numbers(3)+0.1, energies2(3)-3, strcat(num2str(findWeight), ' lbs'), 'FontSize', 12);
        end
        
        leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
        title(leg, 'Battery Energy Density')
        leg.FontSize = 10;
        
    elseif strcmp(test_case, 'speed')
        
        % CRUISE SPEED
        figure(2)
        speed1 = 28:120;
        speed2 = 25:120;
        speed3 = 25:120;
        plot(speed1, energies(1:length(speed1), 1), ':k', 'LineWidth', 2)
        hold on
        plot(speed2, energies(1:length(speed2), 2), '--k', 'LineWidth', 2)
        hold on
        plot(speed3, energies(1:length(speed3), 3), 'k', 'LineWidth', 2)
        box off
        set(gcf,'color','w');
        xlabel('Cruise Speed (kt)', 'Color', 'k', 'FontSize', 14)
        ylabel('Total Energy (kWh)', 'FontSize', 14)
        set(gca, 'linewidth', 2, 'FontSize', 12)
        leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
        title(leg, 'Battery Energy Density')
        leg.FontSize = 10;
        
    elseif strcmp(test_case, 'distance')
        
        % DISTANCE
        figure(3)
        dist1 = 5:58;
        dist2 = 5:141;
        dist3 = 5:257;
        
        plot(dist1, energies(1:length(dist1), 1), ':r', 'LineWidth', 2)
        hold on
        plot(dist2, energies(1:length(dist2), 2), '--r', 'LineWidth', 2)
        hold on
        plot(dist3, energies(1:length(dist3), 3), 'r', 'LineWidth', 2)
        
        box off
        set(gcf,'color','w');
        xlabel('Distance (miles)', 'FontSize', 14)
        ylabel('Total Energy (kWh)', 'FontSize', 14)
        set(gca, 'linewidth', 2, 'FontSize', 12)
        
        grossweights = [3000 6000 9000 12000 15000];
        
        for k = 1:length(grossweights)
            findWeight = grossweights(k);
            numAtWeight = interp1(weights(1:length(dist1)), dist1, findWeight);
            energyAtNum = interp1(dist1, energies(1:length(dist1)), numAtWeight);
            
            numAtWeight2 = interp1(weights(10:length(dist2),2), dist2(10:end), findWeight);
            energyAtNum2 = interp1(dist2, energies(1:length(dist2),2), numAtWeight2);
            
            numAtWeight3 = interp1(weights(43:length(dist3),3), dist3(43:end), findWeight);
            energyAtNum3 = interp1(dist3, energies(1:length(dist3),3), numAtWeight3);
            
            numbers = [numAtWeight numAtWeight2 numAtWeight3];
            energies2 = [energyAtNum energyAtNum2 energyAtNum3];
            
            hold on
            
            plot(numbers, energies2, 'k', 'LineWidth', 1.5)
            text(numbers(3)+3, energies2(3)-5, strcat(num2str(findWeight), ' lbs'), 'FontSize', 12);
        end
        
        leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
        title(leg, 'Battery Energy Density')
        leg.FontSize = 10;
        
    elseif strcmp(test_case, 'hover')
        
        % HOVER TIME
        figure(4)
        time1 = 10:10:1080;
        time2 = 10:10:3330;
        time3 = 10:10:6150;
        plot(time1, energies(1:length(time1), 1), ':r', 'LineWidth', 2)
        hold on
        plot(time2, energies(1:length(time2), 2), '--r', 'LineWidth', 2)
        hold on
        plot(time3, energies(1:length(time3), 3), 'r', 'LineWidth', 2)
        
        box off
        set(gcf,'color','w');
        xlabel('Distance (miles)', 'FontSize', 14)
        ylabel('Total Energy (kWh)', 'FontSize', 14)
        set(gca, 'linewidth', 2, 'FontSize', 12)
        
        grossweights = [3000 6000 9000 12000 15000];
        
        for k = 1:length(grossweights)
            findWeight = grossweights(k);
            numAtWeight = interp1(weights(1:length(time1)), time1, findWeight);
            energyAtNum = interp1(time1, energies(1:length(time1)), numAtWeight);
            
            numAtWeight2 = interp1(weights(1:length(time2),2), time2, findWeight);
            energyAtNum2 = interp1(time2, energies(1:length(time2),2), numAtWeight2);
            
            numAtWeight3 = interp1(weights([1:99 105:length(time3)],3), time3([1:99 105:end]), findWeight);
            energyAtNum3 = interp1(time3, energies(1:length(time3),3), numAtWeight3);
            
            numbers = [numAtWeight numAtWeight2 numAtWeight3];
            energies2 = [energyAtNum energyAtNum2 energyAtNum3];
            
            hold on
            
            plot(numbers, energies2, 'k', 'LineWidth', 1.5)
            text(numbers(3)+3, energies2(3)-5, strcat(num2str(findWeight), ' lbs'), 'FontSize', 12);
            
            
        end
        
        leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
        title(leg, 'Battery Energy Density')
        leg.FontSize = 10;
        
    end
elseif strcmp(heli_type, 'electric')
    % Check sweep test case and output plot
    if strcmp(test_case, 'passenger')
        
        % PASSENGERS
        figure(1)
        pass1 = 1:10;
        pass2 = 1:10;
        pass3 = 2:10;
        plot(pass1, energies(1:length(pass1), 1), ':k', 'LineWidth', 3);
        hold on
        plot(pass2, energies(1:length(pass2), 2), '--k', 'LineWidth', 3);
        hold on
        plot(pass3, energies(1:length(pass3), 3), 'k', 'LineWidth', 3);
        box off
        set(gcf,'color','w');
        xlabel('Number of Passengers', 'FontSize', 17, 'FontWeight', 'bold')
        ylabel('Total Energy (kWh)', 'FontSize', 17, 'FontWeight', 'bold')
        ylim([0 250])
        set(gca, 'linewidth', 2, 'FontSize', 15)
        
        grossweights = [1000 2000 3000 4000 5000];
        for k = 1:length(grossweights)
            findWeight = grossweights(k);
            numAtWeight = interp1(weights(1:length(pass1)), pass1, findWeight);
            energyAtNum = interp1(pass1, energies(1:length(pass1)), numAtWeight);
            
            numAtWeight2 = interp1(weights(1:length(pass2),2), pass2, findWeight);
            energyAtNum2 = interp1(pass2, energies(1:length(pass2),2), numAtWeight2);
            
            numAtWeight3 = interp1(weights(1:length(pass3),3), pass3, findWeight);
            energyAtNum3 = interp1(pass3, energies(1:length(pass3),3), numAtWeight3);
            
            numbers = [numAtWeight numAtWeight2 numAtWeight3];
            energies2 = [energyAtNum energyAtNum2 energyAtNum3];
            
            hold on
            plot(numbers, energies2, 'k', 'LineWidth', 1.5)
            text(numbers(3)+0.1, energies2(3)-3, strcat(num2str(findWeight), ' lbs'), 'FontSize', 13);
        end
        text(numbers(2)+0.1, energies2(2)-3, '5000 lbs', 'FontSize', 13);
        leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
        title(leg, 'Battery Energy Density', 'FontSize', 12)
        leg.FontSize = 13;
        
    elseif strcmp(test_case, 'speed')
        
        % CRUISE SPEED
        figure(2)
        speed1 = 50:120;
        speed2 = 35:120;
        speed3 = 25:120;
        plot(speed1, energies(1:length(speed1), 1), ':k', 'LineWidth', 2)
        hold on
        plot(speed2, energies(1:length(speed2), 2), '--k', 'LineWidth', 2)
        hold on
        plot(speed3, energies(1:length(speed3), 3), 'k', 'LineWidth', 2)
        box off
        set(gcf,'color','w');
        xlabel('Cruise Speed (kt)', 'Color', 'k')
        ylabel('Total Energy (kWh)')
        leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
        title(leg, 'Battery Energy Density')
        
    elseif strcmp(test_case, 'distance')
        
        % DISTANCE
        figure(3)
        dist1 = 5:61;
        dist2 = 5:170;
        dist3 = 5:309;
        
        
        plot(dist1, energies(1:length(dist1), 1), ':r', 'LineWidth', 3)
        hold on
        plot(dist2, energies(1:length(dist2), 2), '--r', 'LineWidth', 3)
        hold on
        plot(dist3, energies(1:length(dist3), 3), 'r', 'LineWidth', 3)
        box off
        set(gcf,'color','w');
        
        xlabel('Distance (miles)', 'FontSize', 17)
        ylabel('Total Energy (kWh)', 'FontSize', 17)
        set(gca, 'linewidth', 2, 'FontSize', 15)
        
        xlabel('Distance (miles)')
        ylabel('Total Energy (kWh)')
        
        grossweights = [2000 3000 4000 5000];
        for k = 1:length(grossweights)
            findWeight = grossweights(k);
            numAtWeight = interp1(weights(1:length(dist1)), dist1, findWeight);
            energyAtNum = interp1(dist1, energies(1:length(dist1)), numAtWeight);
            
            numAtWeight2 = interp1(weights(16:length(dist2),2), dist2(16:end), findWeight);
            energyAtNum2 = interp1(dist2, energies(1:length(dist2),2), numAtWeight2);
            
            numAtWeight3 = interp1(weights(111:length(dist3),3), dist3(111:end), findWeight);
            energyAtNum3 = interp1(dist3, energies(1:length(dist3),3), numAtWeight3);
            
            numbers = [numAtWeight numAtWeight2 numAtWeight3];
            energies2 = [energyAtNum energyAtNum2 energyAtNum3];
            
            hold on
            plot(numbers, energies2, 'k')
            text(numbers(3)+3, energies2(3)-5, strcat(' ', num2str(findWeight), ' lbs'), 'FontSize', 13);
        end
        
        leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
        title(leg, 'Battery Energy Density')
        leg.FontSize = 13;
        
    elseif strcmp(test_case, 'hover')
        
        % HOVER TIME
        figure(4)
        time1 = 10:10:1060;
        time2 = 10:10:3520;
        time3 = 10:10:6630;
        plot(time1, energies(1:length(time1), 1), ':r', 'LineWidth', 3)
        hold on
        plot(time2, energies(1:length(time2), 2), '--r', 'LineWidth', 3)
        hold on
        plot(time3, energies(1:length(time3), 3), 'r', 'LineWidth', 3)
        box off
        set(gcf,'color','w');
        xlabel('Hover Time (sec)', 'FontSize', 17)
        ylabel('Total Energy (kWh)', 'FontSize', 17)
        set(gca, 'linewidth', 2, 'FontSize', 15)
        
        grossweights = [2000 3000 4000 5000];
        grossweights = 6000;
        for k = 1:length(grossweights)
            findWeight = grossweights;
            numAtWeight = interp1(weights(1:length(time1)), time1, findWeight);
            energyAtNum = interp1(time1, energies(1:length(time1)), numAtWeight);
            
            numAtWeight2 = interp1(weights(1:length(time2),2), time2, findWeight);
            energyAtNum2 = interp1(time2, energies(1:length(time2),2), numAtWeight2);
            
            numAtWeight3 = interp1([weights(1:111, 3); weights(118:length(time3), 3)], [time3(1:111) time3(118:end)], findWeight);
            energyAtNum3 = interp1(time3, energies(1:length(time3),3), numAtWeight3);
            
            numbers = [numAtWeight numAtWeight2 numAtWeight3];
            energies2 = [energyAtNum energyAtNum2 energyAtNum3];
            
            hold on
            plot(numbers, energies2, 'k', 'LineWidth', 1.5)
            text(numbers(3)+3, energies2(3)-5, strcat(' ', num2str(findWeight), ' lbs'), 'FontSize', 13);
            
        end
        
        leg = legend('144 Wh/kg', '250 Wh/kg', '400 Wh/kg', 'Location', 'NW');
        title(leg, 'Battery Energy Density')
        leg.FontSize = 13;
        
    end
end

end