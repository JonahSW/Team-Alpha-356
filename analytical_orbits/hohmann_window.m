%jjs280
%03/16/2021
%Implements a function that can calculate the departure date window for a hohmann transfer
%Returns an array of possible departure days as both indexes (days from 01/01/2040) and dates

function [departure_days, departure_dates] = hohmann_window(body1,body2, plotme)
    
    %Internally Defined Constants
    g0 = 9.807;%Earth Average gravity [m/s^2]
    mu_sun = 1.32712e11;%Solar Gravitational Parameter [km^3/sec^2]
     
    body1_data = read_ephemeris(body1);
    body2_data = read_ephemeris(body2);
    
    a1 = body1_data(3,:);
    a2 = body2_data(3,:);
    theta1 = body1_data(2,:);
    theta2 = body2_data(2,:);
    
    departure_days = zeros(1,1);
    L12 = zeros(1,1);
    dtheta = zeros(1,1);
    num_windows = 1;
    for i = 1:1:length(a1)
        L12(i) = pi*(1-((a1(i)+a2(i))/(2*a2(i)))^(3/2));
        dtheta(i) = theta1(i) - theta2(i);
        if (dtheta(i) > (L12(i)-L12(i)*0.01)) & (dtheta(i) < (L12(i)+L12(i)*0.01))
            departure_days(num_windows) = i;
            %disp('Found One!');
            num_windows = num_windows+1;
        end
    end
    
    %Find time as dates
    t1 = datetime(2040,1,1,0,0,0);
    for i = 1:1:length(departure_days)
        departure_dates(i) = t1 + days(departure_days(i));
    end
    
    if plotme == 1
        figure()
        plot(1:1:length(theta1),theta1, 1:1:length(theta2),theta2);
        title('True Anomaly Over Time');
        xlabel('Time (days)');
        ylabel('True Anomaly (rad)');
        legend('Body 1','Body 2')
        
        figure()
        plot(1:1:length(dtheta),dtheta);
        hold on
        plot(1:1:length(L12),L12);
        title('Difference in True Anomaly Over Time');
        xlabel('Time (days)');
        legend('Difference in True Anomaly','Hohmann Transfer Alignment')
        
        for i = 1:1:length(departure_dates)
           xline(departure_days(i),'-',datestr(departure_dates(i)),'HandleVisibility','off'); 
        end
    end
end

