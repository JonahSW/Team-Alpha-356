%jjs280
%03/24/2021
%Implements a function that can calculate the departure date window an orbital transfer given the
%required difference in longitudes
%Returns an array of possible departure days as both indexes (days from 01/01/2040) and dates
%1 is Venus, 2 is Earth, 3 is Mars, 4 is Ceres, 5 in Jupiter
function [departure_days, departure_dates] = transfer_window(body1,body2,L12,plotme)

    body1_data = read_ephemeris(body1);
    body2_data = read_ephemeris(body2);
    
    a1 = body1_data(3,:);
    a2 = body2_data(3,:);
    theta1 = body1_data(2,:);
    theta2 = body2_data(2,:);
    
    departure_days = zeros(1,1);
    dtheta = zeros(1,1);
    num_windows = 1;
    for i = 1:1:length(a1)
        dtheta(i) = theta2(i) - theta1(i);
        if (dtheta(i) > (L12-L12*0.01)) & (dtheta(i) < (L12+L12*0.01))
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