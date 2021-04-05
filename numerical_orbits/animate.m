%jjs280
%03/31/2021
%Function for animating orbits

%In Progress
function animate(position_vector)

    figure()
    
    polarplot(position_vector(1,1),position_vector(2,1));
    
    hold on
    
    for i = 2:1:length(position_vector)
       polarplot(position_vector(1,i),position_vector(2,i),'o');       
    end
    
    hold off
end