%function EGA215_954313_A2_Q3 
    clc
    
    tspan = 25200;
    
    [t , p] = ode45(@Q2, [0,tspan], [-1.497449337113264e+04; -2.436887970840462e+02; 0; 0.099407653982047; (-4.073387410672351 - 2); 0]);

    R = 0; %intital value for radius to start loop
    V = 0; %intital value for velocity to start loop
    y1 = 0; %intital value for y-coordinate to start loop
    x1 = 0; %intital value for x-coordinate to start loop
    for i = 1:length(p)
        x1(i) = p(i,1);
        y1(i) = p(i,2);
        x = p(i,1);
        y = p(i,2);
        z = p(i,3);
        u = p(i,4);
        v = p(i,5);
        w = p(i,6);
        R(i) = sqrt((x^2)+(y^2)+(z^2));
        V(i) = sqrt((u^2)+(v^2)+(w^2));
        R_max = max(R);
        R_min = min(R);
       
    end
    t_max = t(22,1); %Time at max altitude
    t_min = t(1,1); %Time at min altiude
    V_max = V(1,22); %Velocity at max altitude
    V_min = V(1,1); %Velocity at min altitude
    x1_max = sqrt(R_max^2 - ((y1(1,1))^2)); %using pythagoras to work out x-coordinate
    y1_max = sqrt(R_max^2 - ((x1(1,22))^2)); %using pythagoras to work out y-coordinate
    
    figure(7)
    hold on
    plot(p(:,1),p(:,2)) %Satellite rajectory
%     title('Satellites Trajectory'), xlabel('Displacement relative to Earth (km)'), ylabel('Displacement relative to Earth (km)')
%     plot(0, -R_min, '*r') %Coordinates of minimum altitude
%     plot(y1_max, x1_max, '*c') %Coordinates of maximum altitude
%     plot(0,0,'*g') %Coordinates on earth
    
    legend('Orbit','Minima','Maxima','Earths position')  
    
    fprintf('Maximum altitude reached: %f Km. Velocity at Max altitude: %f Km/s\nTime at this position: %f hours \n\n',R_max, V_max, t_max)
    fprintf('Minimum altitude reached: %f Km. Velocity at Min altitude: %f Km/s\nTime at this position: %f hours \n', R_min, V_min, t_min)
    
 function dp = Q2(t,p)

    muy = 398600;    
    x = p(1);
    y = p(2);
    z = p(3);
    u = p(4);
    v = p(5);
    w = p(6);
    
    r = sqrt(x^2+y^2+z^2); 
    
    dx = u;
    dy = v;
    dz = w;
    du = -(muy*x)/(r^3);
    dv = -(muy*y)/(r^3);
    dw = -(muy*z)/(r^3);
         
    dp = [dx; dy; dz; du; dv; dw];   
 end




%end