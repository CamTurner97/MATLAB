%function EGA215_954313_A1_Q4
    clc  

    Isp = 350;      %specific impulse
    epsilon = 0.15; %structural ratio
    pi_PL = 0.05;   %payload fraction
    g = 9.81;       %gravity constant
   
    vbo_infinity = Isp * g * (1-epsilon) * log(1/pi_PL);%burnout velocity with infinite stages
    vbo = 0.98*vbo_infinity; % 98% of burnout velocity with infinite stages

    fprintf('Velocity for an infinite number of stages: %f m/s\nVelocity at 98 percent infinite number of stages: %f m/s\n',vbo_infinity, vbo)

    v = 0;  %intital velocity 
    N = 1;  %first number of stages

        while v < vbo
            v = Isp * g * N * log((1/( (pi_PL^(1/N))*(1-epsilon) + epsilon)));
            N = N + 1;
            hold on
            stairs(N,v, 'r*') 
        end
     
    fprintf('At rocket stage %f the velocity goes over %f m/s and is %f m/s\n', N , vbo, v)
    fprintf('12 stages of the rocket are needed\n')
    v = zeros(12:1);
    
        for N = 1:12
            v(N) = Isp * g * N * log(1/( (pi_PL^(1/N))*(1-epsilon) + epsilon)); 
        end
    
    stairs(1:12,v)
    xlabel('N stages') 
    ylabel('Velocity (m/s)')
    xlim([1 13])
    grid on
    
%end





