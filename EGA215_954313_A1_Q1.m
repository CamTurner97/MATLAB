function EGA215_954313_A1_Q1
    clc 
% Question 1 part (a)
    m01 = 249.5; %kg
    m02 = 113.4; %kg

    mf1 = 170.1; %kg
    mf2 = 58.97; %kg

    me1 = 11.1; %kg/s
    me2 = 4.05; %kg/s

    Isp1 = 235; %s
    Isp2 = 310; %s

    g = 9.81; %m/s2

    c1 = Isp1 * g; %effective exhaust velocity
    c2 = Isp2 * g; %effective exhaust velocity

    tb01 = (m01-mf1)/me1 ; %time to burnout of the first stage 
    vb01 = c1*log(m01/mf1)- g * tb01; %velocity after burnout
    hb01 = (c1/me1) * (mf1*log(mf1/m01) + m01 - mf1) - 0.5* g *((m01-mf1)/me1)^2; %altitude after burnout
    
    tc1 = 4; %coast time sec
    tot_time_1 = tb01 + tc1; %time taken for first burnout and coast
    hc1 = hb01 + vb01 * tc1 - (g * (tc1^2))/2;  %altitude after coasting
    vc1 = vb01 - g * tc1; %velocity after coasting 
    
    tb2dash = (m02-mf2)/me2;
    tb02 = tb01 +tc1 + tb2dash;  %time to burnout of the second stage  
    vb02 = vc1 + c2*log(m02/mf2) - g*tb2dash; %velocity after burnout at second stage
    
    hb02 = hc1 + c2*((mf2/me2)*log(m02/mf2) + tb02) - (g/2)* tb2dash^2; %altitude after burnout of seconf stage
    tc2dash = (vb02/g); %time it takes for momentum to be zero after second stage burnout
    tc2 = tb01 + tc1 + tb2dash + tc2dash; 
    hc2 = hb02 + vb02 * tc2dash - (g/2)* (tc2dash)^2; %altitude after second stage burnout
    
    fprintf("Q1 a)\nAltitude of first stage rocket before coasting: %f m, \nafter coasting it is: %f m \n" , hb01, hc1)
    fprintf("Altitude of second stage rocket before coasting: %f m, \nafter coasting it is: %f m \n" , hb02, hc2)
    
% Part (b)

    x = [0,tb01,tot_time_1,tb02,tc2];
    y = [0,hb01,hc1,hb02,hc2];
    hold on
    grid on
    xlabel('Time (sec)') 
    ylabel('Altitude (m)')
    title('Altitude vs Time For Rocket')
    plot(x,y, '-b', tb01,hb01,'r^', tot_time_1,hc1,'r^',tb02,hb02,'r^',tc2,hc2,'r^')
    
% Part (c)

    rho0 = 1.225;  %sea level atmospheric density
    Cd = 0.8;      %coefficient of drag
    para = 5;      %parachute radius (m)
    S = pi*para^2; %parachute area
    M = (m01-mf1)+(m02-mf2); %mass of rocket 
    W = M * g;     %weight of rocket
    v = sqrt(W/(Cd*0.5*rho0*S)); %velocity from rearranging drag equation (m/s)
    time = hc2/v;  %time taken for rocket to fall
    
    fprintf('Q1 c) \nAltitude: %f m\nVelocity Falling: %f m/s\nTime taken to return to earth: %f sec\n', hc2, v, time)
    
end

