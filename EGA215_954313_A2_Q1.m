function EGA215_954313_A2_Q1 
    clc
    format long
    
    mpl = 6000; %payload kg
    orbital_sp = 8.200; % speed km/s
    
    Isp1 = 290; %specific impulse sec
    Isp2 = 320;
    Isp3 = 350;
    
    g = 9.81/1000; %gravity
    
    epsilon1 = 0.12; %structural ratio 
    epsilon2 = 0.16;
    epsilon3 = 0.2;
    
    c1 = Isp1 * g; % m/s 
    c2 = Isp2 * g;
    c3 = Isp3 * g;
    
    f = @(x) c1.*log(((x.*c1)-1)./(x.*epsilon1.*c1)) + c2.*log(((x.*c2)-1)./(x.*epsilon2.*c2))+ c3.*log(((x.*c3)-1)./(x.*epsilon3.*c3)) - orbital_sp; 

    x_lower = 0.4;
    x_upper = 1;

    x_mid = (x_lower + x_upper)/2;

    while abs(f(x_mid)) >0.00001
        if (f(x_mid)*f(x_upper)) < 0 
            x_lower = x_mid;
        else
            x_upper = x_mid;
        end
        x_mid = (x_lower + x_upper)/2;
    
    end
    fprintf('The root is: %g\n',x_mid)
    
    lag_multi = x_mid;
    fprintf('lagrange miltiplier is: %f\n\n',lag_multi)
    
    n1 = (lag_multi*c1 - 1)/(lag_multi*epsilon1*c1);
    n2 = (lag_multi*c2 - 1)/(lag_multi*epsilon2*c2);
    n3 = (lag_multi*c3 - 1)/(lag_multi*epsilon3*c3);
    
    pi_pl_1 =(1-epsilon1*n1)/((1-epsilon1)*n1); %Payload fraction 
    pi_pl_2 =(1-epsilon2*n2)/((1-epsilon2)*n2);
    pi_pl_3 =(1-epsilon3*n3)/((1-epsilon3)*n3);
    
    m03 = mpl/pi_pl_3; 
    m02 = m03/pi_pl_2;
    m01 = m02/pi_pl_1;
    
    ratio_1 = m02/m01;
    ratio_2 = m03/m02;
    ratio_3 = mpl/m03;
    fprintf('Mass ratios are: \nStage 1-2: %f\nStage 2-3: %f\nStage 3-payload: %f\n\n',ratio_1,ratio_2,ratio_3)
    
    m1 = m01 - m02; %step masses
    m2 = m02 - m03;
    m3 = m03 - mpl;
    fprintf('The step masses are; \nm1 = %f kg\nm2 = %f kg\nm3 = %f kg\n\n',m1,m2,m3)
    
    total_mass = m1+m2+m3+mpl;
    fprintf('Tatal mass of the rocket: %f kg\n\n',total_mass)
    
    me1 = epsilon1*m1; %structural masses
    me2 = epsilon2*m2;
    me3 = epsilon3*m3;
    fprintf('Empty masses are; \nme1 = %f kg\nme2 = %f kg\nme3 = %f kg\n\n',me1,me2,me3)
    
    mp1 = m1-me1; %propellant masses
    mp2 = m2-me2;
    mp3 = m3-me3;
    fprintf('Propellant masses are; \nmp1 = %f kg\nmp2 = %f kg\nmp3 = %f kg\n\n',mp1,mp2,mp3)
    
    psi = ((epsilon1^2)/((1-epsilon1*n1)^2) - 1/(n1^2))*((c2/n2)^2) +...
          ((epsilon2^2)/((1-epsilon2*n2)^2) - 1/(n2^2))*((c3/n3)^2) +...
          ((epsilon3^2)/((1-epsilon3*n3)^2) - 1/(n3^2))*((c1/n1)^2);
    fprintf('psi = %f \nTherefore a minima.\n', psi)
end

















