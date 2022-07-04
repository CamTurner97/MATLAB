function EGA215_954313_A1_Q2
    clc

    m0 = 12000;        %lift off mass (kg)
    n = 15;            %mass ratio
    T2W = 1.4;         %thrust to weight ratio
    Isp = 350;         %specific impulse (sec)
    rho0 = 1.225;      %sea level density of atmosph 
    hscale = 7500;     %density scale height (m)
    g = 9.81;          %gravity constant
    
    W = m0 * g;        %Weight intital (N)
    T = T2W * W;       %thrust at initial weight
    mdot = T/(Isp * g);%exhaust mass flow rate with initial thrust
    mp = m0 - m0/n;    %propellent mass
    tb0 = mp/mdot;     %time to burnout
    
    Re = 6378;         %earths radius at the equator (km)
    
    
    [t , u] = ode45(@Q2a ,[0,tb0], [0; 90; 0; 0]);
    
    figure(1)
    plot(t,u(:,3))
    xlabel('Time (sec)') 
    ylabel('Altitude (m)')
    title('Altitude vs Time For Rocket')
    grid on
    
    figure(2)
    plot(t,u(:,1))
    xlabel('Time (sec)') 
    ylabel('Velocity (m/s)')
    title('Velocity vs Time For Rocket')
    grid on
    
    figure(3)
    plot(t,u(:,2))
    xlabel('Time (sec)') 
    ylabel('Gamma (deg)')
    title('Gamma vs Time For Rocket')
    grid on
    
    figure(4)
    plot(t,u(:,4))
    xlabel('Time (sec)') 
    ylabel('Downrange distance (m)')
    title('Downrange distance vs Time For Rocket')
    grid on
    
    %-----------------------------------%
    %Q2 d)
    %-----------------------------------%
    
    figure(5)
    plot(u(:,4),u(:,3)/1000)
    xlabel('Downrange distance (m)') 
    ylabel('Altitude (km)')
    title('Downrange distance vs Altitude (Rocket Trajectory)')
    axis equal
    grid on
    
    %-----------------------------------%
    %Q2 e)
    %-----------------------------------%

    %Burn out altitude (km)
    hbo = max(u(:,3))/1000;
    
    %burn out velocity 
    vbo = max(u(:,1));
    
    %burn out pitch angle 
    gammab0 = max(u(:,2));
    
    %minimal orbital velocity
    mu = 398600;  %km3/s2
    vorbit = (sqrt(mu/(Re + hbo))) * 1000;
    
    fprintf('Q2 e)\nBurnout Velocity: %f m/s \nBurnout Altitude: %f km \nMinimal orbital velocity: %f m/s \nBurnout Pitch Angle: %f deg\n',vbo,hbo,vorbit,gammab0)
    %-----------------------------------%
    
    rho = rho0*exp((-(u(:,3))/hscale)); %density variation equation 
    v = u(:,1); %velocity values
    q = 0.5.*rho.*v.^2; %dynamic pressure equation
    
    figure(6)
    plot(t,q)
    xlabel('Time (sec)') 
    ylabel('Dynamic pressure (Pa)')
    title('Dynamic pressure vs Time')
    grid on
    
    qmax = max(q(:,1));
    h = 9782.85196000901/1000;
    
    h_qmax = 9782.85196000901; %altitude (m) took from data
    v_qmax = 490.825600514135; %velocity (m/s) took from data
    rho_qmax = rho0*exp((-(h_qmax)/hscale));
    Dynamic_P = 0.5*rho_qmax*v_qmax^2; 
   
    fprintf('Q2 f)\nMax Q: %f Pa\nAltitude: %f km\n', qmax, h)
    fprintf('Numerically determined Q max: %f Pa \n', Dynamic_P)
end

    function du = Q2a(t,u)
        v = u(1);
        gamma = u(2);
        h = u(3);
        x = u(4);
        
        g = 9.81; 
        m0 = 12000;
        rho0 = 1.225;     
        hscale = 7.5; 
        Re = 6378;
        d = 1.2;        %vehicle diameter (m)
        A = pi*d^2/4;   %area of rocket (m2)
        Cd = 0.5;       %drag coefficient
        T2W = 1.4;
        W = m0 * g;         
        T = T2W * W;  
        Isp = 350;
        htower = 110;   %tower height (m)
        mdot = T/(Isp * g);
        
        rho = rho0*exp(-h/hscale);
        gh = g/((1 + h/Re)^2); %gravity variation with altitude
        
        D = A * Cd *0.5* rho * v^2; %Drag equation
        m = m0-mdot*t;
        
        if h<htower
            dv = (T-D)/m - gh*sind(gamma);
            dgamma = (v/(Re+h)-(gh/(max(v,1e-6))))*cosd(gamma);
            dh = v*sind(gamma);
            dx = Re/(Re+h)*v*cosd(gamma);
        else
            dv = (T-D)/m - gh*sind(89.75);
            dgamma = (v/(Re+h)-(gh/(max(v,1e-6))))*cosd(89.75);
            dh = v*sind(89.75);
            dx = Re/(Re+h)*v*cosd(89.75);
        end
        
        du = [dv; dgamma; dh; dx];
       
    end