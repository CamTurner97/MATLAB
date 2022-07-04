function EGA215_954313_A1_Q3
    clc
    
    m0 = 100000;        %lift off mass
    T2W = 1.2;          %thrust to weight ratio
    Isp = 390;          %specific impulse
    g = 9.81;           %gravity constant
        
    n = 9;              %mass ratio
 
    Re = 6378;          %earths radius at the equator (km)
    
    W = m0 * g;         %Weight intital 
    T = T2W * W;        %thrust at initial weight
    mdot = T/(Isp * g); %echaust mass flow rate with initial thrust
    mp = m0 - m0/n;     %propellent mass
    tbo = mp/mdot;      %time to burnout

    [t , u] = ode45(@Q3a ,[0,tbo], [0.01; 90; 0; 0]);
     
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
    ylabel('Downrange Distance (m)')
    title('Downrange Distance vs Time For Rocket')
    grid on
        
    hbo1 = max(u(:,3))/1000;
    vbo1 = max(u(:,1))/1000; % (km/s)
    gammabo1 = min(u(:,2));
    
    mu = 398600; %km3/s2
    vorbit = sqrt(mu/(Re + hbo1)); %minmal orbit velocity
    
    fprintf('Gamma at burnout: %f deg \nVelocity at burnout: %f km/s \nAltitude at burnout: %f km\nMinimal orbit velocity: %f km/s \n', gammabo1, vbo1, hbo1, vorbit)
    fprintf('The mass ratio (n) value: %f \nThe Gimbal angle (alpha): -0.62 deg\n', n)

    %-----------------------------------%
    %Q3 e)
    %-----------------------------------%
    
    figure(5)
    plot(u(:,4),u(:,3))
    xlabel('Downrange distance (m)') 
    ylabel('Altitude (km)')
    title('Downrange distance vs Altitude (Rocket Trajectory)')
    axis equal
    grid on
    
    %-----------------------------------%   
end

function du = Q3a(t,u) 
    v = u(1);
    gamma = u(2);
    h = u(3);
    x = u(4);
  
    m0 = 100000;
    T2W = 1.2;
    Isp = 390;
    Cd = 0.35;
    d = 2.5;
    rho0 = 1.225;
    hscale = 7.5; %m
    g = 9.81;
    
    alpha = -0.62; %deg
    
    Re = 6378;
    
    W = m0 * g;         %Weight intital 
    T = T2W * W;        %thrust at initial weight
    mdot = T/(Isp * g); %echaust mass flow rate with initial thrust
    
    rho = rho0*exp(-h/hscale);
    gh = g/((1 + h/Re)^2);
    
    A = pi*d^2/4;
    D = A * Cd *0.5* rho * v^2;
    m = m0-mdot*t;
    
    
    dv = ((T*cosd(alpha))-D)/m - gh*sind(gamma);
    dgamma = (v/(Re+h)- gh/v)*cosd(gamma) + (T*sind(alpha))/m;
    dh = v*sind(gamma);
    dx = (Re*v*cosd(gamma))/(Re+h);
  
    du = [dv; dgamma; dh; dx];

end