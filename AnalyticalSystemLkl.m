function [lkl,dlkl] = AnalyticalSystemLkl(x,data)

    % Analytical calculation of likelihood for stiffnesses given some measurement
    %  of natural frequencies of a 3-DOF-System
    % gradient calculated with autodiff
    %                                                                   /
    %                                                                   /
    %           _______            _______            _______           /
    %          |       |          |       |          |       |          /
    %          |   m1  |          |   m2  |          |   m3  |          /
    %    k1    |       |    k2    |       |    k3    |       |    k4    /
    %---vvvv---|_______|---vvvv---|_______|---vvvv---|_______|---vvvv---/
    %___________0____0_____________0____0_____________0____0____________/
    %////////////////////////////////////////////////////////////////////

    sigma = 1e-1; % assumed std of measurement

    % spring stiffnesses
    k1 = x(1);
    k2 = x(2);
    k3 = x(3);
    k4 = x(4);
    
    % masses
    m = 100;
    m1 = m;
    m2 = 2*m;
    m3 = 3*m;

    K = [k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3+k4];
    
    % some voodoo to solve third order equations
    a = -m1*m2*m3;
    b = K(1,1)*m2*m3+m1*K(2,2)*m3 + m1*m2*K(3,3);
    c = -K(1,1)*K(2,2)*m3 - K(1,1)*K(3,3)*m2 - K(2,2)*K(3,3)*m1+K(2,3)*K(3,2)*m1+K(1,2)*K(2,1)*m3;
    d = K(1,1)*K(2,2)*K(3,3)-K(1,1)*K(2,3)*K(3,2)-K(1,2)*K(2,1)*K(3,3);
    
    p = 3*a*c-b^2; q = 2*b^3-9*a*b*c+27*a^2*d;
    phi = acos(-q/(2*sqrt(-p^3)));
    y = sqrt(-p)*2*cos(phi/3+(0:2)*2*pi/(3));
    aex = (y-b)/(3*a);

    lkl = -sum((aex-data).^2/sigma.^2);

    if nargout>1
        dlkl = dlgradient(lkl,x);
    end
    
end

