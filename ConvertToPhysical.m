function [xPhys,lb,ub] = ConvertToPhysical(x)
    % transform from std normal to physical space with inverse cdf
    lb = 100;
    ub = 3000;
    xPhys = lb + normcdf(x).*(ub-lb);
end

