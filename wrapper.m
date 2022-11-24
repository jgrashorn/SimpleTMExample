function [post,dpost] = wrapper(x,data,dof,k0)

    % Handling of converting to physical space, setting desired stiffnesses
    % to update etc, setting up autodiff for gradient calculation
    % Also, I couldn't figure out how to vectorize MATLAB's dlfeval so
    % used a for-loop for now. But both can at least be parfor-loops!
    
    x0 = k0;

    post = zeros(1,size(x,2));
    dpost = zeros(length(dof),size(x,2));

    xPhys = x0.*ones(4,size(x,2));
    [xPhys_,lbound,ubound] = ConvertToPhysical(x);
    xPhys(dof,:) = xPhys_;

    if nargout==1
        parfor i=1:size(x,2) % could be parfor if desired
            post(i) = AnalyticalSystemLkl(xPhys(:,i),data)-.5*sum(x(:,i)'.^2);
        end
    else
        parfor i=1:size(x,2) % could also be parfor
            x_ = dlarray(xPhys(:,i),'BC');
            [lkl_,dlkl_] = dlfeval(@AnalyticalSystemLkl,x_,data);
            post(i) = double(extractdata(lkl_))-.5*sum(x(:,i)'.^2);
            dpost(:,i) = double(extractdata(dlkl_(dof))).*normpdf(x(:,i)')*(ubound-lbound)-x(:,i)';
        end
    end
    
end

