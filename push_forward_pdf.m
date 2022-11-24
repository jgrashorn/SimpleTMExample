function pdf=push_forward_pdf(tri_map,x)
    % pushforward density
%     if ~any(size(x)==tri_map.inputDim)
%         error("Check map dimensions!");
%     elseif size(x,1) ~= tri_map.inputDim
%         x = x';
%     end
    xinv = tri_map.Inverse(x,x);
    log_det_grad_x_inverse = - tri_map.LogDeterminant(xinv);
    log_pdf = log(mvnpdf(xinv'))+log_det_grad_x_inverse;
    pdf = log_pdf;
end