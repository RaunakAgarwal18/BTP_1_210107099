%   Dr.Hernan Peraza    hperaza@ipn.mx
%****************************************************
function [ o ] =  alpha_melanophore(fit, min, max )
    o = zeros(size(fit));
    for i=1:size(fit,2)
         o(i)= (max-fit(i))/(max-min);
    end
end
%****************************************************
