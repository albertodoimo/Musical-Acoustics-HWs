function [out] = parallel(Z1, Z2)
out = (1./Z1+1./Z2).^(-1);
end