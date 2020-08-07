function [dff] = calc_dff_mw(F,wsz)
% calculate dF/F with a moving window

dff = zeros(size(F));
for n = 1:length(dff)
    baseline = mean(F(max([1,n-wsz]):n));
    dff(n) = F(n)/baseline;
end

% baseline = quantile(F,0.3,2);
% dff = F./(baseline*ones(1,size(F,2)));

end