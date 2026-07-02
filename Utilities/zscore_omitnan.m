function z = zscore_omitnan(x)

x = double(x);

mu = mean(x, 'omitnan');
sd = std(x, 0, 'omitnan');

if sd == 0 || isnan(sd)
    z = x * NaN;
else
    z = (x - mu) ./ sd;
end

end