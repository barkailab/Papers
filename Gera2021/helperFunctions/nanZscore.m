function out = nanZscore(vec)
    nanPos = isnan(vec);
    out = nan(size(vec));
    out(~nanPos) = zscore(vec(~nanPos));
end