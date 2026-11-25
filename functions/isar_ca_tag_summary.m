function tbl = isar_ca_tag_summary(tbl, info)
% The CSV remains interpretable even when separated from the run's MAT file.
names = {'runId','recordMode','rangeModel','velocityModel','accelerationSource', ...
    'configuredRangeSigmaM','configuredVelocitySigmaMps','appliedRangeSigmaM', ...
    'appliedVelocitySigmaMps'};
for i = 1:numel(names)
    key = names{i};
    tbl.(key) = repmat(info.(key),height(tbl),1);
end
end
