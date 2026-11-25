%TEST_RVCAL_FPGA_ALIGNMENT Verify lib RV peak detection and zero-based delay convention.

clear;
thisFile = mfilename('fullpath');
repoRoot = char(java.io.File(fullfile(fileparts(thisFile),'..')).getCanonicalPath());
run(fullfile(repoRoot,'startup_project.m'));

nSample = 1024;
sampleIndex = 1:nSample;
power = ones(1,nSample);
expectedPeaks = [250,520,780];
for peakIndex = 1:numel(expectedPeaks)
    power = power + (120-10*peakIndex).* ...
        exp(-0.5.*((sampleIndex-expectedPeaks(peakIndex))./3).^2);
end
trainingCells = 64;
guardCells = 32;
alpha = radar.cfarth_alpha(trainingCells,1e-4,'P');
[candidatePeaks,~] = radar.tm_cfar(power,trainingCells,guardCells,0, ...
    alpha,'MEAN',0,round(2*trainingCells*0.8),'NO',0);

for topKCount = [16,3,2]
    [packagePeaks,~,packageThreshold] = radar.find_Peak_exact(power, ...
        candidatePeaks,30,'YES',topKCount,0.5);
    assert(numel(packagePeaks{1}) == numel(expectedPeaks));
    peakPower = sort(power(round(packagePeaks{1})),'descend');
    actualCount = min(numel(peakPower),topKCount);
    expectedThreshold = sum(peakPower(1:actualCount))/actualCount*0.5;
    assert(abs(packageThreshold-expectedThreshold) < 1e-12);
end

[emptyPackagePeaks,~,emptyPackageThreshold] = radar.find_Peak_exact( ...
    ones(1,nSample),{[]},30,'YES',16,0.5);
assert(isempty(emptyPackagePeaks{1}));
assert(emptyPackageThreshold == 0 && isfinite(emptyPackageThreshold));

fs = 200e6;
fc1 = 9.9e9;
fc2 = 10.0e9;
k1 = 2.0e12;
k2 = 2.1e12;
oldDelay1 = 100/fs;
oldDelay2 = 110/fs;
newDelay1 = (100-1)/fs;
newDelay2 = (110-1)/fs;
[oldRange,oldVelocity] = radar.Cal_rv_mat(oldDelay1,oldDelay2,fc1,fc2,k1,k2);
[newRange,newVelocity] = radar.Cal_rv_mat(newDelay1,newDelay2,fc1,fc2,k1,k2);
assert(abs((newRange-oldRange)+0.75) < 1e-9);
assert(abs(newVelocity-oldVelocity) < 1e-9);

rvText = fileread(fullfile(repoRoot,'functions','rvpath_solve_rv.m'));
assert(contains(rvText,"'peakIndexToDelayConvention'") && ...
    contains(rvText,"'(peakIndex-1)/fs'"));
for bandIndex = 1:3
    expectedFormula = sprintf("rvDelay%d = (double(rvPeak%d(:).') - 1) ./ fs;", ...
        bandIndex,bandIndex);
    assert(contains(rvText,expectedFormula));
end

fprintf('RVCal FPGA alignment test passed.\n');
