function out = isar_fit_tip_motion_mode(t, y, modelName, opts)
%ISAR_FIT_TIP_MOTION_MODE Blind bounded fit of one micro-motion hypothesis.
%   OUT = ISAR_FIT_TIP_MOTION_MODE(T,Y,MODELNAME,OPTS) fits the candidate
%   structure MODELNAME ('H0'..'H3', see isar_tip_motion_bank) to the
%   extracted tip Doppler track Y sampled at slow times T (seconds).
%
%   Two ranking modes, selected by OPTS.ranking:
%
%     'raw'   (default) rank the seeds by their own residual, then polish the
%             best few distinct basins. The v3 two-group split accounts for
%             phase mismatch in the raw residual (see below).
%     'short' rank the seeds by the residual after a short Nelder-Mead run
%             each, then polish the best basins with a single ranking. It is a
%             stronger basin finder and fixes A = 6 deg, but it also lands on
%             the swing-at-fsw=fp twin far more often, so it must be used
%             together with isar_fit_tip_parsimonious.
%
%   A damped Gauss-Newton search was tried first and rejected: on the
%   composite model the first Gauss-Newton step overshoots by a factor of four
%   in squared residual (a 0.7 rad phase correction is being linearised), and
%   every damping level that rescues the step is far from the undamped one, so
%   the search terminated at the seed.
%
%   Frequencies are bounded below the 5 Hz Nyquist limit of the 10 Hz
%   slow-time sampling: above it the track cannot distinguish fp from an
%   alias, so allowing fp > 5 Hz would return arbitrary aliased solutions.
%
%   OPTS fields: cfg (isar_tip_motion_cfg), verbose, seed (warm start),
%   ranking ('raw'|'short'), seedPolishEvals (default 50).
%
%   OUT: model, description, active, theta (1x8), rss, rmsHz, k, bic, sd,
%        tref, residual, swingComponentHz, spinIdentifiable (always false).

if nargin < 4 || isempty(opts)
    opts = struct();
end
cfg = fieldOr(opts, 'cfg', isar_tip_motion_cfg());
bank = isar_tip_motion_bank();
idx = find(strcmp({bank.name}, char(modelName)), 1);
if isempty(idx)
    error('RadarProject:TipMode:UnknownModel', 'Unknown model %s.', char(modelName));
end
model = bank(idx);

t = double(t(:));
y = double(y(:));
n = numel(t);
tref = mean(t);
lb = [cfg.minFrequencyHz, cfg.minAngleDeg, cfg.minAngleDeg, cfg.minFrequencyHz, -pi, -pi, -inf, -inf];
ub = [cfg.maxFrequencyHz, cfg.maxAngleDeg, cfg.maxAngleDeg, cfg.maxFrequencyHz, pi, pi, inf, inf];
active = model.active;
if strcmp(model.name, 'H1')
    % The periodic phase has no wall at pi; canonicalise it after refinement.
    lb(5) = -inf;
    ub(5) = inf;
end
% Motion parameters that this hypothesis does not fit are pinned to their
% null value, so H0 really is static and H2 really is a pure swing. Leaving
% them at a seed value would silently fit a different physical structure.
nullIdx = setdiff(1:6, active(active <= 6));

pred = @(th) isar_tip_motion_doppler(t, forceNull(th, nullIdx), cfg, ...
    th(7) + th(8) * (t - tref));
resid = @(th) isar_fold_doppler(y - pred(th), cfg.prfHz);

% ---------------------------------------------------------------- seeding
c0 = mean(y);
p = polyfit(t - tref, y, 1);
base = [1, 5, 0, 1, 0, 0, c0, p(1)];
switch char(modelName)
    case 'H1'
        seeds = seedsPrecession(t, y, cfg, base);
    case 'H2'
        seeds = seedsSwing(base);
    case 'H3'
        seeds = seedsComposite(t, y, base);
    otherwise
        seeds = base;
end
if isfield(opts, 'seed') && ~isempty(opts.seed)
    % Warm start from a supplied solution: used by the noise study, where the
    % question is the local uncertainty rather than global basin selection.
    seeds = reshape(opts.seed, 1, []);
end
seeds = seeds(all(isfinite(seeds), 2), :);
if isempty(seeds)
    error('RadarProject:TipMode:InvalidSeed', 'The fit needs at least one finite seed.');
end
seeds = min(max(seeds, lb), ub);

% --------------------------------------------- rank seeds, refine the best
% Two selection schemes, both of which pool basins rather than trust a single
% ranking:
%
%   'raw'   the v3 scheme: rank by the seed's own residual, but split the
%           seeds into fw <= 2.4 Hz and fw > 2.4 Hz and take 12 basins from
%           each. The raw residual is dominated by phase mismatch, so a
%           single ranking lets accidental matches crowd out the true basin.
%   'short' rank by a short refinement (which measures basin quality and
%           fixes A = 6-10 deg) AND by the raw residual (which still reaches
%           basins the short run misses, e.g. A = 0.5 deg), then pool twelve
%           distinct basins from each ranking.
%
% Neither scheme splits by frequency in the 'short' case: the pooling is
% mutual coverage of two rankings that fail on different cases.
rssRaw = zeros(size(seeds, 1), 1);
for s = 1:size(seeds, 1)
    rssRaw(s) = localRss(resid, seeds(s, :));
end
[~, orderRaw] = sort(rssRaw);
slotsPerRanking = 12;
picked = zeros(0, 1);
keySeen = zeros(0, 3);
if strcmpi(fieldOr(opts, 'ranking', 'raw'), 'raw')
    shortTheta = seeds;
    useStrongSearch = false;
    isExtension = seeds(:, 4) > 2.41;
    for g = [false, true]
        nPass = 0;
        for q = 1:numel(orderRaw)
            i = orderRaw(q);
            if isExtension(i) ~= g
                continue;
            end
            key = round(seeds(i, [1, 4, 2]) * 20) / 20;
            if isempty(keySeen) || all(any(abs(keySeen - key) > 1e-9, 2))
                keySeen(end + 1, :) = key; %#ok<AGROW>
                picked(end + 1, 1) = i; %#ok<AGROW>
                nPass = nPass + 1;
            end
            if nPass >= slotsPerRanking
                break;
            end
        end
    end
else
    rssShort = zeros(size(seeds, 1), 1);
    shortTheta = zeros(size(seeds));
    for s = 1:size(seeds, 1)
        [th, rss] = nmFit(resid, seeds(s, :), lb, ub, active, ...
            fieldOr(opts, 'seedPolishEvals', 50), true);
        rssShort(s) = rss;
        shortTheta(s, :) = th;
    end
    [~, orderShort] = sort(rssShort);
    useStrongSearch = true;
    orders = {orderShort, orderRaw};
    for pass = 1:numel(orders)
        order = orders{pass};
        nPass = 0;
        for q = 1:numel(order)
            i = order(q);
            key = round(shortTheta(i, [1, 4, 2]) * 20) / 20;
            if isempty(keySeen) || all(any(abs(keySeen - key) > 1e-9, 2))
                keySeen(end + 1, :) = key; %#ok<AGROW>
                picked(end + 1, 1) = i; %#ok<AGROW>
                nPass = nPass + 1;
            end
            if nPass >= slotsPerRanking
                break;
            end
        end
    end
end
[bestRss, bestSeed] = min(rssRaw);
bestTheta = seeds(bestSeed, :);
for q = 1:numel(picked)
    [th, rss] = nmFit(resid, shortTheta(picked(q), :), lb, ub, active, 2500, useStrongSearch);
    if rss < bestRss
        bestRss = rss;
        bestTheta = th;
    end
end
% second pass from the winner removes any simplex-collapse residue
[refinedTheta, refinedRss] = nmFit(resid, bestTheta, lb, ub, active, 3000, useStrongSearch);
if refinedRss < bestRss
    bestTheta = refinedTheta;
end
bestTheta = forceNull(bestTheta, nullIdx);
if strcmp(model.name, 'H1')
    bestTheta(5) = atan2(sin(bestTheta(5)), cos(bestTheta(5)));
end

residual = resid(bestTheta);
bestRss = residual' * residual;
k = numel(active);
rmsHz = sqrt(bestRss / n);
bic = n * log(max(bestRss, 1e-30) / n) + k * log(n);
J = localJacobian(resid, bestTheta, active);
s2 = bestRss / max(n - k, 1);
sd = nan(1, 8);
if all(isfinite(J(:)))
    cov = s2 * pinv(J' * J);
    sd(active) = sqrt(max(diag(cov)', 0));
end

out = struct('model', char(modelName), 'description', model.description, ...
    'active', active, 'theta', bestTheta, 'rss', bestRss, 'rmsHz', rmsHz, ...
    'k', k, 'bic', bic, 'sd', sd, 'tref', tref, 'residual', residual, ...
    'swingComponentHz', swingComponent(bestTheta, cfg), ...
    'spinIdentifiable', false, 'seedCount', size(seeds, 1));
if fieldOr(opts, 'verbose', false)
    fprintf('  %s %-17s rms %9.4f Hz  fp %6.3f  a %7.4f  A %7.4f  fw %6.4f\n', ...
        out.model, out.description, out.rmsHz, bestTheta(1), bestTheta(2), ...
        bestTheta(3), bestTheta(4));
end
end

% ======================================================================
function seeds = seedsPrecession(t, y, cfg, base)
% One-dimensional fp scan; (amplitude, phase) solved linearly, alpha from
% the closed-form amplitude relation including the dwell attenuation.
tt = t - mean(t);
seeds = base;
bestRss = inf;
% seedScale converts the fitted Doppler amplitude back to a half angle
seedScale = 2 / cfg.lambdaM * cfg.tipLeverM / sqrt(2) * 2 * pi;
offset = (((0:cfg.dwellSubSamples-1) + 0.5) / cfg.dwellSubSamples - 0.5) * cfg.dwellSec;
for fp = cfg.minFrequencyHz:0.01:cfg.maxFrequencyHz
    M = [sin(2 * pi * fp * t), cos(2 * pi * fp * t), ones(numel(t), 1), tt];
    coef = M \ y;
    gain = seedScale * fp * mean(cos(2*pi*fp*offset));
    tone = hypot(coef(1), coef(2));
    alpha = asind(min(tone / max(gain, eps), 1));
    alpha = min(max(alpha, cfg.minAngleDeg), cfg.maxAngleDeg);
    phase = atan2(coef(2), coef(1));
    physicalTone = gain * sind(alpha) * sin(2*pi*fp*t + phase);
    nuisance = M(:,3:4) \ (y - physicalTone);
    r = isar_fold_doppler(y - physicalTone - M(:,3:4)*nuisance, cfg.prfHz);
    rss = r' * r;
    if rss < bestRss
        bestRss = rss;
        seeds = base;
        seeds(1) = fp;
        seeds(2) = alpha;
        seeds(5) = phase;
        seeds(7:8) = nuisance.';
    end
end
extra = repmat(seeds(1, :), 3, 1);
extra(:, 5) = seeds(1, 5) + [-pi / 2; 0; pi / 2];
seeds = [seeds; extra];
end

function seeds = seedsSwing(base)
seeds = zeros(0, 8);
for fw = 0.40:0.10:4.50
    for a0 = [1, 3, 5, 8]
        s = base;
        s(3) = a0;
        s(4) = fw;
        seeds(end + 1, :) = s; %#ok<AGROW>
    end
end
end

function seeds = seedsComposite(t, y, base)
% Swing-frequency grid covering the whole admissible band: the earlier
% 2.40 Hz ceiling left fw = 3 Hz with no nearby seed and the optimiser slid
% into the fw = fp degeneracy. The swing amplitude stays at a single 3 deg
% level: seeding 2/5/10 deg was tested and fixed A = 10 deg but dropped
% A = 3 deg from 1.00 to 0.40, so it was reverted.
seeds = zeros(0, 8);
for fp = 0.40:0.40:4.80
    for fw = 0.40:0.20:4.80
        for a0 = [2, 6]
            for aSw = 3
                s = base;
                s(1) = fp;
                s(2) = a0;
                s(3) = aSw;
                s(4) = fw;
                seeds(end + 1, :) = s; %#ok<AGROW>
            end
        end
    end
end
% data-driven seeds: the dominant slow-time tone is the precession carrier,
% so anchor fp on the spectral peaks as well as on the coarse grid
tt = t - mean(t);
yd = y - polyval(polyfit(tt, y, 1), tt);
nfft = 8192;
spec = abs(fft(yd .* hannLocal(numel(yd)), nfft));
spec = spec(1:nfft / 2 + 1);
f = (0:nfft / 2)' / (nfft * mean(diff(t)));
band = f > 0.4 & f < 4.8;
[~, ib] = sort(spec .* band, 'descend');
nPeak = min(3, sum(band));
fpeaks = f(ib(1:nPeak));
for r = 1:numel(fpeaks)
    for fw = [0.5, 1.0, 1.5, 2.0, 3.0, 4.5]
        for a0 = [1, 3, 5, 8]
            for aSw = 3
                s = base;
                s(1) = fpeaks(r);
                s(2) = a0;
                s(3) = aSw;
                s(4) = fw;
                seeds(end + 1, :) = s; %#ok<AGROW>
            end
        end
    end
end
end

% ======================================================================
function [th, rss] = nmFit(resid, theta, lb, ub, active, maxEval, strong)
% Bounded Nelder-Mead on the active subset, optimised in normalised
% coordinates so the initial simplex has comparable edge lengths.
if nargin < 7
    strong = false;
end
scale = [1, 5, 5, 1, 1, 1, 5, 1];
theta = reshape(min(max(theta(:), lb(:)), ub(:)), 1, []);
start = theta(active) ./ scale(active);
obj = @(v) rssAt(resid, theta, active, scale, lb, ub, v);
if strong
    % Explicit simplex. fminsearch always sizes its simplex at five percent
    % of each parameter, which is useless for a phase that starts at zero.
    step = [0.30, 0.30, 0.30, 0.30, 1.20, 1.20, 1.00, 0.30];
    v = nmSearch(obj, start, step(active), maxEval, ...
        lb(active) ./ scale(active), ub(active) ./ scale(active));
else
    opt = optimset('MaxFunEvals', maxEval, 'MaxIter', maxEval, ...
        'TolFun', 1e-14, 'TolX', 1e-10, 'Display', 'off');
    v = fminsearch(obj, start, opt);
end
th = theta;
th(active) = min(max(v .* scale(active), lb(active)), ub(active));
rss = localRss(resid, th);
end

function [x, fval] = nmSearch(fun, x0, step, maxEval, lb, ub)
% Nelder-Mead with an explicit initial simplex and box projection.
x0 = reshape(x0, 1, []);
n = numel(x0);
X = repmat(x0, n + 1, 1);
for i = 1:n
    X(i + 1, i) = X(i + 1, i) + step(i);
end
X = min(max(X, lb), ub);
F = zeros(n + 1, 1);
for i = 1:n + 1
    F(i) = fun(X(i, :));
end
neval = n + 1;
while neval < maxEval
    [F, ord] = sort(F);
    X = X(ord, :);
    spread = max(abs(X(2:end, :) - X(1, :)), [], 1);
    if all(spread < 1e-9 * (1 + abs(X(1, :))))
        break;
    end
    xbar = mean(X(1:n, :), 1);
    xr = min(max(xbar + (xbar - X(end, :)), lb), ub);
    fr = fun(xr);
    neval = neval + 1;
    if fr < F(1)
        xe = min(max(xbar + 2 * (xbar - X(end, :)), lb), ub);
        fe = fun(xe);
        neval = neval + 1;
        if fe < fr
            X(end, :) = xe;
            F(end) = fe;
        else
            X(end, :) = xr;
            F(end) = fr;
        end
    elseif fr < F(n)
        X(end, :) = xr;
        F(end) = fr;
    else
        if fr < F(end)
            xc = min(max(xbar + 0.5 * (xr - xbar), lb), ub);
        else
            xc = min(max(xbar + 0.5 * (X(end, :) - xbar), lb), ub);
        end
        fc = fun(xc);
        neval = neval + 1;
        if fc < min(fr, F(end))
            X(end, :) = xc;
            F(end) = fc;
        else
            for i = 2:n + 1
                X(i, :) = min(max(X(1, :) + 0.5 * (X(i, :) - X(1, :)), lb), ub);
                F(i) = fun(X(i, :));
            end
            neval = neval + n;
        end
    end
end
[fval, i] = min(F);
x = X(i, :);
end

function f = rssAt(resid, theta, active, scale, lb, ub, v)
th = theta;
th(active) = min(max(v .* scale(active), lb(active)), ub(active));
f = localRss(resid, th);
end

function f = localRss(resid, th)
r = resid(th);
r = r(:);
f = r' * r;
end

function J = localJacobian(resid, th, active)
base = resid(th);
base = base(:);
J = zeros(numel(base), numel(active));
for j = 1:numel(active)
    idx = active(j);
    step = 1e-6 * max(1, abs(th(idx)));
    tp = th;
    tp(idx) = tp(idx) + step;
    J(:, j) = (resid(tp) - base) / step;
end
end

function th = forceNull(th, idx)
% Pin parameters outside the active set to zero (H0 static, H2 pure swing).
if ~isempty(idx)
    th(idx) = 0;
end
end

function amp = swingComponent(theta, cfg)
% Leading-order single-sided amplitude (Hz) of the fitted swing sideband.
if theta(3) <= 0 || theta(4) <= 0
    amp = 0;
    return;
end
amp = 2 / cfg.lambdaM * cfg.tipLeverM / sqrt(2) * cosd(theta(2)) ...
    * deg2rad(theta(3)) * 2 * pi * theta(4);
end

function w = hannLocal(n)
w = 0.5 - 0.5 * cos(2 * pi * (0:n - 1)' / max(n - 1, 1));
end

function v = fieldOr(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = fallback;
end
end
