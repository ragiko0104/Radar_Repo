function out = isar_identify_tip_mode(t, y, cfg, opts)
%ISAR_IDENTIFY_TIP_MODE Single decision path over the four tip hypotheses.
%   OUT = ISAR_IDENTIFY_TIP_MODE(T,Y,CFG,OPTS) fits H0..H3 and decides the
%   structure with one nested F ladder, replacing the earlier mix of a BIC
%   ranking plus a fixed 0.2 deg amplitude threshold on one side and a
%   separate H1/H3 test on the other.
%
%   The hypotheses are nested as
%
%       H0  static            c0, c1                        k = 2
%       H1  precession        fp, alpha, psiP               k = 5
%       H2  swing             A, fsw, psiW                  k = 5
%       H3  precession+swing  all six                       k = 8
%
%   The four empirical F statistics each compare three added parameters:
%
%       precession present?   F1 = ((RSS0 - RSS1)/3) / (RSS1/(n-5))
%       swing alone?          F2 = ((RSS0 - RSS2)/3) / (RSS2/(n-5))
%       swing given prec.?    F3 = ((RSS1 - RSS3)/3) / (RSS3/(n-8))
%       prec. given swing?    F4 = ((RSS2 - RSS3)/3) / (RSS3/(n-8))
%
%   By default each comparison uses CFG.parsimonyF. CFG.swingPromotionF
%   can override the F3 cutoff for promoting H1 to H3. The structure is
%   then H3 if precession and swing, H1 if precession only, H2 if the
%   swing-only model beats the static one and precession is absent, H0
%   otherwise. With n >= 18, F3 above threshold and F4 at or below
%   threshold select H2 first. The full ladder requires n > 8.
%
%   Why F and not BIC or a fixed amplitude threshold: the tracks are close to
%   deterministic, so BIC over-selects on sub-0.1 Hz extraction artefacts, and
%   a swing placed at fsw = fp reproduces the precession tone exactly
%   (A = 57.6*sin(alpha)) with an equal residual. That twin is a property of
%   the model; the ladder scores the residual improvement from extra parameters.
%
%   OPTS: ranking for H3 ('short' default here, 'raw' otherwise), verbose.
%
%   OUT fields: model, mode, swing, precession, selected (the chosen fit
%   structure from isar_fit_tip_motion_mode), fits (cell H0..H3), ladder
%   (all residuals and F values), rssH1, rssH3.

if nargin < 3 || isempty(cfg)
    cfg = isar_tip_motion_cfg();
end
if nargin < 4 || isempty(opts)
    opts = struct();
end
verbose = fieldOr(opts, 'verbose', false);
thr = cfg.parsimonyF;
swingPromotionF = fieldOr(cfg, 'swingPromotionF', thr);
validateattributes(swingPromotionF, {'numeric'}, {'scalar','real','finite','positive'}, ...
    mfilename, 'swingPromotionF');
n = numel(y);
if n <= 8
    error('RadarProject:TipMode:InsufficientSamples', ...
        'The four-model ladder needs at least 9 samples; received %d.', n);
end
names = {'H0', 'H1', 'H2', 'H3'};
fits = cell(1, 4);
for m = 1:4
    % H3 needs the stronger search: the swing basins are narrow and the raw
    % ranking misses them. The other three are low-dimensional enough that the
    % raw ranking is reliable and much cheaper.
    if m == 4
        rank = fieldOr(opts, 'ranking', 'short');
    else
        rank = 'raw';
    end
    fits{m} = isar_fit_tip_motion_mode(t, y, names{m}, ...
        struct('cfg', cfg, 'ranking', rank));
end
rss0 = fits{1}.rss;
rss1 = fits{2}.rss;
rss2 = fits{3}.rss;
rss3 = fits{4}.rss;
F1 = ((rss0 - rss1) / 3) / max(rss1 / (n - 5), realmin);
F2 = ((rss0 - rss2) / 3) / max(rss2 / (n - 5), realmin);
F3 = ((rss1 - rss3) / 3) / max(rss3 / (n - 8), realmin);
% Reverse direction of F3: does the precession still earn its three
% parameters once the swing is already in the model? Needed because a swing
% can mimic a precession tone (its Doppler is not a pure sinusoid, so H3 fits
% better than H1 even when there is no precession), which makes F1 fire on a
% pure-swing scene.
F4 = ((rss2 - rss3) / 3) / max(rss3 / (n - 8), realmin);
swingAlone = F2 > thr;
precessionGivenSwing = F4 > thr;
% Downgrade gate. Measured behaviour: a pure-swing scene gives F1 ~ 8e3 (>>
% threshold) because H1 places its single tone at fsw and soaks up the
% swing's fundamental, then F3 also passes because H3 fits the swing far
% better than H1 does; the old rule therefore reported "precession+swing"
% with alpha = 0.004 deg. If the swing alone suffices (F3 passes) and adding
% the precession on top of it does not help (F4 fails, measured 1.6 there),
% the precession term is a surrogate and the swing-only structure is the
% selected answer for that example.
% The gate needs usable residual degrees of freedom. With a 1 s aperture the
% track has only 10 samples, so only two residual degrees of freedom remain.
% The observed score fell below the threshold
% on a genuine 1 deg precession and the gate wrongly reported swing-only. The
% gate is therefore disabled when n - 8 < 10.
gateReliable = (n - 8) >= 10;
swingAloneSuffices = gateReliable && ~precessionGivenSwing && F3 > thr;
if swingAloneSuffices
    chosen = 3;
elseif F1 > thr
    if F3 > swingPromotionF
        chosen = 4;
    else
        chosen = 2;
    end
elseif swingAlone
    chosen = 3;
else
    chosen = 1;
end
out = struct();
out.model = string(names{chosen});
out.selected = fits{chosen};
out.mode = string(out.selected.description);
out.swing = chosen == 3 || chosen == 4;
out.precession = chosen == 2 || chosen == 4;
out.fits = fits;
out.ladder = struct('rssH0', rss0, 'rssH1', rss1, 'rssH2', rss2, ...
    'rssH3', rss3, 'Fprecession', F1, 'FswingAlone', F2, ...
    'FswingGivenPrecession', F3, 'FprecessionGivenSwing', F4, ...
    'threshold', thr, 'swingPromotionF', swingPromotionF, 'downgradeEnabled', gateReliable);
out.rssH1 = rss1;
out.rssH3 = rss3;
out.Fstat = F3;
if verbose
    fprintf(['  ladder: F_prec %.1f  F_swingOnly %.1f  F_swingOnPrec %.1f' ...
        '  F_precOnSwing %.1f  -> %s\n'], F1, F2, F3, F4, out.mode);
end
end

function v = fieldOr(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = fallback;
end
end
