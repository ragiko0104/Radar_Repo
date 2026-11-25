function out = isar_fit_tip_parsimonious(t, y, cfg, opts)
%ISAR_FIT_TIP_PARSIMONIOUS Selected tip fit under the unified mode decision.
%   OUT = ISAR_FIT_TIP_PARSIMONIOUS(T,Y,CFG,OPTS) is a thin wrapper over
%   ISAR_IDENTIFY_TIP_MODE: it runs the four-hypothesis F ladder and returns
%   the selected fit structure, so callers that only need parameters do not
%   have to know about the decision.
%
%   The returned structure is the one produced by
%   isar_fit_tip_motion_mode for the chosen hypothesis, with:
%     .model / .theta / .rss / .rmsHz / .bic / .sd   as usual
%     .swingComponentHz   zero when the chosen structure has no swing
%     .parsimony          summary of the ladder (F statistics, threshold)
%
%   OPTS.ranking is forwarded to H3; the other hypotheses always use the
%   cheaper raw ranking inside isar_identify_tip_mode.
%
%   See also ISAR_IDENTIFY_TIP_MODE, ISAR_FIT_TIP_MOTION_MODE.

if nargin < 3 || isempty(cfg)
    cfg = isar_tip_motion_cfg();
end
if nargin < 4 || isempty(opts)
    opts = struct();
end
id = isar_identify_tip_mode(t, y, cfg, opts);
out = id.selected;
out.mode = id.mode;
out.parsimony = struct('rssH1', id.rssH1, 'rssH3', id.rssH3, ...
    'Fstat', id.Fstat, 'keptSwing', id.swing, ...
    'threshold', id.ladder.swingPromotionF, 'ladder', id.ladder);
end
