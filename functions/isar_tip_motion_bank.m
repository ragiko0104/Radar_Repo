function bank = isar_tip_motion_bank()
%ISAR_TIP_MOTION_BANK Candidate micro-motion structures for the tip track.
%   Each hypothesis keeps the same motion primitive and differs only in which
%   terms are active:
%
%       H0 static            no micro motion
%       H1 precession        fp, alpha
%       H2 swing             fsw, A
%       H3 precession+swing  fp, alpha, fsw, A
%
%   Every hypothesis also carries the two nuisance terms c0 + c1*(t - tref).
%   They are not micro-motion parameters: c0 absorbs the unobservable mean of
%   the LOS projection and the systematic offset of the gated peak, and c1
%   absorbs the residual RVP/range-coupling drift that the CA prior leaves
%   behind (about 0.14 Hz/s in the 14-scene set).
%
%   The spin rate fs is deliberately absent. The tip lies exactly on the spin
%   axis, so Ry(2*pi*fs*t) fixes it and no value of fs can change the track;
%   fs is unidentifiable from this observable, not merely hard to estimate.
%
%   Parameter order, shared by every hypothesis:
%       1 fp [Hz]   2 alpha [deg]   3 A [deg]   4 fsw [Hz]
%       5 psiP [rad]  6 psiW [rad]  7 c0 [Hz]  8 c1 [Hz/s]

bank = struct( ...
    'name', {'H0', 'H1', 'H2', 'H3'}, ...
    'description', {'static', 'precession', 'swing', 'precession+swing'}, ...
    'active', {[7 8], [1 2 5 7 8], [3 4 6 7 8], 1:8});
end
