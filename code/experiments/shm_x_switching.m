function E = shm_x_switching(D, opts)
%SHM_X_SWITCHING Stage E7: the analytic model against the gate pattern.
%   E = SHM_X_SWITCHING(D, opts) compares the closed-form harmonics with the
%   exact Fourier integral and an FFT of the gate pattern at MA, sweeps the
%   dead time (TD_US) with the output clamped by the current direction, and
%   reports the share of the range whose optimum loses a pulse at each dead time.
%
%   Options: MA ([0.80 1.00 1.20]), TD_US ([0 1 2 4 8 16 32 64]), FS (2e6),
%            PHI (30 deg), verbose

if nargin < 2, opts = struct(); end
O.MA    = [0.80 1.00 1.20];
O.TD_US = [0 1 2 4 8 16 32 64];
O.F1    = 50;
O.FS    = 2e6;
O.PHI   = 30;
O.verbose = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

P = D.P;  grid = D.grid;

hdr = {'Ma','dead_time_us','dead_time_deg','narrowest_pulse_deg','THD49_model_pct', ...
       'THD49_edge_exact_pct','THD49_sampled_FFT_pct','model_vs_exact_rel_dev', ...
       'THD49_with_dead_time_pct','WTHD49_pct','b1','pulses_lost'};
rows = {};
for m = O.MA
    [dg, gi] = min(abs(grid - m));
    if dg > 1e-9 || ~isfinite(D.env(gi)), continue; end
    th = D.envTH(gi, :);
    for td = O.TD_US
        R = shm_switchcheck(th, P, struct('td_us', td, 'f1', O.F1, 'Fs', O.FS, ...
                                          'phi_deg', O.PHI, 'verbose', false));
        rows(end+1, :) = {grid(gi), td, 360*O.F1*td*1e-6, R.narrowest_pulse_deg, ...
            100*R.THD49_model, 100*R.THD49_exact, 100*R.THD49_fft, R.max_rel_dev, ...
            100*R.THD49_dead, 100*R.WTHD_dead, R.b1_dead, R.pulses_lost}; %#ok<AGROW>
    end
    if O.verbose
        r0 = rows(end - numel(O.TD_US) + 1, :);
        fprintf('  Ma %.2f: model %.4f%%, edge-exact %.4f%%, FFT %.4f%%, deviation %.1e\n', ...
            r0{1}, r0{5}, r0{6}, r0{7}, r0{8});
    end
end

fin = find(isfinite(D.env));
narrow = nan(numel(fin), 1);
for t = 1:numel(fin)
    x = sort(D.envTH(fin(t), :)) * 180/pi;
    edges = [x, 180 - fliplr(x)];
    w = diff([0, edges, 180]);
    w = [2 * w(1), w(2:end-1)];
    narrow(t) = min(w(w > 1e-9));
end
tdl = O.TD_US(O.TD_US > 0);
frac = zeros(size(tdl));
for i = 1:numel(tdl)
    frac(i) = mean(narrow < 360 * O.F1 * tdl(i) * 1e-6);
end

E.name   = 'E7_switching_validation';
E.header = hdr;
E.data   = rows;
E.loss   = struct('td_us', tdl, 'frac_points_losing_pulses', frac, ...
                  'narrowest_pulse_deg', narrow, 'Ma', grid(fin), 'f1', O.F1);
devs = cell2mat(rows(:, 8));
E.note = sprintf(['Switching-level check of the ideal model. The closed-form spectrum matches the ' ...
    'exact Fourier integral of the realised gate pattern to %.1e relative. Dead time is applied by ' ...
    'clamping the output to the load-current direction for t_d after each commutation (load angle ' ...
    '%g deg). A pulse narrower than 360*f1*t_d is lost outright; the loss table gives the fraction of ' ...
    'the envelope whose narrowest pulse falls below that threshold.'], max(devs), O.PHI);
end
