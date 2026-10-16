function d = shm_device(Tvj)
%SHM_DEVICE Datasheet data of the IGBT and diode used in the non-ideal models.
%   d = SHM_DEVICE(Tvj) returns typical values at junction temperature Tvj
%   (25 or 175 C; default 25) for the Infineon IKW20N65ET7 (650 V, 20 A),
%   datasheet Rev. 1.20. Conduction: least-squares line v = V0 + r*i through
%   the typical output and forward characteristics over 0.5-6 A. Switching
%   times are not used by the voltage models.

if nargin < 1 || isempty(Tvj), Tvj = 25; end
d.Iref = 20;   d.Vref = 400;   d.kv = 1.0;
d.E_I      = [4     5     6     8     10    12    15    20    25    30    40    50    60];
d.Eon_g    = [0.254 0.295 0.336 0.420 0.508 0.598 0.740 1.000 1.298 1.644 2.525 3.587 4.698];
d.Eoff_g   = [0.257 0.307 0.357 0.455 0.551 0.645 0.782 1.000 1.205 1.405 1.827 2.279 2.739];
d.ki       = 1.0;

switch Tvj
    case 25
        d.name     = 'Infineon IKW20N65ET7 (650 V, 20 A, IGBT7 + diode), typical, Tvj = 25 C';
        d.Vce0     = 0.700;   d.rce = 0.0467;
        d.Vf0      = 0.841;   d.rf  = 0.0655;
        d.Eon_ref  = 0.36e-3;  d.Eoff_ref = 0.36e-3;
        d.Erec_ref = 0.25 * 0.44e-6 * 400;
        d.ki_rec   = log(0.44/0.28)/log(2);
        d.Tf       = 20e-9;
    case 175
        d.name     = 'Infineon IKW20N65ET7 (650 V, 20 A, IGBT7 + diode), typical, Tvj = 175 C';
        d.Vce0     = 0.482;   d.rce = 0.0749;
        d.Vf0      = 0.549;   d.rf  = 0.0797;
        d.Eon_ref  = 0.58e-3;  d.Eoff_ref = 0.61e-3;
        d.Erec_ref = 0.25 * 1.18e-6 * 400;
        d.ki_rec   = log(1.18/0.85)/log(2);
        d.Tf       = 75e-9;
    otherwise
        error('shm_device:Tvj', 'datasheet values exist for Tvj = 25 or 175 C only');
end
d.Tt  = 2*d.Tf;
d.Tvj = Tvj;
end
