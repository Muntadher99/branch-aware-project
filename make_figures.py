"""Draw Figs. 2-7 of the paper (column width, 8.2 cm) from paper_results.mat.

    python make_figures.py                          results/paper_results.mat -> output/figures/manuscript
    python make_figures.py output/paper_results.mat  -> output/figures/manuscript
    python make_figures.py <results.mat> <folder>   PNG (600 dpi) and PDF into <folder>

Requires Python 3 with numpy, scipy and matplotlib. Fig. 1 is a circuit diagram.
"""
import os
import sys

import numpy as np
import scipy.io as sio
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Liberation Serif', 'TeX Gyre Termes', 'Times New Roman', 'DejaVu Serif'],
    'mathtext.fontset': 'custom', 'mathtext.rm': 'serif', 'mathtext.it': 'serif:italic',
    'mathtext.bf': 'serif:bold', 'mathtext.sf': 'serif', 'mathtext.tt': 'serif', 'mathtext.cal': 'serif',
    'font.size': 10, 'axes.labelsize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'legend.fontsize': 9, 'axes.linewidth': 0.6, 'xtick.direction': 'out', 'ytick.direction': 'out',
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'savefig.dpi': 600, 'axes.grid': True, 'grid.color': '#dcdcd6', 'grid.linewidth': 0.5,
})
C = dict(blue='#2a78d6', orange='#eb6834', aqua='#1baf7a', violet='#4a3aa7', grey='#898781',
         light='#c3c2b7', faint='#e1e0d9', ink='#0b0b0b', ink2='#52514e')
W = 8.2 / 2.54
MA = r'$M_a$'


def load(path):
    return sio.loadmat(path, squeeze_me=True, struct_as_record=False)['R']


def thd49(TH, P):
    """THD49 of each row of TH (rad), closed form, no DC-link imbalance."""
    TH = np.sort(np.atleast_2d(TH), axis=1)
    s = np.asarray(P.s, float)
    n = np.asarray(P.ORD, float)
    b = 4 / (n[None, :] * np.pi) * (np.cos(n[None, :, None] * TH[:, None, :]) @ s)
    b1 = 4 / np.pi * (np.cos(TH) @ s)
    return np.sqrt(np.sum(b ** 2, axis=1)) / np.abs(b1)


def owner_changes(G, own):
    ok = np.isfinite(own[1:]) & np.isfinite(own[:-1])
    ch = np.where((own[1:] != own[:-1]) & ok)[0]
    return 0.5 * (G[ch] + G[ch + 1])


def save(f, folder, name):
    f.savefig(os.path.join(folder, name + '.png'), dpi=600)
    f.savefig(os.path.join(folder, name + '.pdf'))
    plt.close(f)


def fig2_envelope(R, out):
    D = R.D
    G = D.grid
    THD = 100 * D.THD
    sel = set(int(k) - 1 for k in np.atleast_1d(D.sel))
    xc = owner_changes(G, D.own)
    f, a = plt.subplots(figsize=(W, 2.15), constrained_layout=True)
    for k in range(THD.shape[1]):
        a.plot(G, THD[:, k], '-', color=C['grey'] if k in sel else C['light'], lw=0.7 if k in sel else 0.45, zorder=1)
    a.plot(G, 100 * D.env, '-', color=C['blue'], lw=1.7, zorder=3)
    a.plot(G, 100 * R.E2.cont.thd, '--', color=C['aqua'], lw=1.2, zorder=2)
    a.plot(xc, 0.25 * np.ones_like(xc), '^', color=C['ink'], ms=3.5, zorder=4, clip_on=False)
    a.set_ylim(0, 10)
    a.set_xlim(G[0], G[-1])
    a.set_ylabel(r'THD$_{49}$ (%)')
    a.set_xlabel('modulation index ' + MA)
    h = [Line2D([], [], color=C['blue'], lw=1.7), Line2D([], [], color=C['aqua'], ls='--', lw=1.2),
         Line2D([], [], color=C['grey'], lw=0.7), Line2D([], [], color=C['ink'], marker='^', ls='', ms=3.5)]
    a.legend(h, ['per-$M_a$ optimum (lower envelope)', 'continuity-constrained path',
                 '%d branches (%d stored: dark)' % (THD.shape[1], len(sel)), '%d owner changes' % len(xc)],
             loc='upper left', ncol=1, handlelength=1.8, borderaxespad=0.2, labelspacing=0.25,
             frameon=True, facecolor='w', edgecolor='none', framealpha=0.9)
    save(f, out, 'fig2_envelope')


def fig3_angles(R, out):
    D = R.D
    G = D.grid
    f, a = plt.subplots(figsize=(W, 1.8), constrained_layout=True)
    for x in owner_changes(G, D.own):
        a.axvline(x, color=C['faint'], lw=0.7, zorder=0)
    a.plot(G, D.envTH * 180 / np.pi, '.', color=C['blue'], ms=2.2)
    a.set_ylim(0, 90)
    a.set_xlim(G[0], G[-1])
    a.set_yticks([0, 30, 60, 90])
    a.grid(False, axis='x')
    a.set_ylabel('angle (deg)')
    a.set_xlabel('modulation index ' + MA)
    save(f, out, 'fig3_angles')


def fig4_floor(R, out):
    S = R.E2b.split
    Rr = S[:, 0]
    f, a = plt.subplots(figsize=(W, 2.1), constrained_layout=True)
    fl = 1e-7
    s3 = np.isfinite(S[:, 6])
    s2 = np.isfinite(S[:, 5])
    for i in np.where(s3)[0]:
        a.plot([Rr[i] / 1.04] * 2, np.maximum([S[i, 8], S[i, 10]], fl), '-', color=C['orange'], lw=0.9)
    for i in np.where(s2)[0]:
        a.plot([Rr[i] * 1.04] * 2, np.maximum([S[i, 7], S[i, 9]], fl), '-', color=C['orange'], lw=0.9)
    h1, = a.plot(Rr[s3] / 1.04, S[s3, 6], '-o', color=C['orange'], mfc=C['orange'], ms=4.5, lw=1.2)
    h2, = a.plot(Rr[s2] * 1.04, S[s2, 5], '--s', color=C['orange'], mfc='w', ms=4.5, lw=1.2)
    r0, y0 = Rr[s2][0], S[s2, 5][0]
    rr = np.array([Rr[s2][0], Rr[s2][-1]])
    a.plot(rr, y0 * (rr / r0) ** -2, ':', color=C['grey'], lw=0.9)
    a.text(rr[-1] * 0.62, y0 * (rr[-1] / r0) ** -2 * 2.2, r'slope $-2$', color=C['ink2'], ha='left', va='bottom',
           fontsize=9)
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim(8, 125)
    a.set_ylim(1e-7, 1e4)
    a.set_xticks([10, 25, 50, 100])
    a.set_xticklabels(['10', '25', '50', '100'])
    a.set_yticks([1e-6, 1e-4, 1e-2, 1, 1e2])
    a.minorticks_off()
    a.set_xlabel('table rows $R$')
    a.set_ylabel(r'THD$_{49}$ excess (pp)')
    a.legend([h1, h2], ['rows on different branches', 'rows on one branch'], loc='upper right',
             ncol=1, handlelength=2.2, borderaxespad=0.1)
    save(f, out, 'fig4_floor')


def fig5_midpoints(R, out):
    D = R.D
    T = R.E2b.test
    P = D.P
    G = D.grid
    env = 100 * T.env
    ta = 100 * thd49(T.Ya, P)
    tb = 100 * thd49(T.Yb, P)
    xc = owner_changes(G, D.own)
    nb = int(R.TIM.storage.ba_sets) if hasattr(R, 'TIM') else 99
    f, a = plt.subplots(figsize=(W, 2.2), constrained_layout=True)
    a.plot(T.Ma, env, '-', color=C['light'], lw=3.2, solid_capstyle='round', zorder=1)
    a.plot(T.Ma, ta, '-', color=C['orange'], lw=0.9, zorder=2)
    a.plot(T.Ma, tb, '-', color=C['blue'], lw=0.9, zorder=3)
    a.plot(xc, 0.115 * np.ones_like(xc), '^', color=C['ink'], ms=3.2, zorder=4, clip_on=False)
    a.set_yscale('log')
    a.set_ylim(0.1, 200)
    a.set_xlim(0.70, 1.20)
    a.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
    a.set_yticks([0.1, 1, 10, 100])
    a.set_yticklabels(['0.1', '1', '10', '100'])
    a.minorticks_off()
    a.set_xlabel('modulation index $M_a$')
    a.set_ylabel(r'THD$_{49}$ (%)')
    h = [Line2D([], [], color=C['orange'], lw=0.9), Line2D([], [], color=C['blue'], lw=0.9),
         Line2D([], [], color=C['light'], lw=3.2), Line2D([], [], color=C['ink'], marker='^', ls='', ms=3.2)]
    a.legend(h, ['plain table, %d rows' % int(T.R), 'branch-aware, %d sets' % nb, 'optimum', 'owner change'],
             loc='lower left', bbox_to_anchor=(-0.03, 0.99), ncol=2, handlelength=1.6, borderaxespad=0.1,
             columnspacing=0.9, handletextpad=0.5, labelspacing=0.2)
    save(f, out, 'fig5_midpoints')


def fig6_spectrum(R, out):
    S = R.SIM
    spec = {w.name: w for w in np.atleast_1d(S.spec)}
    hdr = list(S.header)
    thc = {row[0]: row[hdr.index('THD49_circuit_pct')] for row in S.data}
    pairs = [('lut100', 'table, 100 rows', C['orange'], -0.45),
             ('branch_at_lut100', 'branch-aware table', C['blue'], 0.45)]
    f, a = plt.subplots(figsize=(W, 2.05), constrained_layout=True)
    base = 1e-2
    hs = []
    for nm, lab, col, off in pairs:
        n = np.asarray(spec[nm].n, float)
        c = np.asarray(spec[nm].circuit, float)
        m = np.asarray(spec[nm].model, float)
        bd = (n > 1) & (n % 2 == 1) & (n % 3 != 0)
        a.bar(n[bd] + off, np.maximum(c[bd], base) - base, bottom=base, width=0.85, color=col, lw=0)
        a.plot(n[bd] + off, np.maximum(m[bd], base), 'o', color=C['ink'], mfc='w', ms=2.6, mew=0.7)
        hs.append(Patch(color=col, label='%s (%.2f %%)' % (lab, thc[nm])))
    hs.append(Line2D([], [], marker='o', color=C['ink'], mfc='w', ls='', ms=2.6, mew=0.7, label='closed form'))
    a.set_yscale('log')
    a.set_ylim(base, 1e3)
    a.set_xlim(3, 51)
    a.set_xticks([5, 11, 17, 23, 29, 35, 41, 47])
    a.set_xlabel('harmonic order $n$')
    a.set_ylabel(r'$|V_n|/|V_1|$ (%)')
    a.legend(handles=hs, loc='upper right', ncol=1, handlelength=1.2, borderaxespad=0.1, labelspacing=0.25)
    a.grid(False, axis='x')
    save(f, out, 'fig6_spectrum')


def fig7_switch(R, out):
    w = {x.name: x for x in np.atleast_1d(R.SIM.step_wave)}['step_jump']
    t = 1e3 * np.asarray(w.t, float)
    i = np.asarray(w.i, float)
    r = np.asarray(w.iref, float)
    I1 = np.max(np.abs(r))
    d = 100 * (i - r) / I1
    f, a = plt.subplots(figsize=(W, 1.75), constrained_layout=True)
    a.axvline(0, color=C['ink2'], ls=':', lw=0.8)
    a.axhline(0, color=C['ink2'], lw=0.5)
    a.plot(t, d, '-', color=C['blue'], lw=1.0)
    a.set_ylim(-3, 3)
    a.set_xlim(-20, 40)
    a.set_xlabel('time after the branch switch (ms)')
    a.set_ylabel(r'$(i_a-i_a^{new})/I_1$ (%)')
    pk = d[t >= 0][:200].max()
    a.annotate('%.2f %%' % pk, xy=(0.4, pk), xytext=(6, 1.6), fontsize=9,
               arrowprops=dict(arrowstyle='-', lw=0.6, color=C['ink2']))
    a.text(-19, 2.75, 'old pattern', fontsize=9, va='top')
    a.text(39, 2.75, 'new pattern', fontsize=9, va='top', ha='right')
    save(f, out, 'fig7_switch')


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, 'results', 'paper_results.mat')
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, 'output', 'figures', 'manuscript')
    os.makedirs(out, exist_ok=True)
    R = load(src)
    for fig in (fig2_envelope, fig3_angles, fig4_floor, fig5_midpoints, fig6_spectrum, fig7_switch):
        fig(R, out)
        print('written', os.path.join(out, fig.__name__ + '.png'))


if __name__ == '__main__':
    main()
