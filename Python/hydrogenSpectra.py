'''
Generate Table "tbl: hydrogen spectra" of Latex/VarianceArithmetic_Short.tex.

Each hydrogen spectral series is one data set whose j-th subset is a single
observed line, with mean alpha_j, within-subset deviation delta(alpha_j), and
weight N_j given by the NIST relative intensity.  Formulas (set mean) and
(set covariance) of the paper then split the total variance of alpha into a
within-subset part and an in-between part; the latter is what the table
exposes on real experimental data.

Per line, the observed vacuum wavelength lambda gives the transition energy
E_m - E_n = h c / lambda, which normalized by (1/n^2 - 1/m^2) mu_e c^2 / 2
gives alpha^2, so

    alpha_j        = sqrt( h c / lambda_j / ( (1/n^2 - 1/m^2) mu_e c^2 / 2 ) )
    delta(alpha_j) = alpha_j delta(lambda_j) / (2 lambda_j)

All atomic input is NIST data held as constants below; nothing is fetched at
run time.  Run with no argument to print the LaTeX tabular, or with -v to also
print the per-line breakdown.
'''

import math
import sys


# ---------------------------------------------------------------------------
# SOURCE 1: NIST Atomic Spectra Database (ASD), version 5.11, NIST Standard
#   Reference Database 78, https://physics.nist.gov/asd
#   A. Kramida, Yu. Ralchenko, J. Reader, and NIST ASD Team (2023).
#   Cited in the paper as \citep{NIST_ASD}.
#
#   Lines query : spectrum "H I", "Observed" wavelength data, wavelengths in
#                 vacuum for all wavelength bands, tab-delimited output with
#                 the columns obs_wl_vac(nm), unc_obs_wl, intens, conf_i,
#                 term_i, conf_k, line_ref.  Retrieved 2026-08-18.
#   Selection   : a line is used only when the database supplies all three of
#                 an observed vacuum wavelength, its uncertainty, and a
#                 relative intensity; Ritz-only lines are excluded.
#
#   Columns below, verbatim from that query:
#       m         upper level, from conf_k
#       lambda    obs_wl_vac(nm)   observed vacuum wavelength, nm
#       dlambda   unc_obs_wl       its uncertainty, nm
#       N         intens           relative intensity, used as the weight N_j
#       ref       line_ref         NIST bibliographic key of the measurement
# ---------------------------------------------------------------------------
NIST_ASD_LINES = {
    # n = 1, Lyman
    1: (
        # m   lambda (nm)  dlambda (nm)        N   ref
        ( 2,  121.56701,   0.00021,      840000, 'L12020'),
        ( 3,  102.5728,    0.0003,       250000, 'L3545'),
        ( 4,   97.2517,    0.0014,        83000, 'L8834'),
        ( 5,   94.9742,    0.0004,        33000, 'L3545'),
        ( 6,   93.7814,    0.0014,        19000, 'L8834'),
        ( 7,   93.0751,    0.0014,         5600, 'L8834'),
        ( 8,   92.6249,    0.0014,        11000, 'L8834'),
        ( 9,   92.3148,    0.0014,        12000, 'L8834'),
        (10,   92.0947,    0.0014,         7400, 'L8834'),
        (11,   91.9342,    0.0014,         6100, 'L8834'),
        (12,   91.8125,    0.0013,         5600, 'L8834'),
    ),
    # n = 2, Balmer
    2: (
        ( 3,  656.46,      0.003,        500000, 'L7400c29'),
        ( 4,  486.271,     0.005,        180000, 'L7439c30'),
        ( 5,  434.1692,    0.0006,        90000, 'L7436c29'),
        ( 6,  410.2892,    0.0006,        70000, 'L7436c29'),
        ( 7,  397.1198,    0.0006,        30000, 'L7436c29'),
        ( 8,  389.0166,    0.0006,        70000, 'L7436c29'),
        ( 9,  383.6485,    0.0006,        30000, 'L7436c29'),
        (10,  379.8987,    0.0006,        17000, 'L7436c29'),
        (11,  377.1704,    0.0006,         9000, 'L7436c29'),
        (12,  375.1217,    0.0006,         8000, 'L7436c29'),
        (13,  373.5431,    0.0006,         6000, 'L7436c29'),
        (14,  372.3005,    0.0006,         5000, 'L7436c29'),
        (15,  371.3034,    0.0006,         3300, 'L7436c29'),
        (16,  370.4913,    0.0006,         2800, 'L7436c29'),
        (17,  369.8209,    0.0006,         2300, 'L7436c29'),
        (18,  369.2602,    0.0006,         2300, 'L7436c29'),
        (19,  368.7881,    0.0006,         2000, 'L7436c29'),
        (20,  368.3871,    0.0006,         1700, 'L7436c29'),
        (21,  368.0418,    0.0006,         1700, 'L7436c29'),
        (22,  367.7423,    0.0006,         1400, 'L7436c29'),
        (23,  367.486,     0.003,          1300, 'L7400c29'),
        (24,  367.237,     0.01,           1200, 'L7400c29'),
        (25,  367.05,      0.003,          1300, 'L7400c29'),
        (26,  366.877,     0.003,          1000, 'L7400c29'),
        (27,  366.712,     0.003,           900, 'L7400c29'),
        (28,  366.569,     0.01,           1000, 'L7400c29'),
        (29,  366.445,     0.01,            900, 'L7400c29'),
        (30,  366.326,     0.01,           1100, 'L7400c29'),
        (31,  366.231,     0.003,           800, 'L7400c29'),
        (32,  366.136,     0.01,            700, 'L7400c29'),
        (33,  366.045,     0.003,           700, 'L7400c29'),
        (34,  365.969,     0.003,           700, 'L7400c29'),
        (35,  365.908,     0.01,            700, 'L7400c29'),
        (36,  365.829,     0.003,           700, 'L7400c29'),
        (37,  365.769,     0.003,           700, 'L7400c29'),
    ),
    # n = 3, Paschen
    3: (
        ( 5, 1282.1578,    0.0008,        32000, 'L7421'),
        ( 6, 1094.117,     0.013,         14000, 'L7421'),
        (12,  875.286,     0.003,          2200, 'L7400c29'),
        (14,  860.075,     0.003,          2300, 'L7400c29'),
        (17,  846.959,     0.003,           700, 'L7400c29'),
        (18,  844.027,     0.003,           700, 'L7400c29'),
        (19,  841.563,     0.003,           700, 'L7400c29'),
        (20,  839.471,     0.003,           700, 'L7400c29'),
        (21,  837.678,     0.003,           600, 'L7400c29'),
        (22,  836.13,      0.003,           500, 'L7400c29'),
        (23,  834.783,     0.003,           500, 'L7400c29'),
        (24,  833.607,     0.003,           370, 'L7400c29'),
        (25,  832.571,     0.003,           290, 'L7400c29'),
        (26,  831.655,     0.003,           290, 'L7400c29'),
        (27,  830.838,     0.003,           250, 'L7400c29'),
        (28,  830.111,     0.003,           210, 'L7400c29'),
        (29,  829.458,     0.003,           170, 'L7400c29'),
        (30,  828.87,      0.003,           140, 'L7400c29'),
        (31,  828.341,     0.003,           140, 'L7400c29'),
        (32,  827.859,     0.003,           140, 'L7400c29'),
        (33,  827.421,     0.003,           140, 'L7400c29'),
        (34,  827.021,     0.003,           110, 'L7400c29'),
        (35,  826.655,     0.003,           110, 'L7400c29'),
        (36,  826.316,     0.01,            110, 'L7400c29'),
        (37,  826.016,     0.01,            110, 'L7400c29'),
        (38,  825.726,     0.01,             90, 'L7400c29'),
        (39,  825.456,     0.01,             90, 'L7400c29'),
        (40,  825.226,     0.01,             90, 'L7400c29'),
    ),
    # n = 4, Brackett
    4: (
        ( 5, 4052.279,     0.008,         11000, 'L7421'),
        ( 7, 2166.1178,    0.0023,         8000, 'L7421'),
        ( 9, 1817.92,      0.004,          2800, 'L7421c31'),
        (10, 1736.6885,    0.0021,         3100, 'L7421'),
        (11, 1681.11,      0.003,          1600, 'L7421'),
        (12, 1641.136,     0.013,          1400, 'L7421c31'),
        (16, 1556.046,     0.006,           800, 'L7421c31'),
    ),
}

SERIES_NAME = {1: 'Lyman', 2: 'Balmer', 3: 'Paschen', 4: 'Brackett'}

# ---------------------------------------------------------------------------
# SOURCE 1 (continued): NIST ASD "Levels" query for spectrum "H I", column
#   Level (cm-1), retrieved 2026-08-18.  Used only for the S-state QFT shift
#   row, L(nS) = E(nS_1/2) - E(nP_1/2), the Lamb shift of the lower level.
# ---------------------------------------------------------------------------
NIST_ASD_LEVEL_S_CM1 = {          # E(nS_1/2), cm^-1
    1:      0.0000000000,
    2:  82258.9543992821,
    3:  97492.221700,
    4: 102823.8530211,
}
NIST_ASD_LEVEL_P12_CM1 = {        # E(nP_1/2), cm^-1;  1P does not exist
    2:  82258.9191133,
    3:  97492.211200,
    4: 102823.8485825,
}
NIST_ASD_LIMIT_CM1 = 109678.77174307    # H I ionization limit, cm^-1

# ---------------------------------------------------------------------------
# SOURCE 2: CODATA 2018 recommended values, https://physics.nist.gov/cuu
#   E. Tiesinga, P. J. Mohr, D. B. Newell, and B. N. Taylor,
#   Rev. Mod. Phys. 93, 025010 (2021).  Cited in the paper as
#   \citep{CODATA_2018}.  h and c are exact by the SI definition.
# ---------------------------------------------------------------------------
PLANCK_J_S      = 6.62607015e-34        # h, exact
LIGHT_SPEED_M_S = 299792458.0           # c, exact
ELECTRON_MASS   = 9.1093837015e-31      # m_e, kg
PROTON_MASS     = 1.67262192369e-27     # m_p, kg
ALPHA_CODATA    = 7.2973525693e-3       # best available alpha_0

# Theoretical 1S Lamb shift, L(1S) = 8172.874 MHz; the 1P level that the other
# series use does not exist, so the QED value is substituted for n = 1.
LAMB_1S_MHZ     = 8172.874
MHZ_PER_CM1     = 29979.2458            # c in MHz cm, exact

# mu_e is the effective electron mass of Formula (atomic hydrogen levels).
REDUCED_MASS = ELECTRON_MASS * PROTON_MASS / (ELECTRON_MASS + PROTON_MASS)
RYDBERG_J    = REDUCED_MASS * LIGHT_SPEED_M_S**2 / 2.0   # mu_e c^2 / 2


def alphaOfLine(n, m, lam, dlam):
    '''alpha and its deviation from one line, n --> m, of vacuum wavelength
       lam +/- dlam in nm.'''
    energy = PLANCK_J_S * LIGHT_SPEED_M_S / (lam * 1e-9)
    gap = 1.0 / n**2 - 1.0 / m**2
    alpha = math.sqrt(energy / (gap * RYDBERG_J))
    return alpha, alpha * dlam / (2.0 * lam)


def qftShift(n):
    '''L(nS)/|E_n|, the fraction of the lower level shifted by quantum field
       theory, from the NIST levels.'''
    if n == 1:
        shift = LAMB_1S_MHZ / MHZ_PER_CM1
    else:
        shift = NIST_ASD_LEVEL_S_CM1[n] - NIST_ASD_LEVEL_P12_CM1[n]
    return shift / (NIST_ASD_LIMIT_CM1 - NIST_ASD_LEVEL_S_CM1[n])


def analyzeSeries(n):
    '''Apply Formulas (set mean) and (set covariance) to the n-th series.'''
    lines = []
    for m, lam, dlam, weight, ref in sorted(NIST_ASD_LINES[n]):
        alpha, dalpha = alphaOfLine(n, m, lam, dlam)
        lines.append(dict(m=m, lam=lam, dlam=dlam, weight=weight, ref=ref,
                          alpha=alpha, dalpha=dalpha))
    total = sum(ln['weight'] for ln in lines)

    def weighted(func):
        return sum(ln['weight'] * func(ln) for ln in lines) / total

    # Formula (set mean): the j-th subset is a single line, so alpha_j is its
    # own mean and N_j is its relative intensity.
    alpha = weighted(lambda ln: ln['alpha'])
    # Formula (set covariance), the two halves of the total variance.
    within = weighted(lambda ln: ln['dalpha']**2)
    between = weighted(lambda ln: (ln['alpha'] - alpha)**2)
    variance = within + between

    return dict(
        n=n, name=SERIES_NAME[n], lines=lines, count=len(lines),
        maxLevel=max(ln['m'] for ln in lines),
        maxIntensity=max(ln['weight'] for ln in lines),
        minIntensity=min(ln['weight'] for ln in lines),
        # Intensity-weighted root mean square of the input uncertainty, both
        # absolute and relative, as the input counterpart of the precision.
        rmsDLambda=math.sqrt(weighted(lambda ln: ln['dlam']**2)),
        rmsRelDLambda=math.sqrt(weighted(lambda ln: (ln['dlam'] / ln['lam'])**2)),
        qft=qftShift(n),
        alpha=alpha, deviation=math.sqrt(variance),
        precision=math.sqrt(variance) / alpha,
        relError=(alpha - ALPHA_CODATA) / ALPHA_CODATA,
        betweenFraction=between / variance,
    )


def kilo(value):
    '''Format an intensity in units of 10^3 with 2 significant digits.'''
    value /= 1e3
    if value >= 10:
        return '%.0f' % value
    return ('%.1f' if value >= 1 else '%.2f') % value


def latexTabular(result):
    def row(label, func):
        return label + ' & ' + ' & '.join(func(r) for r in result) + r' \\'

    out = [r'\begin{tabular}{|l|' + 'r|' * len(result) + '}', r'\hline']
    out.append('series & ' + ' & '.join(r['name'] for r in result) + r' \\')
    out.append(r'\hline')
    out.append(row(r'lower level $n$', lambda r: '%d' % r['n']))
    out.append(row(r'lines used', lambda r: '%d' % r['count']))
    out.append(row(r'maximal level $m_{\max}$', lambda r: '%d' % r['maxLevel']))
    out.append(row(r'strongest intensity ($10^3$)', lambda r: kilo(r['maxIntensity'])))
    out.append(row(r'weakest intensity ($10^3$)', lambda r: kilo(r['minIntensity'])))
    out.append(row(r'intensity weighted $\delta\lambda$ ($10^{-3}$ nm)',
                   lambda r: '%.2f' % (1e3 * r['rmsDLambda'])))
    out.append(row(r'intensity weighted $\delta\lambda/\lambda$ ($10^{-6}$)',
                   lambda r: '%.2f' % (1e6 * r['rmsRelDLambda'])))
    out.append(row(r'$S$-state QFT shift ($10^{-6}$)',
                   lambda r: '%.2f' % (1e6 * r['qft'])))
    out.append(r'\hline')
    out.append(row(r'precision ($10^{-6}$)', lambda r: '%.2f' % (1e6 * r['precision'])))
    out.append(row(r'relative error ($10^{-6}$)', lambda r: '%+.2f' % (1e6 * r['relError'])))
    out.append(row(r'in-between variance percentage',
                   lambda r: r'$%.0f\%%$' % (100 * r['betweenFraction'])))
    out.append(r'\hline')
    out.append(r'\end{tabular}')
    return '\n'.join(out)


def dumpLines(result):
    for r in result:
        print('\n%s series, n=%d, %d lines, alpha = %.9f +/- %.3e'
              % (r['name'], r['n'], r['count'], r['alpha'], r['deviation']))
        print('    %3s %14s %12s %10s %14s %12s %10s'
              % ('m', 'lambda (nm)', 'dlambda', 'N_j', 'alpha_j', 'dalpha_j', 'ref'))
        for ln in r['lines']:
            print('    %3d %14.4f %12.4f %10d %14.9f %12.3e %10s'
                  % (ln['m'], ln['lam'], ln['dlam'], ln['weight'],
                     ln['alpha'], ln['dalpha'], ln['ref']))


def main(argv):
    result = [analyzeSeries(n) for n in sorted(NIST_ASD_LINES)]
    if '-v' in argv or '--verbose' in argv:
        dumpLines(result)
        print('')
    print(latexTabular(result))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
