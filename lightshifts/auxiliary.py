import numpy as np
import lightshifts.atom as atom
from lightshifts.consts import c
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm


def smart_gen_array(func, a, b, sing=[], n=100, eps=1e-17):
    """
    Evaluate function func between a and b and replace singularities
    with np.nan. Every subinterval between the singularities is filled
    with n points and we evaluate up to eps around the singularity.
    """
    # make intervals between singularities
    intervals = []
    sorted_sing = sorted(sing)

    if len(sorted_sing) == 0:
        intervals.append([a, b])
    else:
        for i, s in enumerate(sorted_sing):
            if s < a:
                continue

            # first interval
            if s > a and len(intervals) == 0:
                intervals.append([a, s])

            # last interval
            elif s > b:
                intervals.append([sorted_sing[i-1], b])

            # intermediate intervals
            else:
                intervals.append([sorted_sing[i-1], s])

    if (intervals[-1][1] < b):
        intervals.append([intervals[-1][1], b])

    # solve function within the intervals and stitch with np.nan
    # on the singularities

    xx = []
    yy = []

    for interval in intervals:
        a, b = interval

        # safety margin
        a += eps
        b -= eps

        x = np.linspace(a, b, n)
        y = func(x)

        x = np.append(x, b)
        y = np.append(y, np.nan)

        xx = np.append(xx, x)
        yy = np.append(yy, y)

    return np.array(xx), np.array(yy)


def laser_intensity(laser_power, beam_waist):
    """
    Laser intensity in Watts/cm^2 for laser power in W and and beam waist in m
    """
    return 2/np.pi * laser_power/(beam_waist/1e-2)**2

def _tex_state(state, F=None, show_term=True):
    term, config = state
    try:
        par = atom.parity(state)
        if par=='?':
            return '%s %s'%(term,config)
    except:
        pass
    if F is None:
        Fstr = ""
    else:
        if (F*2)%2==0:
            Fstr = " (F=%d)"%F
        else:
            Fstr = " (F=%d/2)"%(F*2)

    if show_term:
        term_str = '(%s)'%term
    else:
        term_str = ''
    return term_str+r"${}^{%s}\mathrm{%s}_{%s}$"%(config[0], config[1],
            config[2])+Fstr



def plot_total_lightshift_around_hyperfine_state(atom_filename, transitions_filename,
                                       state_f, Ff, q,
                                       Fi=None, df_min=10e9, df_max=10e9, n=200):

    import lightshifts.lightshift_solver as ls_solver

    # make the plot data
    ls = ls_solver(atom_filename, transitions_filename, Fi=Fi)
    fzero = ls.transition_frequency_hyperfine(state_f=state_f, Ff=Ff)
    lam_0 = c/(fzero+df_max)
    lam_1 = c/(fzero-df_min)
    sing = [x['wavelength'] for x in ls.sorted_transitions_dict(hyperfine=True)]
    if q!=0:
        possible_mf = np.arange(-ls.Fi, ls.Fi+1, 1)
    else:
        possible_mf = np.arange(ls.Fi, 0, -1)
    ls_data = []
    for mf in possible_mf:
        ls_data.append(
            smart_gen_array(lambda w: ls.total_lightshift(w, q, mf),
                            lam_0, lam_1, sing, n=n))

    # plot title
    title = {
        -1: '$\sigma^-$ light',
        1: '$\sigma^+$ light',
        0: '$\pi$ light'
    }
    # legend title
    leg_title = {
        -1: '$m_F$',
        1: '$m_F$',
        0: '$|m_F|$'
    }

    # create color map to identify mf component
    fmax = ls.Fi
    norm = mpl.colors.Normalize(vmin=-fmax, vmax=fmax)
    cmap = cm.copper
    m = cm.ScalarMappable(norm=norm, cmap=cmap)

    for i, mf in enumerate(possible_mf):
        w, l = ls_data[i]
        plt.plot((c/w-fzero)/1e9, l, c=m.to_rgba(mf), label=r'{:+1.0f}/2'.format(mf*2))

    for w in sing:
        plt.plot([(c/w-fzero)/1e9, (c/w-fzero)/1e9], [-1e6, 1e6], ':', color='gray', label='')

    plt.gca().axhline(y=0, color='k', lw=0.5)
    plt.xlim(-df_min/1e9,df_max/1e9)
    plt.xlabel(r'detuning from %s $\rightarrow$ %s (GHz)'%(_tex_state(ls.state_i),
                                                           _tex_state(state_f,
                                                               Ff)))

    plt.ylabel(r'$\Delta V_\mathrm{ac}/I_0$ [Hz/($W/\mathrm{cm}^2$)]')
    plt.legend(loc='upper left', title=leg_title[q], frameon=False)
    plt.title(title[q])


def plot_scalar_lightshift(atom_filename, transitions_filename,
              lam_min=200e-9, lam_max=1500e-9, n=200):

    import lightshifts.lightshift_solver as ls_solver

    # make the plot data
    ls = ls_solver(atom_filename, transitions_filename)
    lam_0 = lam_min
    lam_1 = lam_max
    sing = [x['wavelength'] for x in ls.sorted_transitions_dict()]
    w, l = smart_gen_array(lambda w: ls.scalar_lightshift(w), lam_0, lam_1, sing, n=n)

    plt.plot(w*1e9, l, label='(%s)%s'%(ls.state_i[0],ls.state_i[1]))

    for w in sing:
        plt.plot([w*1e9, w*1e9], [-1e6, 1e6], ':', color='gray', label='')
    plt.gca().axhline(y=0, color='k', lw=0.5)
    plt.xlim(lam_0*1e9,lam_1*1e9)
    plt.xlabel(r'wavelength (nm)')
    plt.ylabel(r'$\Delta V_\mathrm{ac}/I_0$ [Hz/($W/\mathrm{cm}^2$)]')
    plt.legend()
    plt.title('scalar light shift')


def frequency_from_wavenumber(wn):
    """
    frequency in Hz from wavenumber in 1/cm^2
    """
    return c*wn*100
