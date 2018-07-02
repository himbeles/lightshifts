import numpy as np
import lightshifts.lightshift_solver as ls_solver
from lightshifts.consts import c
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import json
from sympy.physics.wigner import wigner_6j


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


def plot_total_lightshift_around_hyperfine_state(atom_filename, transitions_filename,
                                       state_f, Ff, q,
                                       Fi=None, df_min=10e9, df_max=10e9, n=200):

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
    plt.xlabel(r'detuning from $%s \rightarrow %s\,(F=%1.0f/2)$ (GHz)'%(ls.state_i[1],
                                                                        state_f[1],
                                                                        2*Ff))
    plt.ylabel(r'$\Delta V_\mathrm{ac}/I_0$ [Hz/($W/\mathrm{cm}^2$)]')
    plt.legend(loc='upper left', title='$m_F$', frameon=False)
    plt.title(title[q])


def plot_scalar_lightshift(atom_filename, transitions_filename,
              lam_min=200e-9, lam_max=1500e-9, n=200):

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


class atom():
    def __init__(self, atom_dict):
        self.atom_dict = atom_dict
        self.states = self.atom_dict['states']

        # flatten hirarchy of states dict
        self.states_flat = {}
        for s0 in self.states:
            for s1 in self.states[s0]:
                self.states_flat[(s0,s1)] = self.states[s0][s1]

    @classmethod
    def from_json(cls, atom_filename):
        atom_dict = cls._import_states_from_json(atom_filename)
        return cls(atom_dict)

    @classmethod
    def _import_states_from_json(cls, atom_filename):
        try:
            atom_dict = cls._import_json(atom_filename)
        except:
            raise Exception('Atom file not found or not in json format.')
        return atom_dict

    @staticmethod
    def _import_json(fn):
        with open(fn, 'r') as f:
                content = json.load(f)
        return content


    @staticmethod
    def _parse_term_SLJ(term):
        S, L, J = term
        S = (int(S) - 1) / 2
        L = ('S', 'P', 'D', 'F', 'G').index(L)
        J = int(J)

        return S, L, J

    @classmethod
    def _dipole_coeff(cls, state_i, state_f, parity_selection=True):
        """fine-structure reduced matrix element, Steck 7.274
        modified with additional delta function for unchanged spin
        and selection rule on wave function parity:
        delta_l (small) = +-1 (parity must change)
        """
        S, L, J = cls._parse_term_SLJ(state_i[1])
        Sf, Lf, Jf = cls._parse_term_SLJ(state_f[1])

        if parity_selection:
            if cls.parity(state_i) == cls.parity(state_f): return 0
        if S != Sf:
            return 0
        A = (-1)**(Jf + L + 1 + S)
        B = np.sqrt((2*Jf + 1) * (2*L + 1))
        C = wigner_6j(L, Lf, 1, Jf, J, S)

        return A*B*C

    @classmethod
    def parity(cls, state):
        """
        get electron configuration parity for two electron system
        used for selection rule for calculating branching ratio.
        """
        def one_electron_parity(t):
            if (t in ('s','d','g')):
                return 1
            elif (t in ('p','f')):
                return -1
            else: return '?'
        try:
            orbitals = [s for s in state[0]][1::2]
            parities = [one_electron_parity(o) for o in orbitals]
            return np.product(parities)
        except:
            return '?'

    def frequency(self, state):
        """
        Frequency of state in Hz
        """
        configuration, term = state[0], state[1]
        return self.states[configuration][term]['frequency']


    def transition_omega(self, state_i, state_f):
        """
        calculate transition angular frequency
        """
        fi = self.frequency(state_i)
        ff = self.frequency(state_f)
        return 2*np.pi * (ff-fi)

    def states_below(self, state):
        """give a dictionary of all states lower in energy than state"""
        return {k: v for k, v in self.states_flat.items() \
                if v['frequency']<self.frequency(state)}

    def states_above(self, state):
        """give a dictionary of all states larger in energy than state"""
        return {k: v for k, v in self.states_flat.items() \
                if v['frequency']>self.frequency(state)}


    def branching_ratios_LS_dict(self, state_i, parity_selection=True, verbose=False):
        """
        Estimate branching ratios between LS coupling states. Gives branching
        ration of an initial state state_i into all energetically lower lying states
        in the imported state library.
        Useful when only the lifetime of a state is known and
        not the decay rates into final states.
        Estimates the ratio of dipole matrix elements between a slection of states,
        using the parity selection rule (e.g. 6s6s -> 6s6p is allowed)
        for the electron configuration and angular momentum selection.

        Arguments:
            state_i (config,term)-tuple: initial state that branches into final states
            parity_selection (bool): use parity selection?
            verbose (bool): give more printed output
        Returns:
            dict: of (state: branching ratio) entries
        """
        # all the states energetically below state_i
        fs = list((s for s in self.states_below(state_i)))
        ratios = np.zeros(len(fs))

        for j, f in enumerate(fs):
            omega = self.transition_omega(state_i, f)
            try:
                M = self._dipole_coeff(state_i, f, parity_selection=parity_selection)
                ratios[j] = (omega**3 * M**2)
            except:
                pass

        branching_ratios = ratios/np.sum(ratios)
        if verbose:
            for j, f in enumerate(fs):
                b = branching_ratios[j]
                if b>0: print('branching ratio into ', f, '=', b)
        return dict(zip(fs, branching_ratios))

    def branching_ratio_LS(self, state_i, state_f, parity_selection=True, verbose=False):
        """
        Estimate branching ratio between LS coupling state state_i and all energetically
        lower lying states in the imported state library.
        Useful when only the lifetime of a state is known and
        not the decay rates into final states.
        Estimates the ratio of dipole matrix elements between a slection of states,
        using the parity selection rule (e.g. 6s6s -> 6s6p is allowed)
        for the electron configuration and angular momentum selection.

        Arguments:
            state_i (config,term)-tuple: initial state that branches into state_f
            state_f (config,term)-tuple: final state
            parity_selection (bool): use parity selection?
            verbose (bool): give more printed output
        Returns:
            float: branching ratio
        """

        ratios = self.branching_ratios_LS_dict(state_i, parity_selection=parity_selection,
                verbose=verbose)
        br = ratios[state_f]
        if verbose:
            print('%s -> %s: %1.4f'%(state_i, state_f, br))
        return br
