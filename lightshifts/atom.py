import numpy as np
import json
from .consts import eps0, hbar, h, c
from sympy.physics.wigner import wigner_6j


class atom():
    def __init__(self, atom_dict):
        self.atom_dict = atom_dict
        self.states = self.atom_dict['states']
        self.I = self.atom_dict['I']

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
            if len(orbitals)==0:
                return '?'
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
                ratios[j] = np.abs(omega**3 * M**2)
            except:
                pass

        if np.sum(ratios)==0:
            branching_ratios = np.array([np.nan for r in ratios])
        else:
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
