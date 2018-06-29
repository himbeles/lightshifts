import numpy as np
import json
from .consts import eps0, hbar, h, c
from sympy.physics.wigner import wigner_6j


class lightshift_solver():
    def __init__(self, atom_filename, transitions_filename):
        self.atom, self.transitions = self.import_states_and_transitions(atom_filename, transitions_filename)
        self.states = self.atom['states']
        self.I = self.atom['I']
        self.Fi = self.I
        self.state_i = self.transitions[0]['state_i']

    def import_states_and_transitions(self, atom_filename, transitions_filename):
        try:
            atom = self._import_json(atom_filename)
        except:
            raise Exception('Atom file not found or not in json format.')
        try:
            transitions = self._import_json(transitions_filename)
        except:
            raise Exception('Transitions file not found or not in json format.')

        return atom, transitions

    @staticmethod
    def _import_json(fn):
        with open(fn, 'r') as f:
                content = json.load(f)
        return content

    def _calc_Fs(self, J):
        a = np.abs(J-self.I)
        b = J + self.I
        r = int(b-a)+1

        return [a + f for f in range(r)]

    def _query_state_prop(self, state, prop):
        configuration, term = state[0], state[1]
        return self.states[configuration][term][prop]

    def _query_state_prop_exist(self, state, prop):
        configuration, term = state[0], state[1]
        return prop in self.states[configuration][term]

    def final_states(self):
        return [s['state_f'] for s in self.transitions]

    def frequency(self, state):
        return self._query_state_prop(state, 'frequency')

    def transition_frequency(self, state_f):
        """
        calculate transition frequency between initial state and state_f
        """
        fi = self.frequency(self.state_i)
        ff = self.frequency(state_f)
        return (ff-fi)

    def transition_frequency_hyperfine(self, state_f, Ff):
        """
        calculate transition frequency between initial state (F=self.Fi) and
        state_f (F=Ff), including hyperfine shift
        """
        hfs_i = self.hyperfine_shift(self.state_i, self.Fi)
        hfs_f = self.hyperfine_shift(state_f, Ff)

        fi = self.frequency(self.state_i) + hfs_i
        ff = self.frequency(state_f) + hfs_f
        return (ff-fi)

    def sorted_transitions_dict(self, hyperfine=False):
        d = []
        if hyperfine:
            for s in self.final_states():
                for F in self._calc_Fs(self._J_from_state(s)):
                    d.append({
                        'state':s,
                        'F': F,
                        'frequency': self.transition_frequency_hyperfine(s,F)})
        else:
            for s in self.final_states():
                    d.append({
                        'state':s,
                        'frequency': self.transition_frequency(s)})
        sd = sorted(d, key=lambda x: x['frequency'])
        for s in sd:
            s['wavelength'] = c/s['frequency']
        return sd

    def _reduced_mat_el_sqd(self, Fi, Ff, Ji, Jf, omega_Jf_Ji, Gamma):
        x = Gamma*3*np.pi*eps0*hbar*c**3/(abs(omega_Jf_Ji)**3)
        y = (2*Jf+1) * (2*Ff+1)
        z = (wigner_6j(Ji, Jf, 1, Ff, Fi, self.I))**2

        return x*y*z

    def _J_from_state(self, state):
        config, term = state
        try:
            J = self.states[config][term]['J']
            return J
        except:
            if len(term)==3:
                S, L, J = term
                J = int(J)
                return J
            else:
                raise Exception('State (%s, %s) has no J provided.'%(config,term))

    def scalar_polarizability(self, lam):
        omega = 2*np.pi*c/lam

        p = 0
        for trans in self.transitions:
            Ji = self._J_from_state(trans['state_i'])
            Jf = self._J_from_state(trans['state_f'])
            omega_Ff_Fi = 2*np.pi*self.transition_frequency(trans['state_f'])
            Gamma = trans['Gamma']

            for Ff in self._calc_Fs(Jf):
                a = 2*omega_Ff_Fi*self._reduced_mat_el_sqd(self.Fi, Ff, Ji, Jf,
                                                           omega_Ff_Fi, Gamma)
                b = 3*hbar*(omega_Ff_Fi**2 - omega**2)

                p = p + a/b

        return p

    def scalar_lightshift(self, lam, laser_intensity=1):
        """
        lightshift in Hz for lambda in nm and laser_intensity in W/cm^2
        """
        laser_intensity_unit = 1/(0.01)**2  # normalization Watt/cm^2
        p = self.scalar_polarizability(lam)
        l = -p/(2*eps0*c)/h*laser_intensity_unit*laser_intensity
        return l

    @staticmethod
    def _hyperfine_shift(F, I, J, hyper_A, hyper_B):
        """
        hyperfine coupling shift in Hz,
        for hyper_A and hyper_B in Hz.
        see Steck (7.134)
        A. Lurio, M. Mandel, and R. Novick. Second-Order Hyperfine and Zeeman Corrections
        for an sl-Configuration. Physical Review, 126(5):1758–1767, June 1962.
        """
        def K(F):
            return F * (F + 1) - I * (I + 1) - J * (J + 1)

        # dipole hyperfine coupling dEa
        dEa = K(F) * hyper_A/2

        # quadrupole hyperfine coupling dEb
        D = 2*I * (2*I - 1) * 2*J * (2*J - 1)
        if D == 0:
            dEb = 0
        else:
            dEb = (3 * K(F) * (K(F)+1) / 2 - 2*I * (I+1) * J * (J+1)) / D * hyper_B
        return (dEa + dEb)

    def hyperfine_shift(self, state, F):
        J = self._J_from_state(state)

        if self._query_state_prop_exist(state, 'hyper_A'):
            hyper_A = self._query_state_prop(state, 'hyper_A')
        else:
            hyper_A = 0
        if self._query_state_prop_exist(state, 'hyper_B'):
            hyper_B = self._query_state_prop(state, 'hyper_B')
        else:
            hyper_B = 0

        return self._hyperfine_shift(F, self.I, J, hyper_A, hyper_B)


    def polarizabilities(self, lam):
        """
        Scalar, vector and tensor polarizability in SI units
        """
        omega = 2*np.pi*c/lam
        Fi = self.Fi

        p0, p1, p2 = (0,0,0)

        for trans in self.transitions:
            state_f = trans['state_f']
            Ji = self._J_from_state(trans['state_i'])
            Jf = self._J_from_state(trans['state_f'])
            Gamma = trans['Gamma']

            for Ff in self._calc_Fs(Jf):
                omega_Ff_Fi = 2*np.pi*self.transition_frequency_hyperfine(state_f, Ff)
                matel = self._reduced_mat_el_sqd(self.Fi, Ff, Ji, Jf,
                                                 omega_Ff_Fi, Gamma)
                f = omega_Ff_Fi*matel/(hbar*(omega_Ff_Fi**2 - omega**2))

                p0 += 2/3 * f
                p1 += (-1)**(Fi+Ff+1) * np.sqrt(6*Fi*(2*Fi+1)/(Fi+1))\
                      * np.float(wigner_6j(1,1,1,Fi,Fi,Ff)) * f
                p2 += (-1)**(Fi+Ff) * np.sqrt(40*Fi*(2*Fi+1)*(2*Fi-1)/(3*(Fi+1)*(2*Fi+3)))\
                      * np.float(wigner_6j(1,1,2,Fi,Fi,Ff)) * f

        return p0,p1,p2


    def lightshifts(self, lam, q, mFi, laser_intensity=1.):
        """
        scalar, vector and tensor lightshift in Hz for lambda in nm and laser_intensity in W/cm^2
        """
        laser_intensity_unit = 1/(0.01)**2  # normalization Watt/cm^2

        Fi = self.Fi
        p0,p1,p2 = self.polarizabilities(lam)

        E02 = laser_intensity*laser_intensity_unit/(2*eps0*c)/h

        # scalar light shift
        l0 = -p0*E02
        l1 = -p1*mFi/Fi*(1)*q*E02
        l2 = -p2*(3*mFi**2-Fi*(Fi+1))/(Fi*(2*Fi-1))*(3*(1-np.abs(q))*E02-E02)/2

        return l0,l1,l2

    def total_lightshift(self, lam, q, mFi, laser_intensity=1.):
        """
        Summed scalar, vector and tensor lightshift in Hz
        for lambda in nm and laser_intensity in W/cm^2
        """
        l0,l1,l2 = self.lightshifts(lam, q, mFi, laser_intensity=laser_intensity)
        return l0+l1+l2
