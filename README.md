# lightshifts

<img src="misc/icon.png" alt="Project Icon" width="64" height="64">

Calculate dynamical scalar, vector and tensor light shifts for atomic states
in the presence of hyperfine coupling.

The module ```lightshifts.lightshift_solver``` solves for the scalar, 
vector and tensor light shifts induced by atomic dipole transitions.

State energies, hyperfine coupling coefficients and atomic transition properties
are provided in the form of json files.


## Installation

Install from PyPI

```sh
pip install lightshifts
```

## Development

For package development, set up development environment with 

```sh
uv sync --dev
```

### Tests

Run tests via 
```sh
uv run pytest 
```


## Usage of ```lightshifts.lightshift_solver```

This is a basic example to calculate the scalar dynamic light shift induced in the groundstate of ytterbium-173 by 100 W/cm^2 light with a  wavelength 670 nm using the example state file [```examples/atom_yb173.json```](examples/atom_yb173.json) and transition file [```examples/transitions_1S0.json```](examples/transitions_1S0.json).

```python
# Import the module using 
import lightshifts.lightshift_solver as ls_solver

# Load the atomic state properties file and the transitions file
ls = ls_solver('examples/atom_yb173.json', 'examples/transitions_1S0.json')

# Calculate the scalar lightshift in Hz
ls.scalar_lightshift(lam=670e-9, laser_intensity=100)
# Out[]: -960.039953941305
```

The total dynamic light shift, including scalar, vector and tensor light shifts, can also be obtained -- here, for a magnetic sublevel $`m_F=-5/2`$ and $`\sigma^+`$ polarized light ($`q=1`$):

```python
# Calculate the total lightshift in Hz
ls.total_lightshift(lam=670e-9, q=1, mFi=-5/2, laser_intensity=100)
# Out[]: -960.038498669505

# or a list holding the scalar, vector and tensor light shift separately
ls.lightshifts(lam=670e-9, q=1, mFi=-5/2, laser_intensity=100)
# Out[]: (-960.039953992523, 0.000817922762313870, 0.000637400256046666)
```

For an example of more advanced module usage, see [```examples/example_lightshifts_yb173.ipynb```](examples/example_lightshifts_yb173.ipynb).

## The state and transition files

State energies, hyperfine coupling coefficients and transition properties for an atom 
are provided in the form of two json files:

-   `atom.json` is a database of (at least) all the states connected by the atomic transitions that induce the dynamic light shift and which are given in the second file. It also includes the `name` of the atom and its nuclear spin `I`. The states are specified and sorted using an electon shell configuration string (e.g. `6s6p`) and the LS coupling name of the state (e.g. `3P1`). From the LS coupling name, the tool infers the relevant quantum number `J` for the given state. If only the quantum number J is known for the state, you can give it an arbitrary name and specify `J` as a property of the state (see [```examples/atom_yb173.json```](examples/atom_yb173.json)). The `frequency` in Hz of all the states must be provided relative to the ground state or *one* reference state). 
   
    Optionally, you can give the hyperfine coupling coefficients A (`hyper_A`) and B (`hyper_B`) in Hz, such that the vector and tensor shifts around these states can be calculated.

    For example:

    ```json
    {
        "name": "yb173",
        "I": 2.5,
        "states": {
            "6s6s": {
                "1S0": {
                    "frequency": 0,
                    "hyper_A": 0,
                    "hyper_B": 0
                }
            },
            "6s6p": {
                "3P0": {
                    "frequency": 518294362279306.1,
                    "_ref_frequency": "NIST",
                    "hyper_A": 0,
                    "hyper_B": 0
                },
                "3P1": {
                    "frequency": 539386800288320.6,
                    "_ref_frequency": "NIST",
                    "hyper_A": -1094328000.0,
                    "hyper_B": -826635000.0,
                    "_ref_hyper": "Pandey et al., PRA 80, 022518 (2009)"
                }
            }
        }
    }
    ```

-   `transitions-stateX.json` lists all the relevant transitions from **one** starting (initial) state `state_i`, that we want to calculate the light shifts for. It should include the transitions most relevant for the wavelength range to be probed, meaning transitions to final states `state_f` with a wavelength near that range and the broader the transition the more important. As the transition strength, specify the decay rate `Gamma` from `state_f` to `state_i` (Einstein A coefficient). For a closed transition, this would be the inverse lifetime of `state_f` or $`\Gamma = 2\pi\gamma`$ with $`\gamma`$ the natural linewidth.

    For example:

    ```json
    [
        {
            "state_i": [
                "6s6s",
                "1S0"
            ],
            "state_f": [
                "6s6p",
                "1P1"
            ],
            "Gamma": 183016105.4172767,
            "_ref_Gamma": "Blagoev-1994"
        },
        {
            "state_i": [
                "6s6s",
                "1S0"
            ],
            "state_f": [
                "6s6p",
                "3P1"
            ],
            "Gamma": 1154601.0853250204,
            "_ref_Gamma": "Blagoev-1994"
        }
    ]
    ````

## Estimate branching ratios

The reduced dipole matrix element can be calculated from a measured transition rate between two LS coupling states. However, sometimes only the lifetime of a state is known experimentally and not the branching ratios into energetically lower lying states. 
The method  ```branching_ratio_LS```  of the class ```lightshifts.atom``` can help by estimating the ratio of dipole matrix elements between a selection of states, using the parity selection rule for the electron configuration and angular momentum selection.

First, import a dictionary of atomic states and their energies (same as the atomic states file above).

```python
import lightshifts.atom as atom
yb = atom.from_json('atom_yb173.json')
```

Then, calculate all branching ratios of an initial state state_i into all energetically lower lying states in the imported state library:

```python
state_i = ('6s5d','3D1')
yb.branching_ratios_LS_dict(state_i, verbose=True)

# Out[]:
#    branching ratio into  ('6s6p', '3P0') = 0.6387527684341578
#    branching ratio into  ('6s6p', '3P1') = 0.3519121426242965
#    branching ratio into  ('6s6p', '3P2') = 0.009335088941545725
#
#    {('6s6s', '1S0'): 0.0,
#     ('6s6p', '3P0'): 0.6387527684341578,
#     ('6s6p', '3P1'): 0.3519121426242965,
#     ('6s6p', '3P2'): 0.009335088941545725}
```

or, calculate the branching into one single final state only: 

```python
state_i = ('6s5d','3D1')
state_f = ('6s6p','3P0')
yb.branching_ratio_LS(state_i, state_f)

# Out[]: 0.6387527684341578
```

An example of how the transition rates were calculated for ytterbium-173 can be found in [```examples/example_transition_collection_yb173.ipynb```](examples/example_transition_collection_yb173.ipynb).
