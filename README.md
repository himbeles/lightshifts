# Lightshifts

Calculate dynamical scalar, vector and tensor light shifts for atomic states
in the presence of hyperfine coupling.

The module ```lightshifts.lightshift_solver``` solves for the scalar, 
vector and tensor light shifts induced by atomic dipole transitions.

State energies, hyperfine coupling coefficients and atomic transition properties
are provided in the form of json files.


## Installation

Install the module using

```bash
pip install git+https://gitlab.physik.uni-muenchen.de/Luis.Riegger/lightshifts.git
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

-   `atom-X.json` is a database of all the states connected by the atomic transitions that induce the dynamic light shift. It also includes the `name` of the atom and its nuclear spin `I`. The states are specified and sorted using an electon shell configuration string (e.g. `6s6p`) and the LS coupling name of the state (e.g. `3P1`). From the LS coupling name, the tool infers the relevant quantum number `J` for the given state. If only the quantum number J is known for the state, you can give it an arbitrary name and specify `J` as a property of the state (see [```examples/atom_yb173.json```](examples/atom_yb173.json)). The `frequency` in Hz of all the states must be provided relative to the ground state or *one* reference state). 
   
   Optionally, you can give the hyperfine coupling coefficients A (`hyper_A`) and B (`hyper_B`) in Hz, such that the vector and tensor shifts around these states can be calculated.


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


