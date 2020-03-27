## Markov Chain Monte Carlo for the Canonical Ensemble

Semi-interactive code for simulations in the *NVT* ensemble according to the Metropolis' acceptance/rejection criteria derived from the detailed balance condition of the ergodic markov chain proposed.

The environment consists of a fixed-size simulation box centered in the origin, filled with spherical molecules based on the density of the system. Their interaction is according to either a Lennard-Jones or Square-Well potential (implementation of additional potentials is straightforward) in the reduced unit's system.

![Acceptance Criteria]{img/AcceptanceNVT.png}

    julia Canonical.jl &rho T