A very basic Molecular Orbital calculator. Kind of utilises the Hartree-Fock approximation.

Works by defining two matrices, H and S. 

    H[x][y] =  integrate wavefunction[x] * hamiltonian(wavefunction[y]) over all space
    S[x][y] =  integrate wavefunction[x] * wavefunction[y] over all space

Then we form a set of relationships;

    det(H-ES) = 0

This is then solved to give a number of solutions for the energy, E. These E values are then used to calculate the orbital coefficients (which are often imprecise).

For example;

    H
    -0.5008	-0.4343	
    -0.4343	-0.5008	
    S
    1.0026	0.9293	
    0.9293	1.0026	
    -0.908244 - [0.728000	0.685577]
    -0.483250 - [-0.999000	0.044710]

Here the model has calculated H and S for two 1s orbitals on hydrogen atoms at a distance of 1.4ao, and has then determined two energies - -0.908 which corresponds to 0.728Sa + 0.686Sb (Sa and Sb are the AOs on the respective atoms) (this is very close to what it should be, eg (Sa + Sb)/sqrt(2)) and -0.483 which it predicts as corresponding to -0.999Sa + 0.045Sb - this is incorrect and will be fixed in a later model.

MIT licence.
