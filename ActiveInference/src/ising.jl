# Import basic math
import StatsBase
import Distributions
import Random

## =============================================================================

"""
    E_rectangular(x̲::Matrix, J::Number, H::Number)

Compute the energy of a configuration in the Ising model on a rectangular
lattice.

# Arguments
- `x̲::Matrix`: a matrix representing the spin configuration. Each element is +1
  or -1.

# Optional Keyword Arguments
- `J`: the coupling constant. Positive for ferromagnetic coupling, negative for
  antiferromagnetic.
- `H`: the external magnetic field.

# Returns
- The energy of the configuration.

# Notes
- The function uses periodic boundary conditions.
- The energy is computed as -J ∑ x̲[i, j] x̲[i', j'] - H ∑ x̲[i, j], where the
  first sum is over all pairs of neighboring spins and the second sum is over
  all spins.
"""
function E_rectangular(x̲::Matrix; J=1.0f0, H=0.0f0)
    # Get the size of the spin configuration matrix
    m, n = size(x̲)
    # Initialize the energy to zero
    energy = 0.0f0

    # Compute the coupling energy between neighboring spins
    for i in 1:m
        for j in 1:n
            # Apply periodic boundary conditions
            # If i > 1, then im1 = i - 1, else im1 = m (the last row)
            im1 = i > 1 ? i - 1 : m
            # If i < m, then ip1 = i + 1, else ip1 = 1 (the first row)
            ip1 = i < m ? i + 1 : 1
            # If j > 1, then jm1 = j - 1, else jm1 = n (the last column)
            jm1 = j > 1 ? j - 1 : n
            # If j < n, then jp1 = j + 1, else jp1 = 1 (the first column)
            jp1 = j < n ? j + 1 : 1

            # Subtract the energy contribution from the neighbor below
            energy -= J * x̲[i, j] * x̲[ip1, j]
            # Subtract the energy contribution from the neighbor to the right
            energy -= J * x̲[i, j] * x̲[i, jp1]
        end # for j
    end # for i

    # Subtract the energy contribution from the external field
    energy -= H * sum(x̲)

    # Return the total energy
    return energy
end

## =============================================================================

"""
    local_field_rectangular(n::Int, x̲::Matrix, J, H)

Compute the local magnetic field at a given spin in the Ising model on a
rectangular lattice.

# Arguments
- `n::Int`: the linear index of the spin.
- `x̲::Matrix`: a matrix representing the spin configuration. Each element is +1
  or -1.

# Optional Keyword Arguments
- `J`: the coupling constant. Positive for ferromagnetic coupling, negative for
  antiferromagnetic.
- `H`: the external magnetic field.

# Returns
- The local magnetic field at the spin.

# Notes
- The function uses periodic boundary conditions.
- The local field is computed as J ∑ x̲[i', j'] + H, where the sum is over the
  neighboring spins of the given spin.
"""
function local_field_rectangular(
    n::Int, x̲::Matrix; J=1.0f0, H=0.0f0
)
    # Get the size of the spin configuration matrix
    m, n = size(x̲)

    # Convert the linear index to row and column indices
    i, j = Tuple(CartesianIndices(x̲)[n])

    # Apply periodic boundary conditions
    # If i > 1, then im1 = i - 1, else im1 = m (the last row)
    im1 = i > 1 ? i - 1 : m
    # If i < m, then ip1 = i + 1, else ip1 = 1 (the first row)
    ip1 = i < m ? i + 1 : 1
    # If j > 1, then jm1 = j - 1, else jm1 = n (the last column)
    jm1 = j > 1 ? j - 1 : n
    # If j < n, then jp1 = j + 1, else jp1 = 1 (the first column)
    jp1 = j < n ? j + 1 : 1

    # Compute the local field at the spin
    # This is the sum of the spins of the neighboring spins, multiplied by the coupling constant, plus the external field
    bₙ = J * (x̲[im1, j] + x̲[ip1, j] + x̲[i, jm1] + x̲[i, jp1]) + H

    # Return the local field
    return bₙ
end

## =============================================================================

"""
    gibbs_sampling!(
        x̲::Matrix{<:Number}, n_iter::Int; 
        J::Number=1, H::Number=0, β::Number=1, 
        local_field::Function=local_field_rectangular
    )

Perform Gibbs sampling on a spin configuration in the Ising model.

# Arguments
- `x̲::Matrix`: a matrix representing the initial spin configuration. Each
  element is +1 or -1.
- `n_iter::Int`: the number of iterations to perform.

## Optional Keyword Arguments
- `J::Number`: the coupling constant. Positive for ferromagnetic coupling,
  negative for antiferromagnetic.
- `H::Number`: the external magnetic field.
- `β::Number`: the inverse temperature.
- `local_field::Function`: a function that computes the local magnetic field at
  a given spin.

# Returns
- The final spin configuration after `n_iter` iterations of Gibbs sampling.

# Notes
- The function modifies the input `x̲` in place.
- The function uses the Gibbs sampling algorithm, which is a type of Markov
  chain Monte Carlo (MCMC) algorithm.
- The local field is computed using the `local_field` function, which should
  take as arguments the linear index of the spin, the spin configuration matrix,
  the coupling constant, and the external field.
"""
function gibbs_sampling!(
    x̲::Matrix,
    n_iter::Int;
    J=1.0f0,
    H=0.0f0,
    β=1.0f0,
    local_field::Function=local_field_rectangular
)
    # Get the size of the spin configuration matrix
    m, n = size(x̲)
    # Compute the total number of spins
    N = m * n

    # Perform Gibbs sampling for `n_iter` iterations
    for iter in 1:n_iter
        # Select a random spin
        spin_idx = rand(1:N)

        # Compute the local field acting on the selected spin
        b_n = local_field(spin_idx, x̲; J=J, H=H)

        # Compute the probability of the spin being +1
        # This is given by the logistic function of the local field, scaled by the inverse temperature
        p_up = 1 / (1 + exp(-2 * β * b_n))

        # Generate a random number
        r = rand()

        # Update the spin according to the generated number and the computed
        # probability If the number is less than the probability, the spin is
        # set to +1, otherwise it is set to -1
        if r < p_up
            x̲[spin_idx] = 1
        else
            x̲[spin_idx] = -1
        end # if r < p_up
    end # for iter

    # Return the final spin configuration
    return x̲
end # function gibbs_sampling!
