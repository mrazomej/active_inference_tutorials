# Import basic math
import StatsBase
import Distributions
import Random
import LinearAlgebra

## =============================================================================
# Monte Carlo simulation of the Ising model on a rectangular lattice
## =============================================================================

@doc raw"""
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

@doc raw"""
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

@doc raw"""
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

## =============================================================================
# Transfer Matrix Method for the Ising Model on a Rectangular Lattice
## =============================================================================

@doc raw"""
    _generate_spin_states(N::Int)

Generate all possible spin states for a 1D system of N spins.

# Arguments
- `N`: the number of spins in the system.

# Returns
- A 2D array of size N x 2^N, where each column represents a possible spin
  state. Each element is +1 or -1, representing spin up or spin down,
  respectively.

# Notes
- The function generates all 2^N possible states of the system by looping over
  the integers from 0 to 2^N - 1 and converting each integer to its binary
  representation. The binary representation is then converted to spin states,
  with 0 representing spin down (-1) and 1 representing spin up (+1).
"""
function _generate_spin_states(N::Int)
    # Define the number of states
    num_states = 2^N
    # Initialize array to store the spin states
    states = Array{Int8,2}(undef, N, num_states)

    # Loop through each state
    for i in 0:(num_states-1)
        # Store the binary representation of the state
        binary = bitstring(i)
        # Loop through each spin
        for j in 1:N
            # Convert the binary representation to spin states
            states[j, i+1] = (
                j <= length(binary) && binary[end-j+1] == '1'
            ) ? 1 : -1
        end # for j
    end # for i

    return states
end

## =============================================================================

@doc raw"""
    generate_pair_matrix(n_row::Int)

Generate a matrix of all possible pairs of column states for a rectangular
lattice of spins.

# Arguments
- `n_row`: the number of rows in the lattice.

# Returns
- A 2D array of matrices, where each matrix represents a pair of column states.
  Each element in the matrices is +1 or -1, representing spin up or spin down,
  respectively.

# Notes
- The function first generates all possible column states using the
  `generate_spin_states` function. Each column state is a vector of length
  `n_row`, where each element is a spin.
- It then creates a matrix where each element is a pair of column states. The
  pairs are generated by taking all possible combinations of the column states.
- The pair of column states is represented as a matrix where the first column is
  the first state and the second column is the second state.
"""
function _generate_pair_matrix(n_row::Int)
    # Generate possible column states
    states = _generate_spin_states(n_row)

    # Get number of column states
    num_states = size(states, 2)

    # Initialize pair matrix to store all possible pair of column states
    pair_matrix = Array{Matrix{Int8}}(undef, num_states, num_states)

    # Loop through rows
    for i in 1:num_states
        # Loop through columns
        for j in 1:num_states
            # Generate pair of column states
            pair_matrix[i, j] = hcat(states[:, i], states[:, j])
        end # for j
    end # for i

    return pair_matrix
end

## =============================================================================

@doc raw"""
    ε_rectangular(pair_state::AbstractMatrix; J::T=1.0) where {T<:Number}

Calculate the energy of two column states in the Ising model on a rectangular
lattice. This is the energy term used to build the transition matrix.

# Arguments
- `pair_state`: a 2D array representing the spin configuration of the lattice.
  Each element is +1 or -1, representing spin up or spin down, respectively.

# Optional Keyword Arguments
- `J`: the coupling constant in the Ising model. It determines the strength of
  the interaction between neighboring spins. Default is 1.0.

# Returns
- The total energy of the spin configuration.

# Notes
- The function calculates the energy by summing over all pairs of neighboring
  spins. For each pair, it adds the product of the spins and the coupling
  constant to the total energy. The interaction with the right and bottom
  neighbors is considered for each spin.
- The energy is negative when neighboring spins are aligned, and positive when
  they are anti-aligned. This is because the Ising model assumes that spins
  prefer to align with their neighbors.
"""
function ε_rectangular(pair_state::AbstractMatrix; J=1.0)
    # Initialize energy to zero
    E = zero(typeof(J))
    # Add the interaction with the right neighbor
    E += sum(J .* pair_state[:, 1] .* pair_state[:, 2])
    # Add the interaction with bottom neighbor
    E += sum(J / 2 .* pair_state[1:end-1, :] .* pair_state[2:end, :])

    return E
end # function

## =============================================================================

@doc raw"""
        transfer_matrix(n_row::Int; β::T=1.0, J::T=1.0) where {T<:Number}

Generate the transfer matrix for a rectangular lattice of spins in the Ising
model.

# Arguments
- `n_row::Int`: the number of rows in the lattice.

# Optional Keyword Arguments
- `β=1.0`: the inverse temperature in the Ising model. Default is 1.0.
- `J=1.0`: the coupling constant in the Ising model. It determines the strength
  of the interaction between neighboring spins. Default is 1.0.

# Returns
- The transfer matrix, a symmetric matrix where each element is the Boltzmann
    factor exp(-βE) of a pair of column states. E is the energy of the pair
    state, calculated using the `ε_rectangular` function.

# Notes
- The function first generates a matrix of all possible pairs of column states
  using the `generate_pair_matrix` function.
- It then calculates the energy of each pair state using the `ε_rectangular`
  function and assigns the Boltzmann factor exp(-βE) to the corresponding
  element of the transfer matrix.
- The transfer matrix is used in the transfer-matrix method, a powerful
  technique for solving statistical mechanics problems on lattices.
"""
function transfer_matrix(n_row::Int; β=1.0, J=1.0)
    # Generate pair matrix
    pair_matrix = _generate_pair_matrix(n_row)

    # Get number of column states
    num_states = size(pair_matrix, 1)

    # Initialize transfer matrix
    M = Matrix{typeof(β)}(undef, num_states, num_states)

    # Loop through rows
    for i in 1:num_states
        # Loop through columns up to current row
        for j in 1:i
            # Calculate energy of pair state
            E = ε_rectangular(pair_matrix[i, j], J=J)
            # Assign energy to transfer matrix
            M[i, j] = exp(-β * E)
        end # for j
    end # for i

    # Return matrix as a symmetric matrix
    return LinearAlgebra.Symmetric(M, :L)
end # function

## =============================================================================

"""
    logZ_rectangular(n_row::Int, n_col::Int; β::T=1.0, J::T=1.0) where {T<:Number}

Calculate the logarithm of the partition function for a rectangular lattice of
spins in the Ising model.

# Arguments
- `n_row`: the number of rows in the lattice.
- `n_col`: the number of columns in the lattice.

# Optional Keyword Arguments
- `β`: the inverse temperature in the Ising model. Default is 1.0.
- `J`: the coupling constant in the Ising model. It determines the strength of
  the interaction between neighboring spins. Default is 1.0.

# Returns
- The logarithm of the partition function.

# Notes
- The function first generates the transfer matrix using the `transfer_matrix`
  function.
- It then computes the partition function by summing the eigenvalues of the
  transfer matrix raised to the power of the number of columns, and takes the
  logarithm of the result.
- The partition function is a key quantity in statistical mechanics. It encodes
  the statistical properties of a system in equilibrium. Its logarithm is
  related to the free energy of the system.
"""
function logZ_rectangular(n_row::Int, n_col::Int; β=1.0, J=1.0)
    # Generate transfer matrix
    M = transfer_matrix(n_row, β=β, J=J)

    # Compute the partition function
    Z = sum(LinearAlgebra.eigen(M).values .^ n_col)

    return log(Z)
end

## =============================================================================

@doc raw"""
    entropy_rectangular(
        n_row::Int, n_col::Int; 
        β::T=1.0, J::T=1.0, k::T=1.0, ε::T=cbrt(eps(Float64))
    ) where {T<:Number}

Calculate the entropy for a rectangular lattice of spins in the Ising model.

# Arguments
- `n_row`: the number of rows in the lattice.
- `n_col`: the number of columns in the lattice.

# Optional Keyword Arguments
- `β`: the inverse temperature in the Ising model. Default is 1.0.
- `J`: the coupling constant in the Ising model. It determines the strength of
  the interaction between neighboring spins. Default is 1.0.
- `k`: the Boltzmann constant. Default is 1.0.
- `ε`: a small number used for finite difference approximation. Default is the
  cube root of the machine epsilon for Float64.

# Returns
- The entropy of the system.

# Notes
- The function first computes the logarithm of the partition function using the
  `logZ_rectangular` function.
- It then computes the derivative of the logarithm of the partition function
  with respect to β using a finite difference approximation.
- The entropy is then calculated using the formula S = k * logZ + k * β *
  ∂logZ∂β, where S is the entropy, k is the Boltzmann constant, logZ is the
  logarithm of the partition function, β is the inverse temperature, and ∂logZ∂β
  is the derivative of logZ with respect to β.
"""
function entropy_rectangular(
    n_row::Int, n_col::Int; β=1.0, J=1.0, k=1.0, ε=cbrt(eps(Float64))
)
    # Compute log partition function
    logZ = logZ_rectangular(n_row, n_col; β=β, J=J)

    # Compute finite difference derivative of log partition function with
    # respect to β
    ∂logZ∂β = (
        logZ_rectangular(n_row, n_col; β=β + ε) -
        logZ_rectangular(n_row, n_col; β=β - ε)
    ) / (2 * ε)

    # Compute entropy
    S = k * logZ - k * β * ∂logZ∂β

    return S
end

## =============================================================================

@doc raw"""
    heat_capacity_rectangular(
        n_row::Int, n_col::Int; 
        β::T=1.0, J::T=1.0, k::T=1.0, ε::T=cbrt(eps(Float64))
    ) where {T<:Number}

Calculate the heat capacity for a rectangular lattice of spins in the Ising
model.

# Arguments
- `n_row`: the number of rows in the lattice.
- `n_col`: the number of columns in the lattice.

# Optional Keyword Arguments
- `β`: the inverse temperature in the Ising model. Default is 1.0.
- `J`: the coupling constant in the Ising model. It determines the strength of
  the interaction between neighboring spins. Default is 1.0.
- `k`: the Boltzmann constant. Default is 1.0.
- `ε`: a small number used for finite difference approximation. Default is the
  cube root of the machine epsilon for Float64.

# Returns
- The heat capacity of the system.

# Notes
- The function first computes the second derivative of the logarithm of the
  partition function with respect to β using a finite difference approximation.
- The heat capacity is then calculated using the formula C = k * β^2 *
  ∂²logZ∂β², where C is the heat capacity, k is the Boltzmann constant, β is the
  inverse temperature, and ∂²logZ∂β² is the second derivative of logZ with
  respect to β.
"""
function heat_capacity_rectangular(
    n_row::Int, n_col::Int; β=1.0, J=1.0, k=1.0, ε=cbrt(eps(Float64))
)
    # Compute second order derivative of log partition function with respect to
    # β
    ∂²logZ∂β² = (
        logZ_rectangular(n_row, n_col; β=β + ε) -
        2 * logZ_rectangular(n_row, n_col; β=β) +
        logZ_rectangular(n_row, n_col; β=β - ε)
    ) / ε^2

    # Compute heat capacity
    C = k * β^2 * ∂²logZ∂β²

    return C
end