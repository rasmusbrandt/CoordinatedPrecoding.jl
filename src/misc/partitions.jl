##########################################################################
# Helper functions for set partitions.

# A block is just a set of integers.
immutable Block
    elements::IntSet
end
Block() = Block(IntSet())

# Inherit part of the interface from IntSet
Base.show(io::IO, b::Block) = showcompact(io, b.elements)
Base.length(b::Block) = length(b.elements)
Base.push!(b::Block, e::Int) = push!(b.elements, e)

# The Partition type describes a partition of integers into blocks.
immutable Partition
    blocks::Set{Block}
end
Partition() = Partition(Set{Block}())

# Inherit part of the interface from Set
Base.show(io::IO, p::Partition) = showcompact(io, p.blocks)
Base.length(p::Partition) = length(p.blocks)

# Create partition from restricted growth string
function Partition(a::Vector; skip_check::Bool=false)
    I = length(a)
    a_max = maximum(a)
    num_blocks = 1 + a_max

    # Consistency checks
    if !skip_check
        all(sort(unique(a)) .== 0:a_max) || Lumberjack.error("Restricted growth string may not skip values.", @Compat.Dict( :a => a ))
    end

    # Build blocks and partition
    blocks = [ Block() for n = 1:num_blocks ]
    for l = 1:I
        push!(blocks[a[l] + 1], l)
    end

    return Partition(Set(blocks))
end

# Create partition from assignment matrix
Partition(A::Matrix; skip_check::Bool=false) =
    Partition(restricted_growth_string(A); skip_check=skip_check)

# Convert from partition to restricted growth string
restricted_growth_string(p::Partition) = restricted_growth_string(logical_matrix(p))

# Convert from partition to assignment matrix
function logical_matrix(p::Partition)
    I = maximum([ maximum(block.elements) for block in p.blocks ])
    A = eye(Int, I, I)

    for block in p.blocks
        elements = collect(block.elements)
        A[elements, elements] = 1
    end

    return A
end

# Convert from assignment matrix to restricted growth string
function restricted_growth_string(A::Matrix)
    I = size(A, 1)

    # Consistency checks
    A' == A || Lumberjack.error("Tried to create partition from asymmetric assignment matrix.", @Compat.Dict(:A => A))
    diag(A) == ones(I) || Lumberjack.error("Elements must belong to their own block.", @Compat.Dict(:A => A))

    # Build restricted growth string
    Ac = copy(A)
    a = zeros(Int, I)
    block_no = 0
    for i = 1:I
        # Find all integers in the block containing integer i.
        elements = find(Ac[i,:] .== 1)

        if length(elements) > 0
            # Consistency check
            for j in elements[2:end]
                find(Ac[j,:] .== 1) == elements || Lumberjack.error("Incorrect structure in assignment matrix.", @Compat.Dict(:A => Ac))
            end

            # Add all elements to block and increment number of blocks.
            a[elements] = block_no
            block_no += 1

            # Do not investigate these elements more
            Ac[elements,elements] = 0
        end
    end

    return a
end

# Convert from restricted growth string to assignment matrix
function logical_matrix(a::Vector)
    I = length(a)
    a_max = maximum(a)
    num_blocks = 1 + a_max

    # Consistency checks
    all(sort(unique(a)) .== 0:a_max) || Lumberjack.error("Restricted growth string may not skip values.", @Compat.Dict(:a => a))

    # Create matrix iteratively
    A = eye(Int, I, I)
    for n = 0:(num_blocks-1)
        elements = find(a .== n)
        A[elements, elements] = 1
    end

    return A
end

# PartitionIterator gives all partitions of the set 1:n. This iterator is
# based on Algorithm H from TAoCP 7.2.1.5. Note that a is zero-based.
immutable PartitionIterator
    n::Int

    PartitionIterator(n) = n >= 2 ? new(n) : error("only applicable for n >= 2.")
end

type PartitionIteratorState
    a::Vector{Int}
    b::Vector{Int}
    m::Int
    terminate::Bool
end
PartitionIteratorState(n) = PartitionIteratorState(zeros(n), ones(n - 1), 1, false)

Base.start(p::PartitionIterator) = PartitionIteratorState(p.n) # H1
Base.done(p::PartitionIterator, s) = s.terminate
function Base.next(p::PartitionIterator, s)
    # H2
    current_partition = Partition(s.a, skip_check=true) # By construction, a is a valid restricted growth string.

    if s.a[p.n] < s.m
        # H3
        s.a[p.n] += 1
    else
        # H4
        j = p.n - 1
        while s.a[j] == s.b[j]
            j -= 1
        end

        # H5
        if j == 1
            s.terminate = true
        else
            s.a[j] += 1

            # H6
            s.m = s.b[j] + (s.a[j] == s.b[j] ? 1 : 0)
            j += 1
            while j < p.n
                s.a[j] = 0
                s.b[j] = s.m
                j += 1
            end
            s.a[p.n] = 0
        end
    end

    # H2 cont'd.
    return current_partition, s
end

##########################################################################
# Stirling and Bell numbers

# Stirling2NumberCoefficientIterator gives all coefficients that when summed
# result in a Stirling number of the second kind. Note that the coefficients
# may be real numbers, but when summed are integers. The coefficients are
# thus calculated as floating point numbers.
immutable Stirling2NumberCoefficientIterator
    n::BigInt
    k::BigInt
    k_fac::BigInt
end
Stirling2NumberCoefficientIterator(n, k) = Stirling2NumberCoefficientIterator(n, k, factorial(k))
Base.start(s::Stirling2NumberCoefficientIterator) = BigInt(0) # state: j
Base.done(s::Stirling2NumberCoefficientIterator, j) = (j > s.k)
Base.next(s::Stirling2NumberCoefficientIterator, j) = ((iseven(s.k - j) ? 1. : -1.)*(binomial(s.k, j)*j^s.n/s.k_fac), j + 1)

# Stirling2NumberIterator gives all Stirling numbers of the second kind
# up to n. The Stirling numbers are integers, but the Stirling number
# coefficients are not. We therefore round the sum of the coefficients,
# in order to end up with an integer.
immutable Stirling2NumberIterator
    n::BigInt
end
Base.start(s::Stirling2NumberIterator) = BigInt(0) # state: k
Base.done(s::Stirling2NumberIterator, k) = (k > s.n)
Base.next(s::Stirling2NumberIterator, k) = (BigInt(round(sum(Stirling2NumberCoefficientIterator(s.n, k)))), k + 1)

# The Bell number can we calculated as the sum of all
# Stirling numbers of the second kind.
bell_number(n) = sum(Stirling2NumberIterator(n))
