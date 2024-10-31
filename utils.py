"""
Some helper functions
"""

from sage.all import matrix, block_matrix, Integer

def pole_number_table(f_pole_numbers, g_pole_numbers):
    """Compute the pole number table"""

    return matrix(len(f_pole_numbers), len(g_pole_numbers), lambda i, j: f_pole_numbers[i] + g_pole_numbers[j])

def choose_function_with_pole(E, n):
    """Construct a function with pole divisor n*Pinf"""

    assert n >= 0, "n has to be in the Weierstrass semigroup"

    # E is assumed to be a proper extension of a rational function field
    y, = E.gens()
    x, = E.base_field().gens()
    x = E(x)

    d = Integer(2 * E.genus() + 1)  # needs to be converted to integer

    Pinf = E.places_infinite(1)[0]

    assert x.divisor_of_poles().multiplicity(Pinf) == 2
    assert y.divisor_of_poles().multiplicity(Pinf) == d

    if n % 2 == 0:
        return x**(n // 2)

    assert n - d >= 0, "n has to be in the Weierstrass semigroup"

    return x**((n - d) // 2) * y

def generator_matrix(basis, places):
    """Construct the canonical generator matrix given a basis and some places"""

    return matrix([[f.evaluate(P) for P in places] for f in basis])

def vector_matrix_multiplication(vec, mat):
    """Compute the vector matrix product"""

    return [sum(v * c for v, c in zip(vec, col)) for col in mat.columns()]

def vsplit(A, K):
    """Split the matrix A to K pieces vertically"""

    Kn, _ = A.dimensions()

    assert Kn % K == 0

    n = Kn // K

    return [A[k*n:(k+1)*n] for k in range(K)]

def hsplit(B, L):
    """Split the matrix B to L pieces horizontally"""

    _, Ln = B.dimensions()

    assert Ln % L == 0

    n = Ln // L

    return [B[:, l*n:(l+1)*n] for l in range(L)]

def assemble_block_matrix_in_column_order(submatrices, K, L):
    """Assemble a sequence of blocks to a block matrix in column major order"""

    blocks = [[None for _ in range(L)] for _ in range(K)]

    i, j = 0, 0
    for block in submatrices:
        blocks[i][j] = block
        i += 1
        if i == K:
            i = 0
            j += 1

    return block_matrix(blocks, subdivide=False)
