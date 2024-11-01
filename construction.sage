from utils import (
    pole_number_table,
    generator_matrix,
    vsplit,
    hsplit,
    vector_matrix_multiplication,
    assemble_block_matrix_in_column_order,
    choose_function_with_pole,
    )

from generate_plot import GASP, A3S, GASP_big, GASP_small

# STEP 1: Set up parameters (these can be changed)
K = L = 4
X = 4
q = 127

assert K % 2 == 0 and K >= 2, "K has to be even and positive"
assert L >= 1, "L has to be positive"
assert X >= 1, "X has to be positive"
assert is_prime_power(q), "q has to be a prime power"

print(f"Given parameters: {K = }, {L = }, {X = }, {q = }")

# STEP 2: Compute parameters for construction
d = K*(L - 1) + 2*X - 1
g = (d - 1) // 2

print(f"Computed parameters: {d = }, {g = }")

# STEP 3: Compute pole numbers
f_pole_numbers = [2*j for j in range(X)] + [d + k for k in range(K)]
g_pole_numbers = [2*j for j in range(X)] + [2*X - 2 + K*(l + 1) for l in range(L)]
h_pole_numbers = sorted(set(f + g for f in f_pole_numbers for g in g_pole_numbers))

print("Pole numbers:")
print(f"{f_pole_numbers = }")
print(f"{g_pole_numbers = }")
print(f"{h_pole_numbers = }")

# STEP 4: Compute pole number table
print("Pole number table:")
print(pole_number_table(f_pole_numbers, g_pole_numbers))

# STEP 5: Create function field
Fq = GF(q)
F.<x> = FunctionField(Fq)

alphas = Fq.list()[:d]
s = prod(x - a for a in alphas)
r = 1 if q % 2 == 0 else 0  # add a nonzero r for fields of characteristic two

_.<T> = F[]
F.<y> = F.extension(T^2 + r*T - s)
x = F(x)

print("Function field genus:", F.genus())

assert 2 * F.genus() + 1 == d

Pinf = F.places_infinite(1)[0]

assert x.divisor_of_poles().multiplicity(Pinf) == 2
assert y.divisor_of_poles().multiplicity(Pinf) == d

# STEP 6: Create basis functions
f_basis = [choose_function_with_pole(F, n) for n in f_pole_numbers]
g_basis = [choose_function_with_pole(F, n) for n in g_pole_numbers]
h_basis = [choose_function_with_pole(F, n) for n in h_pole_numbers]

print(f"{f_basis = }")
print(f"{g_basis = }")
print(f"{h_basis = }")

# STEP 7: Choose places
all_places = F.places()

all_places.remove(Pinf)
places = []
x_coords = []

for Q in all_places:
    x_coord = x.evaluate(Q)
    if x_coord not in x_coords:
        x_coords.append(x_coord)
        places.append(Q)
        if len(places) > 2*K*L + 4*X - 4:
            break

assert len(places) > 2*K*L + 4*X - 4, "not enough rational places"

print(f"Available places: {len(all_places)}")
print(f"Number of chosen places: {len(places)}")

# STEP 8: Create generator matrices
GA = generator_matrix(f_basis, places)
GB = generator_matrix(g_basis, places)
GC = generator_matrix(h_basis, places)

information_set = GC.pivots()
GA = GA[:, information_set]
GB = GB[:, information_set]
GC = GC[:, information_set]

# STEP 9: Compute recovery threshold
recovery_threshold = N = len(information_set)
print(f"Recovery threshold: {recovery_threshold}")
print(f"(Expected: {3 * (K // 2) * L + K // 2 + 3*X - 2})")
print(f"(A3S: {A3S(K, L, X)})")
print(f"(GASP_big: {GASP_big(K, L, X)})")
print(f"(GASP_small: {GASP_small(K, L, X)})")
print(f"(GASP: {GASP(K, L, X)})")

# STEP 10: Test that this works
block_size = 10

# Create the matrices to be multiplied 
A = random_matrix(Fq, K*block_size, block_size)
B = random_matrix(Fq, block_size, L*block_size)

# Split the matrices
A_split = vsplit(A, K)
B_split = hsplit(B, L)

# Create the random matrices
R = [random_matrix(Fq, block_size, block_size) for _ in range(X)]
S = [random_matrix(Fq, block_size, block_size) for _ in range(X)]

# Encode the matrices using the generator matrices
A_enc = vector_matrix_multiplication(R + A_split, GA)
B_enc = vector_matrix_multiplication(S + B_split, GB)

# Compute the products (done by the workers)
C_enc = [a * b for a, b in zip(A_enc, B_enc)]

# Decode
encoded_blocks = vector_matrix_multiplication(C_enc, GC^-1)[-K*L:]  # take the last K*L blocks
C = assemble_block_matrix_in_column_order(encoded_blocks, K, L)

if C == A*B:
    print("Got the correct result!")
else:
    assert False, "got the wrong result"
