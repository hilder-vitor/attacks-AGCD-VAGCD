from sage.modules.free_module_integer import IntegerLattice as Lattice


def orthogonal_lattice(list_of_basis_vectors):
    n = len(list_of_basis_vectors)
    d = len(list_of_basis_vectors[0])
    _B = Matrix(ZZ, n, d, list_of_basis_vectors)
    prod_vec = product([vector(ZZ, bi).norm() for bi in list_of_basis_vectors]) # ||b1|| * ... * ||bn||
    c = ZZ(ceil(2**((d-1)/2 + (d-n)*(d-n-1)/4) * prod_vec))
    B_aux = block_matrix([[c * _B], [1]]) # collums of B_aux are the basis of the auxiliary lattice
    L_aux = Lattice(B_aux.transpose(), lll_reduce = True) # taking transpose because sage spans the lattice using the rows of the given matrix

    basis_orthogonal = [L_aux.reduced_basis[i][-d:] for i in range(d-n)] # This [-d:] is used to take only the last d entries of each vector

    return Lattice(basis_orthogonal, lll_reduce = False) # B_ort is already an LLL-reduced basis



def sample_qi(gamma, eta):
    return ZZ.random_element(0, 2^(gamma - eta))

def sample_ri(rho):
    return ZZ.random_element(0, 2^(rho))

def sample_vec_agcd(p, gamma, eta, rho, m):
    return vector(ZZ, [p*sample_qi(gamma, eta) + sample_ri(rho) for _ in range(m)])

def sample_K_invK(p, m):
    K = Matrix.random(Integers(p), m, m)
    invK = K^-1
    return K.lift(), invK.lift()

def sample_C_bar(p, K, gamma, eta, rho, m, t):
    Cbar = Matrix(ZZ, t, m)
    for i in range(t):
        vi = sample_vec_agcd(p, gamma, eta, rho, m) * K
        for j in range(m):
            Cbar[i, j] = vi[j]

    return Cbar

def get_R_Q(p, invK, Cbar):
    R = (Cbar * invK) % p
    QK = (Cbar -  R*K) / p
    Q = (QK*invK) % p
    return R, Q


gamma = 20
rho = 4
eta = 10
m = 5

# Increasing c improves the probability that the attack works. If c < 1, then the it does not work
c = 2
t = ZZ(ceil(c * gamma * m / (eta - rho)))

print("m = %d" % m)
print("t = %d" % t)


p = random_prime(2^eta,False, 2^(eta-1))
x0 = p * random_prime(2^(gamma - eta), False, 2^(gamma - eta -1))
Zx0 = ZZ.quo(x0)
K, invK = sample_K_invK(p, m)

R = Matrix(ZZ, t, m, [sample_ri(rho) for _ in range(t*m)])
Q = Matrix(ZZ, t, m, [sample_qi(gamma, eta) for _ in range(t*m)])
C = (p*Q + R) % x0
Cbar = (C * K) % x0
B = block_matrix(ZZ, 2, 2, [[x0, 0], [-Cbar[m:t, :]*(Matrix(Zx0, Cbar[0:m, :])^-1).lift(), 1]])
L_ort_x0 = Lattice(B)

assert(L_ort_x0.rank() == t)
assert(L_ort_x0.dimension() == t)

for bi in L_ort_x0.reduced_basis:
    assert((bi * C) % x0 == 0)

L_ort_x0.BKZ(block_size=40, prune=10, fp='rr', precision=100)


basis_L_ort = []
for i in range(t-m):
    bi = L_ort_x0.reduced_basis[i]
    basis_L_ort.append(bi)

L_ort = Lattice(basis_L_ort)
assert(L_ort.rank() == t-m)
print("Constructed lattice L_ort whose vectors are ortogonal to the noise R")


L_R = orthogonal_lattice([bi for bi in L_ort.reduced_basis])

assert(L_R.rank() == m)
for j in range(m):
    rj = vector(R[:, j].list())
#    print("r%d = %s" % (j, rj))
    assert(rj in L_R)
print("Lattice L_R contains R .... OK")


B = L_R.reduced_basis.transpose()
assert(B.nrows() == t)
assert(B.ncols() == m)

B1 = Matrix(ZZ, m, m, B[0:m, :])
B2 = Matrix(ZZ, m, m, B[m:2*m, :])

Cbar1 = Matrix(ZZ, m, m, Cbar[0:m, :])
Cbar2 = Matrix(ZZ, m, m, Cbar[m:2*m, :])

W = B1 * B2^-1
M = Cbar1 * Cbar2^-1

###### just for debugging
R1 = Matrix(ZZ, m, m, R[0:m, :])
R2 = Matrix(ZZ, m, m, R[m:2*m, :])
assert(W == R1 * R2^-1)
###########################################

f = W.characteristic_polynomial()

print(f)
# all the entries of A are multiples of p
A = Matrix(ZZ, m, m, [a_b.numerator() for a_b in f(M).list()])

rec_p = A[0, 0]
for i in range(m):
    for j in range(m):
        rec_p = gcd(rec_p, A[i, j])
        if log(rec_p, 2) < eta:
            print("recovered p: %d" % rec_p)
            print("   secret p: %d" %p)
            exit(0)

