#################################################################################
#       Orthogonal lattice attack on the original AGCD
#   as described in the paper Algorithms for the approximate common
#   divisor problem, by Galbraith, Gebregiyorgis and Murphy.
#
#   No noiseless x0
#   Single prime
#
#       Notice that you must by yourself the block size of BKZ that achieves the 
#   chosen root-Hermite factor.
#################################################################################


from sage.modules.free_module_integer import IntegerLattice as Lattice


def sample_r(rho):
	return ZZ.random_element(-2^rho+1,2^rho)

def sample_vec_r(rho, m):
	return vector(ZZ, [sample_r(rho) for _ in range(m)])

def sample_q(gamma, eta):
	return ZZ.random_element(0, 2^(gamma - eta))

def sample_vec_q(gamma, eta, m):
	return vector(ZZ, [sample_q(gamma, eta) for _ in range(m)])


def sample_basis(gamma, eta, rho, t):
    R = 2^rho
    B = Matrix(ZZ, t, t+1)
    r = sample_vec_r(rho, t)
    q = sample_vec_q(gamma, eta, t)
    c = p*q + r
    for i in range(t):
        B[i, 0] = c[i]
        B[i, i+1] = R
    return r, q, c, B

rho = 47
eta = 50
rhfLLL = 1.014
gamma_discriminant_neg = ZZ(ceil((eta - rho)^2 / (4 * log(rhfLLL,2))));
gamma = ceil(gamma_discriminant_neg / 2)

print("gamma to rule out LLL: %d" % gamma_discriminant_neg)
print("chosen gamma: %d" % gamma)

t = ceil(sqrt( gamma / (4*log(rhfLLL,2)) ))

p = random_prime(2^eta, False, 2^(eta-1))

r, q, c, B = sample_basis(gamma, eta, rho, t) # r and q are used for debugging


L = Lattice(B)
L.BKZ(block_size=30, prune=15) # XXX: you must choose the correct block size for the given root-Hermite factor


total = 0
vecs_ort_q = []
avg_norm = 0
avg_rhf_LLL = 0
for i in range(L.reduced_basis.nrows()):
    vi = L.reduced_basis[i]
    ui = B.solve_left(vi)
    if 0 == ui * q:
        total += 1
        vecs_ort_q.append(ui)

print("Gap: eta - rho = %d" % (eta - rho))

print("")
print("t = %d" % t)
print("total of vectors orthogonal to q: %d" % total)

if total < t-1:
    print("ATK NOT OK")
    exit(1)

U = Matrix(ZZ, t-1, t, vecs_ort_q)
V = U.right_kernel()
rec_q = V.basis()[0]
rec_p = c[0] / rec_q[0]
print("original  p: %d" % p)
print("recovered p: %d" % rec_p)


if eta - 1 <= log(rec_p, 2) <= eta:
    print("ATK OK")
else:
    print("ATK NOT OK")
