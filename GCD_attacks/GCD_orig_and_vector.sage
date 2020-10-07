#################################################################################
#       GCD attacks on the original AGCD and on the VAGCD problems.
#
#   These attacks output an integer g that is equal to p with high probability.
#   If g is not equal to p, then g should be a multiple of p small enough for
# us to be able to factor g and recover p by taking the eta-bit factor.
#
#################################################################################


def genXi(rho,eta,gam,p):
  return p*ZZ.random_element(2^(gam-eta))+ZZ.random_element(2^rho)

def genVecXi(rho,eta,gam,m,p,K):
  v=vector([genXi(rho,eta,gam,p) for i in range(m)])
  return K*v 

def genK(p, m):
  K = Matrix.random(Integers(p), m, m)
  return K.lift()

# Multipoint polynomial evaluation.
# This is the simplest possible recursive implementation.
# Efficiency could be improved (probably by a constant factor) by
# precomputing and storing the product tree corresponding to the list li
def multiPointEval(f,li,R):
  n,x = len(li), R.0
  if n==1:
      return [f(x=li[0])]
  li1, li2 = li[:n//2], li[n//2:]
  f1 = f.quo_rem(prod((x-xi for xi in li1)))[1]
  f2 = f.quo_rem(prod((x-xi for xi in li2)))[1]
  return multiPointEval(f1, li1, R) + multiPointEval(f2, li2, R)

def attackGCDMultiPoint(rho=12,eta=100,gam=1000,verbose=True):
    p=random_prime(2^eta,lbound=2^(eta-1),proof=False)
    if verbose:
        print("rho = %d, eta = %d, gamma = %d, p = %d" % (rho, eta, gam, p))

    t = cputime(subprocesses=True)
    s = rho
    B = floor(2^(1. * rho * (s + 1) / (s - 1)))
    fac_B = factorial(B)

    # defining the polynomial ring R and the variable x
    R.<x> = ZZ['x']

    g = 1
    L1 = [genXi(rho, eta, gam, p) for _ in range(ceil(2^(rho/2+1)))]
    f = prod([(x - ci) for ci in L1])
    for j in range(1, s+1):
        print("iteration: %d / %d" % (j, s))
        L2 = [genXi(rho, eta, gam, p) for _ in range(ceil(2^(rho/2)))]
        eval_f = multiPointEval(f, L2, R)

        print("bitlen(f(c1)) = %d" % eval_f[0].nbits())
        yj = prod(eval_f)
          
        if not p.divides(yj):
            print("_____________ p DOES NOT divide product %d" % j)

        print("bitlen(yj) = %d" % yj.nbits())
        if g == 1:
            g = prime_to_m_part(yj, fac_B)
        else:
            print("bitlen(g) = %d" % g.nbits())
            _g = gcd(g, yj)
            # just update g if the new gcd isn't smaller than p
            # (which can happen if p doesn't divide some y)
            if _g.nbits() >= eta - 1:
        #        g = prime_to_m_part(_g, fac_B)
                g = _g
        print("bitlen(g) = %d" % g.nbits())
        print("does p divides g? %s\n" % (p.divides(g)))

        if eta - 1 <= g.nbits() <= eta:
            break
  
    print("\nFINISHED:\ntime: %s" % cputime(t))
    print("g == p? %s" % (g == p))
    print("does p divides g? %s" % (p.divides(g)))
    print("bitlen(g) = %d" % g.nbits())
    print("g.factor() = %s" % g.factor())


def attackGCDVectorMultiPoint(rho=3, eta=15,gam=45,m=2,verbose=True):
    p=random_prime(2^eta,lbound=2^(eta-1),proof=False)
    if verbose:
        print("m = %d, rho = %d, eta = %d, gamma = %d, p = %d" % (m, rho, eta, gam, p))

    t = cputime(subprocesses=True)
    s = rho
    B = floor(2^(1. * m * rho * (m * rho + 1) / (m * rho - 1)))
    fac_B = factorial(B)

    K = genK(p, m)


    # defining the polynomial ring R and the variable x
    R.<x> = ZZ['x']

    g = 1
    L1 = [genVecXi(rho, eta, gam, m, p, K) for _ in range(ceil(2^(m * rho/2 + 1)))]
    polys = [prod([(x - c[l]) for c in L1]) for l in range(m)]
    for j in range(1, s+1):
        print("iteration: %d / %d" % (j, s))
        L2 = [genVecXi(rho, eta, gam, m, p, K) for _ in range(ceil(2^(m * rho/2)))]

        for l in range(0, m):
            f = polys[l]
            eval_f = multiPointEval(f, [c[l] for c in L2], R)
            y = prod(eval_f)

            if g == 1:
                g = prime_to_m_part(y, fac_B)
            else:
                _g = gcd(g, y)
                # just update g if the new gcd isn't smaller than p
                # (which can happen if p doesn't divide some y)
                if _g.nbits() >= eta - 1:
                    g = prime_to_m_part(_g, fac_B)
                print("bitlen(g) = %d\n" % g.nbits())

                if eta - 1 <= g.nbits() <= eta:
                    break;
        if eta - 1 <= g.nbits() <= eta:
            break

    print("\nFINISHED:\ntime: %s" % cputime(t))
    print("g == p? %s" % (g == p))
    print("does p divides g? %s" % (p.divides(g)))
    print("bitlen(g) = %d" % g.nbits())
    print("g.factor() = %s" % g.factor())



def attackGCD_eurocrypt2012(rho=12, eta=100, gamma=1000, verbose=True):
    p = random_prime(2^eta,lbound=2^(eta-1),proof=False)
    if verbose:
        print("rho = %d, eta = %d, gamma = %d, p = %d" % (rho, eta, gamma, p))
    
    t = cputime(subprocesses=True)
    s = rho
 
    B = floor(2^(1. * rho * (s + 1) / (s - 1)))
    fac_B = factorial(B)
   
    for j in range(1, s+1):
        x = genXi(rho, eta, gamma, p)
        z = prod([x-i for i in range(2^rho)])
        if j == 1:
            g = z
            continue

        g = prime_to_m_part(gcd(g, z), fac_B)

        print("j = %d, gcd size = %d" % (j, g.nbits()))
        if eta - 1 <= g.nbits() <= eta:
            break
 
    print("\nFINISHED:\ntime: %s" % cputime(t))

    print("g == p? %s" % (g == p))
    print("does p divides g? %s" % (p.divides(g)))
    print("bitlen(g) = %d" % g.nbits())
    print("g.factor() = %s" % g.factor())


if __name__== "__main__":
    rho = 5
    m = 2
    eta = 50
    gamma = 500
    attackGCDMultiPoint(rho, eta, gamma)
    print("---------------------------------------------------------")
    attackGCD_eurocrypt2012(rho, eta, gamma)
    print("---------------------------------------------------------")
    attackGCDVectorMultiPoint(rho, eta, gamma, m)
