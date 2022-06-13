def LSP(A):
    # input: matrix A over a field.
    # returns L,S,P,rrp such that A = L S P is the LSP decomposition of A. rrp
    # gives the indices of nonzero rows in S (in increasing order), which is
    # also the row rank profile of A.
    m,n = A.dimensions()
    rrp = []
    pivIndex = 0
    L = Matrix.identity(A.base_ring(),m,m)
    S = copy(A)
    P = [i for i in range(n)]
    for i in range(m):
        #find column with pivot element on row i, if there is some
        pivot = pivIndex
        while (pivot < n and S[i,pivot] == 0):
            pivot += 1
        if pivot < n:
            rrp.append(i)
            S.swap_columns(pivIndex,pivot)
            # simulate P.swap_rows(pivIndex,pivot)
            tmp = P[pivIndex]
            P[pivIndex] = P[pivot]
            P[pivot] = tmp
            for k in range(i+1,m):
                L[k,i] = S[k,pivIndex]/S[i,pivIndex]
                S.add_multiple_of_row(k,i,-L[k,i])
            pivIndex += 1
    return L,S,P,rrp

def left_kernel(A):
    # input: matrix A over a field
    # returns: matrix K, basis of left kernel of A in reduced row echelon form
    # (with pivots = rightmost nonzero entries)
    m,n = A.dimensions()
    L,S,P,rrp = LSP(A)
    Li = L.inverse()
    rk = len(rrp)
    K = Matrix(A.base_ring(), m-rk, m)
    krow = 0
    ind = 0
    i = 0
    while krow < rk:
        if i != rrp[krow]:
            K.set_row(ind, Li[i])
            ind += 1
        else:
            krow += 1
        i += 1
    for ii in range(i,m):
        K.set_row(ii-i+ind, Li[ii])
    return K



## TODO : make refined version with more compact storage of L&S and of kernel
## TODO : deduce nullspace from LSP
## TODO : find out how this differs from Crout-based PLUQ (Algo 1 in DPS15)

def LiSP(A):
    # input: matrix A over a field.
    # returns L,S,P,rrp such that A = L^-1 S P is
    #the LSP decomposition of A. rrp gives the
    #indices of nonzero rows in S, which is also
    # the row rank profile of A.
    m = A.nrows()
    n = A.ncols()
    rrp = []
    pivIndex = 0
    L = Matrix.identity(A.base_ring(),m,m)
    S = copy(A)
    P = [i for i in range(n)]
    for i in range(m):
        # find column, if any, with pivot element on row i
        pivot = pivIndex
        while (pivot < n and S[i,pivot] == 0):
            pivot += 1
        if pivot < n:
            rrp.append(i)
            S.swap_columns(pivIndex,pivot)
            # simulate P.swap_rows(pivIndex,pivot)
            tmp = P[pivIndex]
            P[pivIndex] = P[pivot]
            P[pivot] = tmp
            for k in range(i+1,m):
                cstFactor = -S[k,pivIndex]/S[i,pivIndex]
                S.add_multiple_of_row(k,i,cstFactor)
                L.add_multiple_of_row(k,i,cstFactor)
            pivIndex += 1
    return L,S,P,rrp


def check_many_LSP(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        Li1,S1,P1,nz1 = LiSP(A)
        L2,S2,P2,nz2 = LSP(A)
        P2mat = Matrix(Permutation([i+1 for i in P2])).transpose() # transform list into column permutation
        correct = correct and (Li1.inverse() == L2 and (S1,P1,nz1) == (S2,P2,nz2) and L2*S2*P2mat == A)
        #correct = (expand_LSP(*LSP_compact(A)) == LSP(A))
        i += 1
    if correct:
        print("tall rectangular: ok")
    else:
        print("tall rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,3,6)
        Li1,S1,P1,nz1 = LiSP(A)
        L2,S2,P2,nz2 = LSP(A)
        P2mat = Matrix(Permutation([i+1 for i in P2])).transpose() # transform list into column permutation
        correct = correct and (Li1.inverse() == L2 and (S1,P1,nz1) == (S2,P2,nz2) and L2*S2*P2mat == A)
        i += 1
    if correct:
        print("wide rectangular: ok")
    else:
        print("wide rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,4,4)
        Li1,S1,P1,nz1 = LiSP(A)
        L2,S2,P2,nz2 = LSP(A)
        P2mat = Matrix(Permutation([i+1 for i in P2])).transpose() # transform list into column permutation
        correct = correct and (Li1.inverse() == L2 and (S1,P1,nz1) == (S2,P2,nz2) and L2*S2*P2mat == A)
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A

    return True

def check_many_kernel(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        correct == (A.left_kernel() == (left_kernel(A).row_space()))
        i += 1
    if correct:
        print("tall rectangular: ok")
    else:
        print("tall rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,3,6)
        correct == (A.left_kernel() == (left_kernel(A).row_space()))
        i += 1
    if correct:
        print("wide rectangular: ok")
    else:
        print("wide rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,4,4)
        correct == (A.left_kernel() == (left_kernel(A).row_space()))
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A

    return True

# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
