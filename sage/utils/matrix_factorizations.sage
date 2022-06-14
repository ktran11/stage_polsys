## TODO : try Crout-based PLUQ (Algo 1 in DPS15) ??

#######################################
#  LSP factorization and left kernel  #
#           (naive storage)           #
#######################################

def LSP(A):
    # input: matrix A over a field.
    # returns L,S,P,rrp,nrrp such that
    # - A = L S P is the LSP decomposition of A
    # - rrp gives the indices of nonzero rows in S (in increasing order), which
    # is also the row rank profile of A
    # - nrrp gives the indices of zero rows in S (in increasing order)
    m,n = A.dimensions()
    rrp = []
    nrrp = []
    pivIndex = 0
    L = Matrix.identity(A.base_ring(),m,m)
    S = copy(A)
    P = [i for i in range(n)]
    for i in range(m):
        #find column with pivot element on row i, if there is some
        pivot = pivIndex
        while (pivot < n and S[i,pivot] == 0):
            pivot += 1
        if pivot == n:
            nrrp.append(i)
        else:
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
    return L,S,P,rrp,nrrp

def LSP_to_left_kernel(A):
    # input: matrix A over a field
    # returns: matrix K, basis of left kernel of A in reduced row echelon form
    # (with pivots = rightmost nonzero entries)
    m,n = A.dimensions()
    L,S,P,rrp,nrrp = LSP(A)
    Li = L.inverse()
    K = Matrix(A.base_ring(), len(nrrp), m)
    for i in range(len(nrrp)):
        K.set_row(i, Li[nrrp[i]])
    return K

def LiSP(A):
    # input: matrix A over a field.
    # returns L,S,P,rrp,nrrp such that
    # - A = L^{-1} S P is the LSP decomposition of A
    # - rrp gives the indices of nonzero rows in S (in increasing order), which
    # is also the row rank profile of A
    # - nrrp gives the indices of zero rows in S (in increasing order)
    m = A.nrows()
    n = A.ncols()
    rrp = []
    nrrp = []
    pivIndex = 0
    L = Matrix.identity(A.base_ring(),m,m)
    S = copy(A)
    P = [i for i in range(n)]
    for i in range(m):
        # find column, if any, with pivot element on row i
        pivot = pivIndex
        while (pivot < n and S[i,pivot] == 0):
            pivot += 1
        if pivot == n:
            nrrp.append(i)
        else:
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
    return L,S,P,rrp,nrrp

########################################
#  PLUQ factorization and left kernel  #
#          (compact storage)           #
########################################

def PLUQ(A):
    ## input: matrix A over a field, dimensions m x n.
    # returns LU,P,Q,rank such that 
    #  - LU is L and U compactly stored the usual way
    #  - A = P L U Q is the PLUQ decomposition of A
    #  - with P a permutation (list in [0...m])
    #  - and Q a permutation (list in [0...n])
    #  - rank is the rank of A
    # Rotations are performed to preserve the row rank profile, so that in the
    # end we have P = rrp+nrrp, where rrp is the row rank profile of A and nrrp
    # is the list of rows not in rrp (both rrp and nrrp are increasing)
    ## Using lexico order + row rotations + column rotations
    # ==> preserves rank profile matrix
    # ==> row transposition enough to preserve RPM TODO
    # ==> version with column transposition is enough to compute kernel..? go
    # for that one? what about row transposition? TODO
    m,n = A.dimensions()
    LU = copy(A)
    Q = [i for i in range(n)]
    P = [i for i in range(m)]
    nullity = 0
    rank = 0
    while (rank + nullity < m):
        #find column with pivot element on row rank, if there is some
        pivot = rank
        while (pivot < n and LU[rank,pivot] == 0):
            pivot += 1
        if pivot == n:
            #print(f"---rank{rank}--nullity---\n{LU}\n")
            P = row_rotation(LU,rank,P)
            #P = row_transposition(LU,rank,m-1-nullity,P)
            #print(f"---rank{rank}--nullity---\n{LU}\n")
            nullity += 1
        else:
            #print(f"---rank{rank}--rank---\n{LU}\n")
            #Q = column_rotation(LU,rank,pivot+1,Q)
            Q = column_transposition(LU,rank,pivot,Q)
            for k in range(rank+1,m):
                LU[k,rank] = LU[k,rank]/LU[rank,rank]
                for j in range(rank+1,n):
                    LU[k,j] = LU[k,j] - LU[k,rank]*LU[rank,j]
            #print(f"---rank{rank}--rank---\n{LU}\n")
            rank += 1
    return LU,P,Q,rank

def PLUQ_to_left_kernel(A):
    # input: matrix A over a field
    # returns: matrix K, basis of left kernel of A in reduced row echelon form
    # (with pivots = rightmost nonzero entries)
    m,n = A.dimensions()
    LU,P,Q,rank = PLUQ(A)
    for i in range(rank):
        LU[i,i] = 1
        for j in range(i+1,n):
            LU[i,j] = 0
    K = - LU[rank:,:rank] * LU[:rank,:rank].inverse()
    return K,P[:rank],P[rank:]


###########
#  Utils  #
###########

def column_rotation(A,cstart,cend,C):
    # columns 0, ..., cstart, .., cend-1, ..., n-1
    # are permuted into
    # columns 0, ..., cstart-1, cstart+1, .., cend-1, cstart, cend, ..., n-1
    # If a list C (of length n) is given, it is also
    # permuted accordingly
    n = A.ncols()
    perm_list = list(range(cstart)) + list(range(cstart+1,cend)) + [cstart] + list(range(cend,n))
    perm = Permutation([i+1 for i in perm_list]).inverse()
    A.permute_columns(perm)
    return perm.action(C)

def column_transposition(A,src,tgt,C):
    # columns 0, ..., src, ..., tgt, ..., n-1
    # are permuted into
    # columns 0, ..., tgt, src+1 ..., tgt-1, src, tgt+1, ... n-1
    # If a list R (of length m) is given, it is also
    # permuted accordingly
    A.swap_columns(src,tgt)
    tmp = C[src]
    C[src] = C[tgt]
    C[tgt] = tmp
    return C

def row_rotation(A,row,R):
    # rows 0, ..., row, ..., m-1
    # are permuted into
    # rows 0, ..., row-1, row+1, ..., m-1, row
    # If a list R (of length m) is given, it is also
    # permuted accordingly
    m = A.nrows()
    perm = Permutation([i+1 for i in range(m) if i != row] + [row+1])
    A.permute_rows(perm)
    return perm.action(R)

def row_transposition(A,src,tgt,R):
    # rows 0, ..., src, ..., tgt, ..., m-1
    # are permuted into
    # rows 0, ..., tgt, src+1 ..., tgt-1, src, tgt+1, ... m-1
    # If a list R (of length m) is given, it is also
    # permuted accordingly
    A.swap_rows(src,tgt)
    tmp = R[src]
    R[src] = R[tgt]
    R[tgt] = tmp
    return R

def PLUQ_to_LSP(LU,P,Q,rank):
    ## WARNING: this only works if row rotations were used (and column rotations?)
    m,n = LU.dimensions()
    rrp = P[:rank]
    nrrp = P[rank:]
    L,S = expand_PLUQ(LU,P,Q,rank)
    L.permute_rows(Permutation([i+1 for i in P]).inverse())
    L.permute_columns(Permutation([i+1 for i in P]).inverse())
    S.permute_rows(Permutation([i+1 for i in P]).inverse())
    #Q = Matrix(Permutation([i+1 for i in Q])).transpose()
    return L,S,Q,rrp,nrrp

def expand_PLUQ(LU,P,Q,rank):
    m,n = LU.dimensions()
    # retrieve L
    L = Matrix(LU.base_ring(), m, m)
    for j in range(rank):
        for i in range(j+1,m):
            L[i,j] = LU[i,j]
    for i in range(m):
        L[i,i] = 1
    # retrieve U
    U = Matrix(LU.base_ring(), m, n)
    U[:rank,:] = LU[:rank,:]
    for i in range(rank):
        for j in range(i):
            U[i,j] = 0
    return L,U

def expand_PLUQ_kernel(K,rrp,nrrp):
    nullity,rank = K.dimensions()
    K = K.augment(Matrix.identity(K.base_ring(), nullity, nullity))
    K.permute_columns(Permutation([i+1 for i in rrp+nrrp]).inverse())
    return K


#############
#  Testing  #
#############

def check_triL(L):
    m,n = L.dimensions()
    if m != n:
        return False
    unit = all([L[i,i] == 1 for i in range(m)]) 
    if not unit:
        return False
    tri = all([L[i,j] == 0 for i in range(m) for j in range(i+1,m)])
    if not tri:
        return False
    return True

def check_triU(U,rank):
    m,n = U.dimensions()
    if rank > m or rank > n:
        return False
    inv = all([U[i,i] != 0 for i in range(rank)]) 
    if not inv:
        return False
    tri = all([U[i,j] == 0 for i in range(m) for j in range(min(i,n))])
    if not tri:
        return False
    return True

def check_many_PLUQ(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        LU,P,Q,rank = PLUQ(A)
        L,U = expand_PLUQ(LU,P,Q,rank)
        correct = rank==A.rank()             \
                      and check_triL(L)      \
                      and check_triU(U,rank)
        L.permute_rows(Permutation([i+1 for i in P]).inverse())
        #L.permute_columns(Permutation([i+1 for i in P]).inverse())
        #U.permute_rows(Permutation([i+1 for i in P]).inverse())
        U.permute_columns(Permutation([i+1 for i in Q]).inverse())
        correct = correct and L*U == A
        i += 1
    if correct:
        print("tall rectangular: ok")
    else:
        print("tall rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,3,6)
        LU,P,Q,rank = PLUQ(A)
        L,U = expand_PLUQ(LU,P,Q,rank)
        correct = rank==A.rank()             \
                      and check_triL(L)      \
                      and check_triU(U,rank)
        L.permute_rows(Permutation([i+1 for i in P]).inverse())
        #L.permute_columns(Permutation([i+1 for i in P]).inverse())
        #U.permute_rows(Permutation([i+1 for i in P]).inverse())
        U.permute_columns(Permutation([i+1 for i in Q]).inverse())
        correct = correct and L*U == A
        i += 1
    if correct:
        print("wide rectangular: ok")
    else:
        print("wide rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,4,4)
        LU,P,Q,rank = PLUQ(A)
        L,U = expand_PLUQ(LU,P,Q,rank)
        correct = rank==A.rank()             \
                      and check_triL(L)      \
                      and check_triU(U,rank)
        L.permute_rows(Permutation([i+1 for i in P]).inverse())
        #L.permute_columns(Permutation([i+1 for i in P]).inverse())
        #U.permute_rows(Permutation([i+1 for i in P]).inverse())
        U.permute_columns(Permutation([i+1 for i in Q]).inverse())
        correct = correct and L*U == A
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A

    return True


def check_many_LSP(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        Li1,S1,P1,rrp1,nrrp1 = LiSP(A)
        L2,S2,P2,rrp2,nrrp2 = LSP(A)
        P2mat = Matrix(Permutation([i+1 for i in P2])).transpose() # transform list into column permutation
        correct = correct and (Li1.inverse() == L2 and (S1,P1,rrp1,nrrp1) == (S2,P2,rrp2,nrrp2) and L2*S2*P2mat == A)
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
        Li1,S1,P1,rrp1,nrrp1 = LiSP(A)
        L2,S2,P2,rrp2,nrrp2 = LSP(A)
        P2mat = Matrix(Permutation([i+1 for i in P2])).transpose() # transform list into column permutation
        correct = correct and (Li1.inverse() == L2 and (S1,P1,rrp1,nrrp1) == (S2,P2,rrp2,nrrp2) and L2*S2*P2mat == A)
        i += 1
    if correct:
        print("wide rectangular: ok")
    else:
        print("wide rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,4,4)
        Li1,S1,P1,rrp1,nrrp1 = LiSP(A)
        L2,S2,P2,rrp2,nrrp2 = LSP(A)
        P2mat = Matrix(Permutation([i+1 for i in P2])).transpose() # transform list into column permutation
        correct = correct and (Li1.inverse() == L2 and (S1,P1,rrp1,nrrp1) == (S2,P2,rrp2,nrrp2) and L2*S2*P2mat == A)
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A

    return True

def check_many_PLUQ_to_kernel(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        K = expand_PLUQ_kernel(*PLUQ_to_left_kernel(A))
        correct == (A.left_kernel() == (K.row_space()))
        i += 1
    if correct:
        print("tall rectangular: ok")
    else:
        print("tall rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,3,6)
        K = expand_PLUQ_kernel(*PLUQ_to_left_kernel(A))
        correct == (A.left_kernel() == (K.row_space()))
        i += 1
    if correct:
        print("wide rectangular: ok")
    else:
        print("wide rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,4,4)
        K = expand_PLUQ_kernel(*PLUQ_to_left_kernel(A))
        correct == (A.left_kernel() == (K.row_space()))
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A

def check_many_PLUQker_is_LSPker(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        K1 = expand_PLUQ_kernel(*PLUQ_to_left_kernel(A))
        K2 = LSP_to_left_kernel(A)
        correct == (K1 == K2)
        i += 1
    if correct:
        print("tall rectangular: ok")
    else:
        print("tall rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,3,6)
        K1 = expand_PLUQ_kernel(*PLUQ_to_left_kernel(A))
        K2 = LSP_to_left_kernel(A)
        correct == (K1 == K2)
        i += 1
    if correct:
        print("wide rectangular: ok")
    else:
        print("wide rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,4,4)
        K1 = expand_PLUQ_kernel(*PLUQ_to_left_kernel(A))
        K2 = LSP_to_left_kernel(A)
        correct == (K1 == K2)
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A


def check_many_LSP_to_kernel(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        correct == (A.left_kernel() == (LSP_to_left_kernel(A).row_space()))
        i += 1
    if correct:
        print("tall rectangular: ok")
    else:
        print("tall rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,3,6)
        correct == (A.left_kernel() == (LSP_to_left_kernel(A).row_space()))
        i += 1
    if correct:
        print("wide rectangular: ok")
    else:
        print("wide rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,4,4)
        correct == (A.left_kernel() == (LSP_to_left_kernel(A).row_space()))
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A

    return True

# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
