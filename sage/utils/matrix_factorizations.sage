def LSP(A):
    # input: matrix A over a field.
    # returns L,S,P,nzRows such that A = L S P is
    #the LSP decomposition of A. nzRows gives the
    #indices of nonzero rows in S.
    m = A.nrows()
    n = A.ncols()
    nzRows = []
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
            nzRows.append(i)
            S.swap_columns(pivIndex,pivot)
            # simulate P.swap_rows(pivIndex,pivot)
            tmp = P[pivIndex]
            P[pivIndex] = P[pivot]
            P[pivot] = tmp
            for k in range(i+1,m):
                L[k,i] = S[k,pivIndex]/S[i,pivIndex]
                S.add_multiple_of_row(k,i,-L[k,i])
            pivIndex += 1
    return L,S,P,nzRows



def LiSP(A):
    # input: matrix A over a field.
    # returns L,S,P,nzRows such that A = L^-1 S P is
    #the LSP decomposition of A. nzRows gives the
    #indices of nonzero rows in S.
    m = A.nrows()
    n = A.ncols()
    nzRows = []
    pivIndex = 0
    L = Matrix.identity(A.base_ring(),m,m)
    S = copy(A)
    P = Matrix.identity(A.base_ring(),n,n)
    for i in range(m):
        #find column with pivot element on row i, if there is some
        pivot = pivIndex
        while (pivot < n and S[i,pivot] == 0):
            pivot += 1
        if pivot < n:
            nzRows.append(i)
            S.swap_columns(pivIndex,pivot)
            P.swap_rows(pivIndex,pivot)
            for k in range(i+1,m):
                cstFactor = -S[k,pivIndex]/S[i,pivIndex]
                S.add_multiple_of_row(k,i,cstFactor)
                L.add_multiple_of_row(k,i,cstFactor)
            pivIndex += 1
    return L,S,P,nzRows

def LiSP_compact(A):
    # input: matrix A over a field.
    # returns L,S,P,nzRows such that A = L^-1 S P is
    #the LSP decomposition of A. nzRows gives the
    #indices of nonzero rows in S.
    # Representation is compact: L and S are merged
    # TODO work in progress!
    m = A.nrows()
    n = A.ncols()
    nzRows = []
    pivIndex = 0
    L = Matrix.identity(A.base_ring(),m,m)
    S = copy(A)
    P = [i for i in range(1,n+1)]
    for i in range(m):
        # find column, if any, with pivot element on row i
        pivot = pivIndex
        while (pivot < n and S[i,pivot] == 0):
            pivot += 1
        if pivot < n:
            nzRows.append(i)
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
    return (L,S,P,nzRows)

def expand_LSP(L,S,P,nz):
    return L, S, Permutation(P).inverse().to_matrix(), nz

def check_many(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        Li1,S1,P1,nz1 = LiSP(A)
        L2,S2,P2,nz2 = LSP(A)
        P2 = Matrix(Permutation([i+1 for i in P2])).transpose() # transform list into column permutation
        correct = correct and (Li1.inverse() == L2 and (S1,P1,nz1) == (S2,P2,nz2) and L2*S2*P2 == A)
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
        P2 = Matrix(Permutation([i+1 for i in P2])).transpose() # transform list into column permutation
        correct = correct and (Li1.inverse() == L2 and (S1,P1,nz1) == (S2,P2,nz2) and L2*S2*P2 == A)
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
        P2 = Matrix(Permutation([i+1 for i in P2])).transpose() # transform list into column permutation
        correct = correct and (Li1.inverse() == L2 and (S1,P1,nz1) == (S2,P2,nz2) and L2*S2*P2 == A)
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A

    return True

# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
