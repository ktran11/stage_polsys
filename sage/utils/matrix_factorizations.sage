def LiSP(A):
    # input: matrix A over a field.
    # returns L,S,P,nonzRows such that A = L^-1 S P is
    #the LSP decomposition of A. nonzRows gives the
    #indices of nonzero rows in S.
    m = A.nrows()
    n = A.ncols()
    nonzRows = []
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
            nonzRows.append(i)
            S.swap_columns(pivIndex,pivot)
            P.swap_rows(pivIndex,pivot)
            for k in range(i+1,m):
                cstFactor = -S[k,pivIndex]/S[i,pivIndex]
                S.add_multiple_of_row(k,i,cstFactor)
                L.add_multiple_of_row(k,i,cstFactor)
            pivIndex += 1
    return L,S,P,nonzRows

def LiSP_compact(A):
    # input: matrix A over a field.
    # returns L,S,P,nonzRows such that A = L^-1 S P is
    #the LSP decomposition of A. nonzRows gives the
    #indices of nonzero rows in S.
    # Representation is compact: L and S are merged
    # TODO work in progress!
    m = A.nrows()
    n = A.ncols()
    nonzRows = []
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
            nonzRows.append(i)
            S.swap_columns(pivIndex,pivot)
            # simulate P.swap_columns(pivIndex,pivot)
            tmp = P[pivIndex]
            P[pivIndex] = P[pivot]
            P[pivot] = tmp
            for k in range(i+1,m):
                cstFactor = -S[k,pivIndex]/S[i,pivIndex]
                S.add_multiple_of_row(k,i,cstFactor)
                L.add_multiple_of_row(k,i,cstFactor)
            pivIndex += 1
    return (L,S,P,nonzRows)

def expand_LSP(L,S,P,nz):
    return L, S, Permutation(P).inverse().to_matrix(), nz

def check_many(max_iter=1000):
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(GF(2),6,3)
        correct = (expand_LSP(*LSP_compact(A)) == LSP(A))
        i += 1
    print("tall rectangular: ok")

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(GF(2),3,6)
        correct = (expand_LSP(*LSP_compact(A)) == LSP(A))
        i += 1
    print("wide rectangular: ok")

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(GF(2),4,4)
        correct = (expand_LSP(*LSP_compact(A)) == LSP(A))
        i += 1
    print("square: ok")
    return correct,A

# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
