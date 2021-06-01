import numpy as np 

def inversa(A):



def lu_decomp(A):
    
    #nao checa se esse metodo pode ser aplicado na matriz
    
    n = A.shape[0]
    u = np.zeros((n, n))
    l = np.identity(n)

    for i in range(n):
        for j in range(i, n):
            u[i, j] = A[i, j] - sum([l[i, k]*u[k, j] for k in range(i)]) 
        for j in range(i+1, n):
            l[j, i] = (A[j, i] - sum([l[j, k]*u[k, i] for k in range(i)]))/u[i, i]
    
    return l, u

def gaussian_elimination(A, b):
    
    # nao checa se o primeiro elemento da matriz pivo eh zero ou se o det da mesma eh != 0
    
    n = A.shape[0]
    augA = np.concatenate((A, b.T), axis=1)
    augA = augA.astype(float)
    
    for i in range(n):
        augA[i] = augA[i]/augA[i, i]
        for j in range (i+1, n):
            augA[j] = augA[j] + augA[i]*(-augA[j, i]/augA[i, i])
    
    for i in range(n):
        for j in range(i+1, n):
            augA[-(j+1)] = augA[-(j+1)] + augA[-(i+1)]*(-augA[-(j+1), -(i+2)]/augA[-(i+1), -(i+2)])

    return augA


def cholesky_decomp(A):
    #   não verifica se A é uma Matriz Simetrica Positiva Definida
    n = A.shape[0]
    H = np.zeros((n, n))

    for i in range(n):
        H[i, i] = np.sqrt(A[i, i] - sum([H[i, j]**2 for j in range(i)]))
        for j in range(i+1, n):
            H[j, i] = (A[i, j] - sum([ H[i, k]*H[j, k] for k in range(i)]))/ H[i, i]

    return H

def gauss_jacobi(A, b, x0, e):

    #precisa checar se converge

    D = np.diag(np.diag(A))
    C = np.identity(A.shape[0]) - (np.linalg.solve(D, A))
    g = np.linalg.solve(D, b)

    for i in range(5000):
        if(np.linalg.norm(A*x0 - b) <= e):
            break
        x0 = C@x0 + g

    return x0 
    

def gauss_seidel(A, b, x0, e):
    
    #precisa testar se a sequência converge (criterio de sassenfeld)

    L = np.tril(A)
    R = np.triu(A, 1)
    C = -np.linalg.solve(L, R)
    g = np.linalg.solve(L, b)

    for i in range(5000):
        if(np.linalg.norm(A*x0 - b) <= e):
            break
        x0 = C@x0 + g

    return x0 

def inverse(A):
    det = np.linalg.det(A)
    if det == 0:
        raise Exeption("Erro: Determinante da matrix == 0")
    
    C = np.zeros(A.shape)
    nrows, ncols = C.shape
    for i in range(nrows):
        for j in range(ncols):
            C[i, j] = (-1)**(i+j)*np.linalg.det(np.delete(np.delete(A, i, axis=0), j, axis=1))
    return C.T/det
    


