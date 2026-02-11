#https://open.kattis.com/problems/goatropes
#SIMPLEX (Linear Programming)

import sys
import math
class Simplex:
    def __init__(self, c, A, b):
        # c: objective coefficients (maximize c·x)
        # A: constraint matrix (Ax <= b)
        # b: constraint bounds
        self.m = len(A)  # constraints
        self.n = len(c)  # variables
        self.tab = [[0.0] * (self.n + self.m + 1) for _ in range(self.m + 1)]
        
        # Fill objective row
        for j in range(self.n):
            self.tab[0][j] = -c[j]  # simplex minimizes, so negate for max
        
        # Fill constraint rows
        for i in range(self.m):
            for j in range(self.n):
                self.tab[i+1][j] = A[i][j]
            self.tab[i+1][self.n + i] = 1.0  # slack variable
            self.tab[i+1][-1] = b[i]
        
        # Basis indices for slacks
        self.basis = list(range(self.n, self.n + self.m))
    
    def pivot(self, row, col):
        m, n = self.m, self.n + self.m
        pivot_val = self.tab[row][col]
        
        # Normalize pivot row
        for j in range(n + 1):
            self.tab[row][j] /= pivot_val
        
        # Eliminate column from other rows
        for i in range(self.m + 1):
            if i != row:
                factor = self.tab[i][col]
                if abs(factor) > 1e-10:
                    for j in range(n + 1):
                        self.tab[i][j] -= factor * self.tab[row][j]
        
        # Update basis
        self.basis[row - 1] = col
    
    def solve(self) -> float:
        m, n = self.m, self.n + self.m
        
        while True:
            # Find entering column (most negative in objective row)
            col = -1
            min_val = 0.0
            for j in range(n):
                if self.tab[0][j] < min_val - 1e-10:
                    min_val = self.tab[0][j]
                    col = j
            
            if col == -1:  # optimal
                break
            
            # Find leaving row (minimum ratio)
            row = -1
            min_ratio = float('inf')
            for i in range(1, m + 1):
                if self.tab[i][col] > 1e-10:
                    ratio = self.tab[i][-1] / self.tab[i][col]
                    if ratio < min_ratio - 1e-10:
                        min_ratio = ratio
                        row = i
            
            if row == -1:  # unbounded
                return float('inf')
            
            self.pivot(row, col)
        
        return self.tab[0][-1]  # optimal value (negated back later)

def solve_goat_ropes():
    data = sys.stdin.read().strip().split()
    if not data:
        return
    
    idx = 0
    n = int(data[idx]); idx += 1
    coords = []
    for _ in range(n):
        x = int(data[idx]); idx += 1
        y = int(data[idx]); idx += 1
        coords.append((x, y))
    
    # Compute distances
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        xi, yi = coords[i]
        for j in range(i + 1, n):
            xj, yj = coords[j]
            d = math.hypot(xi - xj, yi - yj)
            dist[i][j] = d
            dist[j][i] = d
    
    # LP: maximize sum(r_i)
    # Constraints: r_i + r_j <= dist[i][j] for i<j
    c = [1.0] * n  # objective coefficients
    
    # Build constraints matrix
    constraints = []
    b = []
    for i in range(n):
        for j in range(i + 1, n):
            row = [0.0] * n
            row[i] = 1.0
            row[j] = 1.0
            constraints.append(row)
            b.append(dist[i][j])
    
    simplex = Simplex(c, constraints, b)
    result = simplex.solve()
    
    print(f"{result:.2f}")

if __name__ == "__main__":
    solve_goat_ropes()