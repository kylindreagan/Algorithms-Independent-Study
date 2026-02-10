#https://open.kattis.com/problems/superdoku
#Max Flow

import sys
from collections import deque

input = sys.stdin.readline

n, k = map(int, input().split())
grid = [list(map(int, input().split())) for _ in range(k)]

used_col = [[False]*(n+1) for _ in range(n)]

for r in range(k):
    seen = set()
    for c in range(n):
        x = grid[r][c]

        if x in seen or used_col[c][x]:
            print("no")
            sys.exit()

        seen.add(x)
        used_col[c][x] = True

for _ in range(n-k):
    grid.append([0]*n)

def hopcroft_karp(adj):
    n_left = len(adj)
    n_right = n

    pairU = [-1]*n_left
    pairV = [-1]*(n_right+1)
    dist = [0]*n_left
    INF = 10**9

    def bfs():
        q = deque()
        for u in range(n_left):
            if pairU[u] == -1:
                dist[u] = 0
                q.append(u)
            else:
                dist[u] = INF

        found = False

        while q:
            u = q.popleft()
            for v in adj[u]:
                pu = pairV[v]
                if pu != -1 and dist[pu] == INF:
                    dist[pu] = dist[u] + 1
                    q.append(pu)
                elif pu == -1:
                    found = True

        return found

    def dfs(u):
        for v in adj[u]:
            pu = pairV[v]
            if pu == -1 or (dist[pu] == dist[u]+1 and dfs(pu)):
                pairU[u] = v
                pairV[v] = u
                return True
        dist[u] = INF
        return False

    matching = 0
    while bfs():
        for u in range(n_left):
            if pairU[u] == -1 and dfs(u):
                matching += 1

    return matching, pairU


for r in range(k, n):

    adj = [[] for _ in range(n)]

    for c in range(n):
        for num in range(1, n+1):
            if not used_col[c][num]:
                adj[c].append(num)

    matching, pairU = hopcroft_karp(adj)

    if matching < n:
        print("no")
        sys.exit()

    for c in range(n):
        num = pairU[c]
        grid[r][c] = num
        used_col[c][num] = True

print("yes")
for row in grid:
    print(*row)