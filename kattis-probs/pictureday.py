#https://open.kattis.com/problems/pictureday
#MAX FLOW

import sys
from collections import deque

INF = 10**18

class MaxFlow:
    def __init__(self, n):
        self.n = n
        self.adj = [[] for _ in range(n)]
    
    def add_edge(self, fr, to, cap):
        forward = [to, cap, None]
        backward = [fr, 0, forward]
        forward[2] = backward
        self.adj[fr].append(forward)
        self.adj[to].append(backward)
    
    def bfs_level(self, s, t):
        level = [-1]*self.n
        q = deque([s])
        level[s] = 0
        while q:
            v = q.popleft()
            for to, cap, rev in self.adj[v]:
                if cap > 0 and level[to] < 0:
                    level[to] = level[v] + 1
                    q.append(to)
        return level
    
    def dfs(self, v, t, f, level, it):
        if v == t:
            return f
        for i in range(it[v], len(self.adj[v])):
            it[v] = i
            to, cap, rev = self.adj[v][i]
            if cap > 0 and level[v] < level[to]:
                ret = self.dfs(to, t, min(f, cap), level, it)
                if ret > 0:
                    self.adj[v][i][1] -= ret
                    rev[1] += ret
                    return ret
        return 0
    
    def max_flow(self, s, t):
        flow = 0
        while True:
            level = self.bfs_level(s, t)
            if level[t] < 0:
                break
            it = [0]*self.n
            while True:
                f = self.dfs(s, t, INF, level, it)
                if f == 0:
                    break
                flow += f
        return flow

def main():
    input_data = sys.stdin.read().strip().split()
    idx = 0
    N = int(input_data[idx]); idx += 1
    G_a = int(input_data[idx]); idx += 1
    G_b = int(input_data[idx]); idx += 1
    a = int(input_data[idx]); idx += 1
    b = int(input_data[idx]); idx += 1
    c = int(input_data[idx]); idx += 1
    
    A_groups = [[] for _ in range(G_a)]
    B_groups = [[] for _ in range(G_b)]
    
    student_to_A = [-1]*(N+1)
    student_to_B = [-1]*(N+1)
    
    for i in range(G_a):
        sz = int(input_data[idx]); idx += 1
        for _ in range(sz):
            s = int(input_data[idx]); idx += 1
            student_to_A[s] = i
            A_groups[i].append(s)
    
    for i in range(G_b):
        sz = int(input_data[idx]); idx += 1
        for _ in range(sz):
            s = int(input_data[idx]); idx += 1
            student_to_B[s] = i
            B_groups[i].append(s)
    
    count = [[0]*G_b for _ in range(G_a)]
    for s in range(1, N+1):
        ai = student_to_A[s]
        bi = student_to_B[s]
        count[ai][bi] += 1
    
    nodes = 2 + G_a + G_b
    source = 0
    sink = nodes - 1
    mf = MaxFlow(nodes)
    
    left_start = 1
    right_start = 1 + G_a
    
    for i in range(G_a):
        mf.add_edge(source, left_start + i, a)
    
    for j in range(G_b):
        mf.add_edge(right_start + j, sink, b)
    
    for i in range(G_a):
        for j in range(G_b):
            if count[i][j] > 0:
                mf.add_edge(left_start + i, right_start + j, c * count[i][j])
    
    result = mf.max_flow(source, sink)
    print(result)

if __name__ == "__main__":
    main()