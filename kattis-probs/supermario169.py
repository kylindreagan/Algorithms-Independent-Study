#WORK IN PROGRESS
#https://open.kattis.com/problems/supermario169
#TSP

import math
import sys

INF = float("inf")

def dist_sq(a, b):
    return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2

def dist(a, b):
    return math.sqrt(dist_sq(a, b))

def solve_cluster_tsp(switch_pos, coins):
    k = len(coins)
    if k == 0:
        return 0.0, [0.0], [switch_pos]
    d0 = [dist(switch_pos, c) for c in coins]
    
    # coin to coin distances
    d = [[0.0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            d_ij = dist(coins[i], coins[j])
            d[i][j] = d_ij
            d[j][i] = d_ij
    
    dp = [[INF] * k for _ in range(1<<k)]
    
    for i in range(k):
        dp[1<<i][i] = d0[i]
    
    for mask in range(1, 1<<k):
        # Find all possible last coins
        for last in range(k):
            if not (mask & (1<<last)) or dp[mask][last] == INF:
                continue
            
            # Try to add a new coin
            for nxt in range(k):
                if mask & (1<<nxt):
                    continue
                new_mask = mask | (1<<nxt)
                new_cost = dp[mask][last] + d[last][nxt]
                if dp[new_mask][nxt] > new_cost:
                    dp[new_mask][nxt] = new_cost
    
    full_mask = (1<<k) - 1
    
    return_to_switch = [dist(coins[i], switch_pos) for i in range(k)]
    
    # Cost to collect all coins and return to switch
    min_return_cost = INF
    for i in range(k):
        total = dp[full_mask][i] + return_to_switch[i]
        if total < min_return_cost:
            min_return_cost = total
    
    end_costs = [dp[full_mask][i] for i in range(k)]
    
    return min_return_cost, end_costs, coins

def main():
    data = sys.stdin.read().strip().split()
    if not data:
        return
    
    it = iter(data)
    n = int(next(it))
    start = (int(next(it)), int(next(it)), int(next(it)))
    
    switches = []
    coin_clusters = []
    cluster_return_costs = []  # cost to collect all coins and return to switch
    cluster_end_costs = []     # for each ending coin, cost to collect all coins
    
    for _ in range(n):
        k = int(next(it))
        sx, sy, sz = int(next(it)), int(next(it)), int(next(it))
        switches.append((sx, sy, sz))
        
        coins = []
        for _ in range(k):
            cx, cy, cz = int(next(it)), int(next(it)), int(next(it))
            coins.append((cx, cy, cz))
        coin_clusters.append(coins)
        
        # Precompute for this cluster
        return_cost, end_costs, _ = solve_cluster_tsp(switches[-1], coins)
        cluster_return_costs.append(return_cost)
        cluster_end_costs.append(end_costs)
    
    # Precompute all distances between switches
    switch_dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            d = dist(switches[i], switches[j])
            switch_dist[i][j] = d
            switch_dist[j][i] = d
    
    # Precompute start to switch distances
    start_to_switch = [dist(start, s) for s in switches]
    
    coin_to_switch = []
    for i in range(n):
        coins = coin_clusters[i]
        k = len(coins)
        cts = [[0.0] * n for _ in range(k)]
        for coin_idx in range(k):
            for j in range(n):
                if i == j:
                    cts[coin_idx][j] = 0.0
                else:
                    cts[coin_idx][j] = dist(coins[coin_idx], switches[j])
        coin_to_switch.append(cts)

    dp = [[[INF] * (max(len(c) for c in coin_clusters) + 1) for _ in range(n)] for _ in range(1<<n)]
    
    for i in range(n):
        # Start -> switch i -> collect coins -> end at switch
        cost = start_to_switch[i] + cluster_return_costs[i]
        dp[1<<i][i][0] = min(dp[1<<i][i][0], cost)  # last_coin_idx = -1 -> index 0
        
        # Start -> switch i -> collect coins -> end at coin j
        k = len(coin_clusters[i])
        for j in range(k):
            cost = start_to_switch[i] + cluster_end_costs[i][j]
            dp[1<<i][i][j+1] = min(dp[1<<i][i][j+1], cost)  # last_coin_idx = j -> index j+1
    
    # DP transitions
    for mask in range(1, 1<<n):
        for last_switch in range(n):
            if not (mask & (1<<last_switch)):
                continue
            
            k_last = len(coin_clusters[last_switch])
            for last_pos in range(k_last + 1):  # 0=switch, 1..k_last=coins
                current_cost = dp[mask][last_switch][last_pos]
                if current_cost == INF:
                    continue
                
                if last_pos == 0:
                    current_pos = switches[last_switch]
                else:
                    current_pos = coin_clusters[last_switch][last_pos-1]
                
                # Try to go to an unvisited switch
                for next_switch in range(n):
                    if mask & (1<<next_switch):
                        continue
                    
                    # Distance from current position to next switch
                    if last_pos == 0:
                        travel = switch_dist[last_switch][next_switch]
                    else:
                        travel = dist(current_pos, switches[next_switch])
                    
                    new_mask = mask | (1<<next_switch)
                    
                    # Option 1: At next switch, collect coins and end at switch
                    cost1 = current_cost + travel + cluster_return_costs[next_switch]
                    dp[new_mask][next_switch][0] = min(dp[new_mask][next_switch][0], cost1)
                    
                    # Option 2: At next switch, collect coins and end at coin j
                    k_next = len(coin_clusters[next_switch])
                    for j in range(k_next):
                        cost2 = current_cost + travel + cluster_end_costs[next_switch][j]
                        dp[new_mask][next_switch][j+1] = min(dp[new_mask][next_switch][j+1], cost2)
    
    # Find minimum cost to visit all switches
    full_mask = (1<<n) - 1
    answer = INF
    
    for last_switch in range(n):
        for last_pos in range(len(coin_clusters[last_switch]) + 1):
            answer = min(answer, dp[full_mask][last_switch][last_pos])
    
    if n == 0:
        print("0.0000000000")
    else:
        print(f"{answer:.10f}")

if __name__ == "__main__":
    main()