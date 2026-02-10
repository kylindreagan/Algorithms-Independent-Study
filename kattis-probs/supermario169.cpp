//https://open.kattis.com/problems/supermario169

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <sstream>

using namespace std;

const double INF = numeric_limits<double>::max();

struct Point {
    int x, y, z;
    Point() : x(0), y(0), z(0) {}
    Point(int x, int y, int z) : x(x), y(y), z(z) {}
};

double dist_sq(const Point& a, const Point& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return dx*dx + dy*dy + dz*dz;
}

double dist(const Point& a, const Point& b) {
    return sqrt(dist_sq(a, b));
}

struct ClusterResult {
    double return_cost;
    vector<double> end_costs;
    vector<Point> coins;
};

ClusterResult solve_cluster_tsp(const Point& switch_pos, const vector<Point>& coins) {
    int k = coins.size();
    ClusterResult result;
    result.coins = coins;
    
    if (k == 0) {
        result.return_cost = 0.0;
        result.end_costs.push_back(0.0);
        return result;
    }
    
    // Distance from switch to each coin
    vector<double> d0(k);
    for (int i = 0; i < k; i++) {
        d0[i] = dist(switch_pos, coins[i]);
    }
    
    // Distance between coins
    vector<vector<double>> d(k, vector<double>(k, 0.0));
    for (int i = 0; i < k; i++) {
        for (int j = i + 1; j < k; j++) {
            double d_ij = dist(coins[i], coins[j]);
            d[i][j] = d_ij;
            d[j][i] = d_ij;
        }
    }
    
    // DP for TSP within cluster
    int full_mask_size = 1 << k;
    vector<vector<double>> dp(full_mask_size, vector<double>(k, INF));
    
    // Initialize with starting at each coin
    for (int i = 0; i < k; i++) {
        dp[1 << i][i] = d0[i];
    }
    
    // DP transitions
    for (int mask = 1; mask < full_mask_size; mask++) {
        for (int last = 0; last < k; last++) {
            if (!(mask & (1 << last)) || dp[mask][last] == INF) {
                continue;
            }
            
            for (int nxt = 0; nxt < k; nxt++) {
                if (mask & (1 << nxt)) {
                    continue;
                }
                int new_mask = mask | (1 << nxt);
                double new_cost = dp[mask][last] + d[last][nxt];
                if (dp[new_mask][nxt] > new_cost) {
                    dp[new_mask][nxt] = new_cost;
                }
            }
        }
    }
    
    int full_mask = (1 << k) - 1;
    
    // Cost to return to switch from each coin
    vector<double> return_to_switch(k);
    for (int i = 0; i < k; i++) {
        return_to_switch[i] = dist(coins[i], switch_pos);
    }
    
    // Minimum cost to collect all coins and return to switch
    double min_return_cost = INF;
    for (int i = 0; i < k; i++) {
        double total = dp[full_mask][i] + return_to_switch[i];
        if (total < min_return_cost) {
            min_return_cost = total;
        }
    }
    result.return_cost = min_return_cost;
    
    // Cost to end at each coin (without returning to switch)
    result.end_costs.resize(k);
    for (int i = 0; i < k; i++) {
        result.end_costs[i] = dp[full_mask][i];
    }
    
    return result;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int n;
    cin >> n;
    
    Point start;
    cin >> start.x >> start.y >> start.z;
    
    vector<Point> switches(n);
    vector<vector<Point>> coin_clusters(n);
    vector<double> cluster_return_costs(n);
    vector<vector<double>> cluster_end_costs(n);
    
    // Read input and precompute cluster solutions
    for (int i = 0; i < n; i++) {
        int k;
        cin >> k;
        cin >> switches[i].x >> switches[i].y >> switches[i].z;
        
        vector<Point> coins(k);
        for (int j = 0; j < k; j++) {
            cin >> coins[j].x >> coins[j].y >> coins[j].z;
        }
        coin_clusters[i] = coins;
        
        // Solve TSP for this cluster
        ClusterResult cluster_result = solve_cluster_tsp(switches[i], coins);
        cluster_return_costs[i] = cluster_result.return_cost;
        cluster_end_costs[i] = cluster_result.end_costs;
    }
    
    if (n == 0) {
        cout << "0.0000000000\n";
        return 0;
    }
    
    // Precompute distances between switches
    vector<vector<double>> switch_dist(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d = dist(switches[i], switches[j]);
            switch_dist[i][j] = d;
            switch_dist[j][i] = d;
        }
    }
    
    // Precompute distances from start to each switch
    vector<double> start_to_switch(n);
    for (int i = 0; i < n; i++) {
        start_to_switch[i] = dist(start, switches[i]);
    }
    
    // Precompute distances from coins to switches
    vector<vector<vector<double>>> coin_to_switch(n);
    for (int i = 0; i < n; i++) {
        int k = coin_clusters[i].size();
        vector<vector<double>> cts(k, vector<double>(n, 0.0));
        for (int coin_idx = 0; coin_idx < k; coin_idx++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    cts[coin_idx][j] = 0.0;
                } else {
                    cts[coin_idx][j] = dist(coin_clusters[i][coin_idx], switches[j]);
                }
            }
        }
        coin_to_switch[i] = cts;
    }
    
    // Find maximum number of coins in any cluster
    int max_coins = 0;
    for (const auto& cluster : coin_clusters) {
        max_coins = max(max_coins, (int)cluster.size());
    }
    
    // DP: dp[mask][last_switch][last_pos]
    // last_pos = 0: at switch, last_pos > 0: at coin (last_pos-1)
    int dp_size = 1 << n;
    vector<vector<vector<double>>> dp(
        dp_size, 
        vector<vector<double>>(
            n, 
            vector<double>(max_coins + 1, INF)
        )
    );
    
    // Initialize DP
    for (int i = 0; i < n; i++) {
        // Start -> switch i -> collect coins -> end at switch
        double cost = start_to_switch[i] + cluster_return_costs[i];
        dp[1 << i][i][0] = min(dp[1 << i][i][0], cost);
        
        // Start -> switch i -> collect coins -> end at coin j
        int k = coin_clusters[i].size();
        for (int j = 0; j < k; j++) {
            double cost = start_to_switch[i] + cluster_end_costs[i][j];
            dp[1 << i][i][j + 1] = min(dp[1 << i][i][j + 1], cost);
        }
    }
    
    // DP transitions
    for (int mask = 1; mask < dp_size; mask++) {
        for (int last_switch = 0; last_switch < n; last_switch++) {
            if (!(mask & (1 << last_switch))) {
                continue;
            }
            
            int k_last = coin_clusters[last_switch].size();
            for (int last_pos = 0; last_pos <= k_last; last_pos++) {
                double current_cost = dp[mask][last_switch][last_pos];
                if (current_cost == INF) {
                    continue;
                }
                
                Point current_pos;
                if (last_pos == 0) {
                    current_pos = switches[last_switch];
                } else {
                    current_pos = coin_clusters[last_switch][last_pos - 1];
                }
                
                // Try to go to an unvisited switch
                for (int next_switch = 0; next_switch < n; next_switch++) {
                    if (mask & (1 << next_switch)) {
                        continue;
                    }
                    
                    // Distance from current position to next switch
                    double travel;
                    if (last_pos == 0) {
                        travel = switch_dist[last_switch][next_switch];
                    } else {
                        travel = dist(current_pos, switches[next_switch]);
                    }
                    
                    int new_mask = mask | (1 << next_switch);
                    
                    // Option 1: At next switch, collect coins and end at switch
                    double cost1 = current_cost + travel + cluster_return_costs[next_switch];
                    if (dp[new_mask][next_switch][0] > cost1) {
                        dp[new_mask][next_switch][0] = cost1;
                    }
                    
                    // Option 2: At next switch, collect coins and end at coin j
                    int k_next = coin_clusters[next_switch].size();
                    for (int j = 0; j < k_next; j++) {
                        double cost2 = current_cost + travel + cluster_end_costs[next_switch][j];
                        if (dp[new_mask][next_switch][j + 1] > cost2) {
                            dp[new_mask][next_switch][j + 1] = cost2;
                        }
                    }
                }
            }
        }
    }
    
    // Find minimum cost to visit all switches
    int full_mask = (1 << n) - 1;
    double answer = INF;
    
    for (int last_switch = 0; last_switch < n; last_switch++) {
        int k_last = coin_clusters[last_switch].size();
        for (int last_pos = 0; last_pos <= k_last; last_pos++) {
            answer = min(answer, dp[full_mask][last_switch][last_pos]);
        }
    }
    
    cout << fixed << setprecision(10) << answer << "\n";
    
    return 0;
}