
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def create_boxplots():
    """Create boxplots for all benchmark results"""
    
    # Part (a) - HW5
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Part (a) - HW5 Benchmark Results', fontsize=16, fontweight='bold')
    
    # Load HW5 data
    hw5_data = pd.read_csv('hw5_benchmark_results.csv')
    hw5_smooth_data = pd.read_csv('hw5_smooth_benchmark_results.csv')
    
    # Success Rate
    ax = axes[0, 0]
    params = [f"({row['n']}, {row['r']})" for _, row in hw5_data.iterrows()]
    success_rates = hw5_data['success_rate'].values * 100
    success_rates_smooth = hw5_smooth_data['success_rate'].values * 100
    
    x = np.arange(len(params))
    width = 0.35
    ax.bar(x - width/2, success_rates, width, label='Without Smoothing', alpha=0.8)
    ax.bar(x + width/2, success_rates_smooth, width, label='With Smoothing', alpha=0.8)
    ax.set_ylabel('Success Rate (%)')
    ax.set_title('Success Rate vs (n, r)')
    ax.set_xticks(x)
    ax.set_xticklabels(params, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Path Length
    ax = axes[0, 1]
    path_lengths = hw5_data['avg_path_length'].values
    path_lengths_smooth = hw5_smooth_data['avg_path_length'].values
    
    ax.bar(x - width/2, path_lengths, width, label='Without Smoothing', alpha=0.8)
    ax.bar(x + width/2, path_lengths_smooth, width, label='With Smoothing', alpha=0.8)
    ax.set_ylabel('Average Path Length')
    ax.set_title('Path Length vs (n, r)')
    ax.set_xticks(x)
    ax.set_xticklabels(params, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Computation Time
    ax = axes[0, 2]
    times = hw5_data['avg_time_ms'].values
    times_smooth = hw5_smooth_data['avg_time_ms'].values
    
    ax.bar(x - width/2, times, width, label='Without Smoothing', alpha=0.8)
    ax.bar(x + width/2, times_smooth, width, label='With Smoothing', alpha=0.8)
    ax.set_ylabel('Average Time (ms)')
    ax.set_title('Computation Time vs (n, r)')
    ax.set_xticks(x)
    ax.set_xticklabels(params, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Effect of n (fixing r)
    ax = axes[1, 0]
    for r_val in [0.5, 1.0, 1.5, 2.0]:
        subset = hw5_data[hw5_data['r'] == r_val]
        ax.plot(subset['n'], subset['success_rate'] * 100, marker='o', label=f'r={r_val}')
    ax.set_xlabel('Number of Samples (n)')
    ax.set_ylabel('Success Rate (%)')
    ax.set_title('Effect of n on Success Rate')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Effect of r (fixing n)
    ax = axes[1, 1]
    for n_val in [200, 500]:
        subset = hw5_data[hw5_data['n'] == n_val]
        ax.plot(subset['r'], subset['avg_path_length'], marker='o', label=f'n={n_val}')
    ax.set_xlabel('Connection Radius (r)')
    ax.set_ylabel('Average Path Length')
    ax.set_title('Effect of r on Path Length')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Smoothing Improvement
    ax = axes[1, 2]
    improvement = ((path_lengths - path_lengths_smooth) / path_lengths) * 100
    ax.bar(x, improvement, alpha=0.8, color='green')
    ax.set_ylabel('Path Length Reduction (%)')
    ax.set_title('Smoothing Improvement')
    ax.set_xticks(x)
    ax.set_xticklabels(params, rotation=45, ha='right')
    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('hw5_benchmarks.png', dpi=300, bbox_inches='tight')
    print("Saved hw5_benchmarks.png")
    
    # Part (b) - HW2 Workspaces
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Part (b) - HW2 Workspace Benchmark Results', fontsize=16, fontweight='bold')
    
    # Load HW2 data
    w1_data = pd.read_csv('w1_benchmark_results.csv')
    w1_smooth_data = pd.read_csv('w1_smooth_benchmark_results.csv')
    w2_data = pd.read_csv('w2_benchmark_results.csv')
    w2_smooth_data = pd.read_csv('w2_smooth_benchmark_results.csv')
    
    params_w = [f"({row['n']}, {row['r']})" for _, row in w1_data.iterrows()]
    x = np.arange(len(params_w))
    
    # W1 Success Rate
    ax = axes[0, 0]
    sr_w1 = w1_data['success_rate'].values * 100
    sr_w1_smooth = w1_smooth_data['success_rate'].values * 100
    
    ax.bar(x - width/2, sr_w1, width, label='Without Smoothing', alpha=0.8)
    ax.bar(x + width/2, sr_w1_smooth, width, label='With Smoothing', alpha=0.8)
    ax.set_ylabel('Success Rate (%)')
    ax.set_title('W1 Success Rate vs (n, r)')
    ax.set_xticks(x)
    ax.set_xticklabels(params_w, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # W2 Success Rate
    ax = axes[0, 1]
    sr_w2 = w2_data['success_rate'].values * 100
    sr_w2_smooth = w2_smooth_data['success_rate'].values * 100
    
    ax.bar(x - width/2, sr_w2, width, label='Without Smoothing', alpha=0.8)
    ax.bar(x + width/2, sr_w2_smooth, width, label='With Smoothing', alpha=0.8)
    ax.set_ylabel('Success Rate (%)')
    ax.set_title('W2 Success Rate vs (n, r)')
    ax.set_xticks(x)
    ax.set_xticklabels(params_w, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Comparison W1 vs W2
    ax = axes[0, 2]
    ax.bar(x - width/2, sr_w1, width, label='W1', alpha=0.8)
    ax.bar(x + width/2, sr_w2, width, label='W2', alpha=0.8)
    ax.set_ylabel('Success Rate (%)')
    ax.set_title('W1 vs W2 Success Rate')
    ax.set_xticks(x)
    ax.set_xticklabels(params_w, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Path lengths
    ax = axes[1, 0]
    pl_w1 = w1_data['avg_path_length'].values
    pl_w1_smooth = w1_smooth_data['avg_path_length'].values
    
    ax.bar(x - width/2, pl_w1, width, label='Without Smoothing', alpha=0.8)
    ax.bar(x + width/2, pl_w1_smooth, width, label='With Smoothing', alpha=0.8)
    ax.set_ylabel('Average Path Length')
    ax.set_title('W1 Path Length vs (n, r)')
    ax.set_xticks(x)
    ax.set_xticklabels(params_w, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Computation times
    ax = axes[1, 1]
    time_w1 = w1_data['avg_time_ms'].values
    time_w2 = w2_data['avg_time_ms'].values
    
    ax.bar(x - width/2, time_w1, width, label='W1', alpha=0.8)
    ax.bar(x + width/2, time_w2, width, label='W2', alpha=0.8)
    ax.set_ylabel('Average Time (ms)')
    ax.set_title('Computation Time: W1 vs W2')
    ax.set_xticks(x)
    ax.set_xticklabels(params_w, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Overall comparison
    ax = axes[1, 2]
    metrics = ['W1 Success', 'W2 Success', 'W1 Time', 'W2 Time']
    # Normalize to 0-100 scale for comparison
    w1_success_norm = np.mean(sr_w1)
    w2_success_norm = np.mean(sr_w2)
    w1_time_norm = 100 - (np.mean(time_w1) / max(np.mean(time_w1), np.mean(time_w2)) * 100)
    w2_time_norm = 100 - (np.mean(time_w2) / max(np.mean(time_w1), np.mean(time_w2)) * 100)
    
    values = [w1_success_norm, w2_success_norm, w1_time_norm, w2_time_norm]
    colors = ['blue', 'orange', 'blue', 'orange']
    
    ax.bar(range(len(metrics)), values, color=colors, alpha=0.8)
    ax.set_ylabel('Normalized Score (0-100)')
    ax.set_title('Overall Performance Comparison')
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('hw2_benchmarks.png', dpi=300, bbox_inches='tight')
    print("Saved hw2_benchmarks.png")
    
    plt.show()

if __name__ == "__main__":
    create_boxplots()