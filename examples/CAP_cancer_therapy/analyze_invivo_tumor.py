#!/usr/bin/env python3
"""
Analyze tumor volume results from the in vivo 200mm³ tumor model simulation.
This script reads the statistics file and plots tumor volume over time.
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def analyze_tumor_volume(results_dir="results_inVivo"):
    """
    Analyze tumor volume from simulation results.
    
    Args:
        results_dir (str): Directory containing simulation results
    """
    
    # Check if results directory exists
    if not os.path.exists(results_dir):
        print(f"Error: Results directory '{results_dir}' not found.")
        print("Please run the simulation first with: make invivo")
        return
    
    # Look for statistics file
    stats_file = os.path.join(results_dir, "stats.csv")
    if not os.path.exists(stats_file):
        print(f"Error: Statistics file '{stats_file}' not found.")
        return
    
    try:
        # Read the statistics file using CSV module
        data = {}
        with open(stats_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                for key, value in row.items():
                    # Strip whitespace from column names
                    clean_key = key.strip()
                    if clean_key not in data:
                        data[clean_key] = []
                    try:
                        # Try to convert to float, keep as string if it fails
                        data[clean_key].append(float(value))
                    except ValueError:
                        data[clean_key].append(value)
        
        # Check if tumor_volume_mm3 column exists
        if 'tumor_volume_mm3' not in data:
            print("Error: tumor_volume_mm3 column not found in statistics file.")
            print("Available columns:", list(data.keys()))
            return
        
        # Convert to numpy arrays for easier handling
        tumor_volume = np.array(data['tumor_volume_mm3'])
        current_time = np.array(data['current_time'])
        n_cells = np.array(data['N_cells'])
        
        # Print basic statistics
        print("=== Tumor Volume Analysis ===")
        print(f"Initial tumor volume: {tumor_volume[0]:.2f} mm³")
        print(f"Final tumor volume: {tumor_volume[-1]:.2f} mm³")
        print(f"Maximum tumor volume: {tumor_volume.max():.2f} mm³")
        print(f"Volume change: {tumor_volume[-1] - tumor_volume[0]:.2f} mm³")
        print(f"Growth factor: {tumor_volume[-1] / tumor_volume[0]:.2f}x")
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('In Vivo Tumor Model Analysis (200mm³)', fontsize=16)
        
        # Plot 1: Tumor volume over time
        axes[0, 0].plot(current_time, tumor_volume, 'b-', linewidth=2)
        axes[0, 0].set_xlabel('Time')
        axes[0, 0].set_ylabel('Tumor Volume (mm³)')
        axes[0, 0].set_title('Tumor Volume Over Time')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Cell count over time
        axes[0, 1].plot(current_time, n_cells, 'r-', linewidth=2)
        axes[0, 1].set_xlabel('Time')
        axes[0, 1].set_ylabel('Number of Cells')
        axes[0, 1].set_title('Cell Count Over Time')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Volume vs Cell count correlation
        axes[1, 0].scatter(n_cells, tumor_volume, alpha=0.6, s=10)
        axes[1, 0].set_xlabel('Number of Cells')
        axes[1, 0].set_ylabel('Tumor Volume (mm³)')
        axes[1, 0].set_title('Volume vs Cell Count')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Growth rate (volume change per time unit)
        if len(tumor_volume) > 1:
            volume_diff = np.diff(tumor_volume)
            time_diff = np.diff(current_time)
            growth_rate = volume_diff / time_diff
            
            axes[1, 1].plot(current_time[1:], growth_rate, 'g-', linewidth=2)
            axes[1, 1].set_xlabel('Time')
            axes[1, 1].set_ylabel('Volume Growth Rate (mm³/time)')
            axes[1, 1].set_title('Tumor Growth Rate')
            axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save the plot
        plot_file = os.path.join(results_dir, "tumor_volume_analysis.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved as: {plot_file}")
        
        # Show the plot
        plt.show()
        
        # Save summary to file
        summary_file = os.path.join(results_dir, "volume_summary.txt")
        with open(summary_file, 'w') as f:
            f.write("=== In Vivo Tumor Volume Analysis Summary ===\n")
            f.write(f"Initial tumor volume: {tumor_volume[0]:.2f} mm³\n")
            f.write(f"Final tumor volume: {tumor_volume[-1]:.2f} mm³\n")
            f.write(f"Maximum tumor volume: {tumor_volume.max():.2f} mm³\n")
            f.write(f"Volume change: {tumor_volume[-1] - tumor_volume[0]:.2f} mm³\n")
            f.write(f"Growth factor: {tumor_volume[-1] / tumor_volume[0]:.2f}x\n")
            f.write(f"Final cell count: {int(n_cells[-1])}\n")
            f.write(f"Simulation time: {current_time[-1]:.2f}\n")
        
        print(f"Summary saved as: {summary_file}")
        
    except Exception as e:
        print(f"Error analyzing results: {e}")

if __name__ == "__main__":
    # Check if a results directory was provided as argument
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = "results_inVivo"
    
    analyze_tumor_volume(results_dir)
