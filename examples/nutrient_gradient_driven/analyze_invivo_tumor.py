#!/usr/bin/env python3
"""
Analyze cell population and oxygen gradient results from the simulation.
This script reads the statistics file and plots cell count and oxygen levels over time.
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def analyze_oxygen_gradient(results_dir="results_inVitro"):
    """
    Analyze cell population dynamics in oxygen gradient.
    
    Args:
        results_dir (str): Directory containing simulation results
    """
    
    # Check if results directory exists
    if not os.path.exists(results_dir):
        print(f"Error: Results directory '{results_dir}' not found.")
        print("Please run the simulation first with: make invitro")
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
        
        # Check if N_cells column exists
        if 'N_cells' not in data:
            print("Error: N_cells column not found in statistics file.")
            print("Available columns:", list(data.keys()))
            return
        
        # Convert to numpy arrays for easier handling
        current_time = np.array(data['current_time'])
        n_cells = np.array(data['N_cells'])
        
        # Print basic statistics
        print("=== Oxygen Gradient Model Analysis ===")
        print(f"Initial cell count: {int(n_cells[0])}")
        print(f"Final cell count: {int(n_cells[-1])}")
        print(f"Maximum cell count: {int(n_cells.max())}")
        print(f"Cell count change: {int(n_cells[-1] - n_cells[0])}")
        if n_cells[0] > 0:
            print(f"Growth factor: {n_cells[-1] / n_cells[0]:.2f}x")
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Oxygen Gradient Model Analysis', fontsize=16)
        
        # Plot 1: Cell count over time
        axes[0, 0].plot(current_time, n_cells, 'b-', linewidth=2)
        axes[0, 0].set_xlabel('Time')
        axes[0, 0].set_ylabel('Number of Cells')
        axes[0, 0].set_title('Cell Population Over Time')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Growth rate (cell count change per time unit)
        if len(n_cells) > 1:
            cell_diff = np.diff(n_cells)
            time_diff = np.diff(current_time)
            growth_rate = cell_diff / time_diff
            
            axes[0, 1].plot(current_time[1:], growth_rate, 'g-', linewidth=2)
            axes[0, 1].set_xlabel('Time')
            axes[0, 1].set_ylabel('Growth Rate (cells/time)')
            axes[0, 1].set_title('Population Growth Rate')
            axes[0, 1].grid(True, alpha=0.3)
            axes[0, 1].axhline(y=0, color='r', linestyle='--', alpha=0.5)
        
        # Plot 3: Log scale cell count (if significant growth)
        if n_cells.max() / n_cells.min() > 2:
            axes[1, 0].semilogy(current_time, n_cells, 'r-', linewidth=2)
            axes[1, 0].set_xlabel('Time')
            axes[1, 0].set_ylabel('Number of Cells (log scale)')
            axes[1, 0].set_title('Cell Population (Log Scale)')
            axes[1, 0].grid(True, alpha=0.3)
        else:
            axes[1, 0].plot(current_time, n_cells, 'r-', linewidth=2)
            axes[1, 0].set_xlabel('Time')
            axes[1, 0].set_ylabel('Number of Cells')
            axes[1, 0].set_title('Cell Population')
            axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Cell count histogram (distribution over time)
        axes[1, 1].hist(n_cells, bins=30, color='purple', alpha=0.7, edgecolor='black')
        axes[1, 1].set_xlabel('Number of Cells')
        axes[1, 1].set_ylabel('Frequency')
        axes[1, 1].set_title('Cell Count Distribution')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save the plot
        plot_file = os.path.join(results_dir, "oxygen_gradient_analysis.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved as: {plot_file}")
        
        # Show the plot
        plt.show()
        
        # Save summary to file
        summary_file = os.path.join(results_dir, "population_summary.txt")
        with open(summary_file, 'w') as f:
            f.write("=== Oxygen Gradient Model Analysis Summary ===\n")
            f.write(f"Initial cell count: {int(n_cells[0])}\n")
            f.write(f"Final cell count: {int(n_cells[-1])}\n")
            f.write(f"Maximum cell count: {int(n_cells.max())}\n")
            f.write(f"Cell count change: {int(n_cells[-1] - n_cells[0])}\n")
            if n_cells[0] > 0:
                f.write(f"Growth factor: {n_cells[-1] / n_cells[0]:.2f}x\n")
            f.write(f"Simulation time: {current_time[-1]:.2f}\n")
        
        print(f"Summary saved as: {summary_file}")
        
    except Exception as e:
        print(f"Error analyzing results: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # Check if a results directory was provided as argument
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = "results_inVitro"
    
    analyze_oxygen_gradient(results_dir)
