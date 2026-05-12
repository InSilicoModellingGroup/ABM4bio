#!/usr/bin/env python
# coding: utf-8

import paramiko
import os
import argparse
import time
import re
import sys

# -------------------------------
# INPUTS FROM C++ (WITH DEFAULTS)
# -------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--step_num", type=int, default=0, help="Current increment")
parser.add_argument("--cell_count", type=int, default=1)
parser.add_argument("--contractile_force", type=float, default=10.0, help="Maximum contractile force applied by each cell")
parser.add_argument("--min_cell_radius", type=float, default=5.0, help="Minimum cell radius")
parser.add_argument("--max_cell_radius", type=float, default=40.0, help="Maximum cell radius")
parser.add_argument("--strut_radius", type=float, default=2.5, help="Radius of each strut - assumed homogeneous")
parser.add_argument("--delta_F", type=float, default=0.1, help="Perturbance force used to calculate k_ce and k_ecm")
parser.add_argument("--perturbance_dist", type=float, default=1000.0, help="Distance between cells being perturbed at the same time")
parser.add_argument("--random_state", type=int, default=0, help="Randomize the cells or repeat the last run positions and attachments")
parser.add_argument("--num_attachments", type=int, default=4, help="The maximum number of points a cell can attach to")
parser.add_argument("--verbose", action="store_true")
parser.add_argument("--private_key_path", type=str, required=True, help="Path to private key (.pem)")
parser.add_argument("--user_name", type=str, required=True, help="HPC username")
parser.add_argument("--host_name", type=str, required=True, help="HPC hostname")

args = parser.parse_args()

# assign to variables
step_num = args.step_num
cell_count = args.cell_count
contractile_force = args.contractile_force
min_cell_radius = args.min_cell_radius
max_cell_radius = args.max_cell_radius
strut_radius = args.strut_radius
delta_F = args.delta_F
perturbance_dist = args.perturbance_dist
random_state = args.random_state
num_attachments = args.num_attachments
verbose = args.verbose
private_key_path = args.private_key_path
user_name = args.user_name
host_name = args.host_name

# translate verbose into Slurm flag
verbose_flag = "--VERBOSE True" if verbose else ""

# Start by loading private key (.pem extension)

key = paramiko.RSAKey.from_private_key_file(private_key_path)

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

# Connect using the private key
ssh.connect(hostname=host_name, username=user_name, pkey=key)

try:
    # Run a simple command to confirm login (returns username)
    stdin, stdout, stderr = ssh.exec_command('whoami')
    output = stdout.read().decode().strip()
    error = stderr.read().decode().strip()

    if output:
        print(f"Logged in as: {output}")
    else:
        print(f"SSH command failed: {error}")
        
        
    # If step_num > 0 transfer the position of the cells to the HPC
    if step_num >= 0:
    	# Set up a SHH file transfer protocal (SFTP)
        sftp = ssh.open_sftp()
        local_get_file = f'./results/cell_positions/cells_t{step_num:04d}.csv'
        hpc_put_file = "/home/"+user_name+f"/FEM_SOLVER_folder/detached_cell_positions/detach_cells_step_{step_num}.csv"

        sftp.put(local_get_file, hpc_put_file)
        print(f"Copied local file {local_get_file} to remote path {hpc_put_file}")
    
        
    env_vars = {
        "STEP_NUM": step_num,
        "CELL_COUNT": cell_count,
        "CONTRACTILE_FORCE": contractile_force,
        "MIN_CELL_RADIUS": min_cell_radius,
        "MAX_CELL_RADIUS": max_cell_radius,
        "STRUT_RADIUS": strut_radius,
        "DELTA_F": delta_F,
        "PERTURBANCE_DIST": perturbance_dist,
        "RANDOM_STATE": random_state,
        "NUM_ATTACHMENTS": num_attachments,
        "VERBOSE": verbose_flag,  # either "--VERBOSE True" or ""
    }

    # Build the export string for sbatch
    env_export = " ".join(f'{key}="{value}"' for key, value in env_vars.items())
  
    # Command to submit SLURM job
    submit_cmd = f"cd FEM_SOLVER_folder && {env_export} sbatch SLURM_submit.sh"


    stdin, stdout, stderr = ssh.exec_command(submit_cmd)
    
    job_output = stdout.read().decode().strip()
    
    # Print stdout line by line in real time
    for line in iter(stdout.readline, ""):
    	print(line, end="", flush=True)

    # Optionally also print stderr in real time
    for line in iter(stderr.readline, ""):
    	print(line, end="", flush=True)
    	
    match = re.search(r"Submitted batch job (\d+)", job_output)
    job_id = match.group(1)
    
    # Wait for a Slurm job to finish by polling squeue.
    start_time = time.time()
    while True:
        stdin, stdout, stderr = ssh.exec_command(f"squeue -j {job_id} -h -o '%T'")
        status = stdout.read().decode().strip()
        
        elapsed = int(time.time() - start_time)

        if not status:
            sys.stdout.write(
                f"\rJob {job_id} finished after {elapsed}s.{' ' * 20}\n"
            )
            sys.stdout.flush()
            break
        
        sys.stdout.write(
            f"\rJob {job_id} status: {status} | elapsed: {elapsed}s"
        )
        sys.stdout.flush()
        
        #print(f"Job {job_id} status: {status} | elapsed: {elapsed}s", flush=True)

        time.sleep(30) # Check every 30 seconds 
 
        
    output_prefix="FEM_SOLVE"
        
    # Check Slurm job output for errors.
    out_file = f"{output_prefix}-{job_id}.out"
    stdin, stdout, stderr = ssh.exec_command(f"cd FEM_SOLVER_folder && tail -n 50 {out_file}")
    last_lines = stdout.read().decode()
    err_lines = stderr.read().decode()

    if err_lines:
        print("Errors found in stderr:")
        print(err_lines)

    # Simple check for Python traceback
    if "Traceback" in last_lines or "Error" in last_lines:
        print(f"Job {job_id} appears to have failed.")
        success = False
    else:
        print(f"Job {job_id} completed successfully.")
        success = True
        
    if not success:
        print("Job failed. Check the Slurm output for details.")
    else:
        print("Job completed successfully")
        out_file = f"FEM_SOLVER_folder/FEM_SOLVE-{job_id}.out"
        ssh.exec_command(f"rm -f {out_file}")
        
   

    # Set up a SHH file transfer protocal (SFTP)
    if step_num == 0:
    	sftp = ssh.open_sftp()
    
    # Make folder for results
    local_dir = f"./results/FEM/step_{step_num}"
    os.makedirs(local_dir, exist_ok=True)

    # Retrieve the lattice mesh file
    hpc_get_file = "/home/"+user_name+f"/FEM_SOLVER_folder/results/step_{step_num}/lattice.1d"
    local_put_file = local_dir + '/lattice.1d'

    sftp.get(hpc_get_file, local_put_file)
    print(f"Copied remote file {hpc_get_file} to local path {local_put_file}")

  
    # Retrieve the cell information file
    hpc_get_file = f"/home/"+user_name+f"/FEM_SOLVER_folder/results/step_{step_num}/cell_mechanics_step_{step_num}.dat"
    local_put_file = local_dir + f'/cell_mechanics_step_{step_num}.dat'

    sftp.get(hpc_get_file, local_put_file)
    print(f"Copied remote file {hpc_get_file} to local path {local_put_file}")


finally:
    sftp.close()
    stdout.channel.close()
    stderr.channel.close()
    ssh.close()






