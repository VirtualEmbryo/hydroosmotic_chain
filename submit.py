#!/usr/bin/env python
# submit.py

"""
submit.py script.py [options] config????

    Python script for simulation submission to a cluster. It creates submission bash files for Slurm.

    Arguments
    ---------
    script.py : python file
        The script to execute for data analysis. Must be a python script (.py).
    config???? : list of directories
        List of directories where the script.py has to be executed.

    Options
    -------
    queue : optional, chosen partition, default : debug
        submit will run the script.py in the specified partition.


    Examples
    --------

    Functions

M. Le Verge--Serandour
Adapted from submit.py in cytosim_public
Creation : 28/09/18
"""


import os, sys
import subprocess
subcmd = 'sbatch'
queue = 'debug'

folder_path = '/share/mathieu.leverge/git/cavitation/network/'

def write_gen(directories, queue=queue) :
    nconfig=len(directories)
    filename='sim.sh'
    f = open(filename, 'w')
    f.write('#!/bin/bash\n')
    f.write('#SBATCH --job-name=net\n')
    f.write('#SBATCH --ntasks=1\n')
    f.write('#SBATCH --nodes=1\n')
    f.write('#SBATCH --cpus-per-task=1\n')
    f.write('#SBATCH --partition=' + queue + '\n')
    f.write('#SBATCH --time=1-0:00\n')
    f.write('#SBATCH --mem=512\n')
    f.write('#SBATCH --signal=INT@60\n')
    f.write('#SBATCH --signal=TERM@120\n')
    f.write('#SBATCH --output=logs/out/out%a.txt\n')
    f.write('#SBATCH --error=logs/err/err%a.txt\n')
    f.write('#SBATCH --array=1-'+str(nconfig)+'\n')
    
    f.write('export OPENBLAS_NUM_THREADS=1\n')
    
    f.write('./job$SLURM_ARRAY_TASK_ID')
    f.close()

    os.chmod(filename, 0700)
    return filename

def write_s(config, n) :
    global folder_path
    filename='job' + str(n)
    f = open(filename, 'w')
    f.write('#!/bin/bash\n')
    f.write('cd ' + config + '\n')
    f.write('/share/mathieu.leverge/git/cavitation/network/simulation.py test.ini')
    
    os.chmod(filename, 0700)
    return 0

def main(args):
    global subcmd, queue
    list_dir = []

    for arg in args :
        if arg.startswith('queue=') :
            queue = arg[len('queue='):]
        else :
            list_dir += [arg]


    try :
        os.makedirs('logs')
        os.makedirs('logs/err')
        os.makedirs('logs/out')
    except : pass
    
    nconfig = len(list_dir)
    f = write_gen(list_dir, queue=queue)
    n = 1
    
    for d in list_dir :
        write_s(d, n)
        n += 1
        
    submitname = 'sim.sh'
    subprocess.call([subcmd + ' ' + submitname], shell = True)
    return ;


if __name__ == "__main__" :

    if len(sys.argv)<1 or sys.argv[1] == 'help' :
        print(__doc__)
    else :
	if len(sys.argv) > 1 :
            args = sys.argv[1:]
            main(args=args)
        else :
            main()