#!/usr/bin/env python
# compress.py

"""
compress.py run???? [options]

    Python script for simulation submission to a cluster. It creates submission bash files for Slurm.

    Arguments
    ---------
    run???? : list of directories
        List of directories to compress.

    Options
    -------
    queue : optional, chosen partition, default : debug
        submit will run the script.py in the specified partition.
    runtime : optional, default : 1-0:00
        maximum running time of the simulation. Syntax is day-hours:minutes:seconds


    Examples
    --------
    1.
        compress.py run????

    Functions

M. Le Verge--Serandour
Adapted from submit.py in chain_lumen
Creation : 01/04/18
"""


import os, sys
import subprocess
subcmd = 'sbatch'
queue = 'debug'
runtime = '1-0:00'
cpu_per_task = 1
nodelist = ''

folder_path = '/share/mathieu.leverge/git/chain_lumen/_ressources/'

def write_gen(directories, queue=queue, runtime=runtime, cpu_per_task=cpu_per_task, nodelist='') :
    nconfig=len(directories)
    filename='compress.sh'
    f = open(filename, 'w')
    f.write('#!/bin/bash\n')
    f.write('#SBATCH --job-name=compress\n')
    f.write('#SBATCH --ntasks=1\n')
    f.write('#SBATCH --nodes=1\n')
    f.write('#SBATCH --cpus-per-task='+ str(cpu_per_task) +'\n')
    if len(nodelist) > 0 :
        f.write('#SBATCH --nodelist=' + nodelist + '\n')
    else :
        f.write('#SBATCH --partition=' + queue + '\n')
    f.write('#SBATCH --time='+runtime+'\n')
    f.write('#SBATCH --mem=512\n')
    f.write('#SBATCH --signal=INT@60\n')
    f.write('#SBATCH --signal=TERM@120\n')
    f.write('#SBATCH --array=1-'+str(nconfig)+'\n')
    
    f.write('./comp$SLURM_ARRAY_TASK_ID\n')
    f.close()

    os.chmod(filename, 0700)
    return filename

def write_s(confname, dirconfig, n) :
    #global confname
    filename='comp' + str(n)
    f = open(filename, 'w')
    f.write('#!/bin/bash\n')
    f.write('tar cf ' + dirconfig + '.tar ' + dirconfig + ' \n')
    f.write('bzip2 ' + dirconfig + '.tar')
    
    os.chmod(filename, 0700)
    return 0

def main(args):
    global subcmd, queue, runtime, cpu_per_task, script, confname, nodelist
    list_dir = []

    for arg in args :
        if arg.endswith('.py') :
            script = os.path.abspath(arg)
        elif arg.startswith('queue=') :
            queue = arg[len('queue='):]
            
        elif arg.startswith('runtime=') :
            runtime = arg[len('runtime='):]
            
        elif arg.startswith('cpu-per-task=') :
            cpu_per_task = int(arg[len('cpu-per-task='):])
            if cpu_per_task < 1 :
                print('Error : too few cpus (0) per task assigned. Variable set to one.')
                cpu_per_task = 1
            elif cpu_per_task > 1 :
                print(str(cpu_per_task) + 'assigned for one task.')

        elif arg.startswith('nodelist=') :
            nodelist = arg[len('nodelist='):]
        else :
            list_dir += [arg]
    
    nconfig = len(list_dir)
    f = write_gen(list_dir, queue=queue, runtime=runtime, cpu_per_task=cpu_per_task, nodelist=nodelist)
    n = 1
    
    for dirconfig in list_dir :
        list_files = os.listdir(dirconfig)
        for elem in list_files :
            if elem.endswith('.conf') :
                confname = elem
        write_s(confname, dirconfig, n)
        n += 1
        
    submitname = 'compress.sh'
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