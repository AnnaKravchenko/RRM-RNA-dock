"""
This scripts assume that each docking directory contains:
  [motif].dat : the final dat file
  frag{}r.pdb : the bound form
  [motif].list : the list of fragments fro mthe library (coarse-grain)
"""
# input: 1. docking_dir with frags etc
# input: 2. $LIBRARY (nalib)
# input: 3. nproc

def run(*popenargs, **kwargs):
    input = kwargs.pop("input", None)
    check = kwargs.pop("handle", False)

    if input is not None:
        if 'stdin' in kwargs:
            raise ValueError('stdin and input arguments may not both be used.')
        kwargs['stdin'] = subprocess.PIPE

    process = subprocess.Popen(*popenargs, **kwargs)
    try:
        stdout, stderr = process.communicate(input)
    except:
        process.kill()
        process.wait()
        raise
    retcode = process.poll()
    if check and retcode:
        raise subprocess.CalledProcessError(
            retcode, process.args, output=stdout, stderr=stderr)
    return retcode, stdout, stderr


import sys, os, subprocess, multiprocessing
docking_dir = sys.argv[1]
if len(sys.argv) > 3:
    nprocesses = int(sys.argv[3])
else:
    nprocesses = -1
    print('Specifiy number of cpu!!!')
library = sys.argv[2]
os.chdir(docking_dir)
assert os.path.exists("boundfrag.list")
boundfrag = [(l.split()[0], l.split()[1]) for l in open("boundfrag.list") if len(l.strip())]
jobs = []
for frag, motif in boundfrag:
    bound = "frag{}r.pdb".format(frag)
    docked_poses = "{}-e7.dat".format(motif)
    ens_list = "{0}/{1}-clust1.0r.list".format(library, motif)
    assert os.path.exists(bound), bound
    assert os.path.exists(ens_list), ens_list
    assert os.path.exists(docked_poses), docked_poses
    #cmd = "conda activate attract"
    cmd = "python2 $ATTRACTDIR/lrmsd.py {0} `head -1 {2}` {1} --ens 2 {2} --allatoms |  cut -d' ' -f2- | nl -w1 -s' '  > frag{3}.lrmsd".format(
        docked_poses, bound, ens_list, frag )
    jobs.append(cmd)

def run_job(cmd):
    print(cmd)
    run(cmd, shell=True)

pool = multiprocessing.Pool(nprocesses)
pool.map(run_job, jobs)
