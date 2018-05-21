#! /usr/bin/env python
import sys
import re
from subprocess import Popen,PIPE
import re

# Slurm states
SLURM_TO_STATE = {
    'COMPLETED': 'success',
    'FAILED': 'failed',
    'CANCELLED': 'failed',
    'TIMEOUT': 'failed',
    'PREEMPTED': 'failed',
    'NODE_FAIL': 'failed',
    
    'PENDING': 'running',
    'RUNNING': 'running',
    'SUSPENDED': 'running',
    'COMPLETING': 'running',
    'CONFIGURING': 'running',
    'REQUEUED': 'running',

    'REVOKED': 'failed',
    'SPECIAL_EXIT': 'failed',
}
DEFAULT_STATE = 'running'

def jobstatus(jobid):
    try:
        o,e = Popen('sacct -u $(whoami) -nXPo state -j ' + sys.argv[1],
                    shell=True, stdout=PIPE, stderr=PIPE).communicate()
        ret = '{}'.format(str(o, 'utf-8').split()[0])
        if ret in SLURM_TO_STATE:
            return (SLURM_TO_STATE[ret], ret)
        else:
            print('Unknown slurm state: {}'.format(ret), file=sys.stderr)
            return (DEFAULT_STATE, ret)
    except Exception as e:
        print(e, file=sys.stderr)
        return (DEFAULT_STATE, None)

if __name__ == '__main__':
    jobid = sys.argv[1].strip()
    if not re.match('\d+', jobid):
        print("Argument is not job ID: {}".format(jobid), file=sys.stderr)
        print(DEFAULT_STATE)
    else:
        state, ret = jobstatus(jobid)     
        print(state)
