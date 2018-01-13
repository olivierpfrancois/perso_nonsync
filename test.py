# multiproc_test.py

import sys
# Block python from writing pyc files
sys.dont_write_bytecode = True

import random
from multiprocessing.dummy import Pool as ThreadPool
import time


def f(x):
    return [x * x, x * 2]


if __name__ == '__main__':
    
    for i in range(1, 10):
        print i
        t0 = time.time()
        
        p = ThreadPool(i)
        test = p.map(f, range(1, 10000000))
        
        t1 = time.time()
    
        total = t1 - t0
        print(total)
        
        print type(test)
        print test[0]
        
        p.close()
        p.join()
    
'''
def list_append(count, id, out_list):
    """
    Creates an empty list and then appends a 
    random number to the list 'count' number
    of times. A CPU-heavy operation!
    """
    for i in range(count):
        out_list.append(random.random())


if __name__ == "__main__":
    size = 10000000  # Number of random numbers to add
    procs = 2  # Number of processes to create

    # Create a list of jobs and then iterate through
    # the number of processes appending each process to
    # the job list 
    jobs = []
    for i in range(0, procs):
        out_list = list()
        process = multiprocessing.Process(target=list_append,
                                          args=(size, i, out_list))
        jobs.append(process)

    # Start the processes (i.e. calculate the random number lists)        
    for j in jobs:
        j.start()

    # Ensure all of the processes have finished
    for j in jobs:
        j.join()

    print "List processing complete."
'''
