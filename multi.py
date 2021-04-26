import multiprocessing
import time
import random
import subprocess

unique_filename = "test_unique"
lock = multiprocessing.Lock()
start = time.perf_counter()

def multiclust(n):
    time.sleep(random.randint(1,3))
    result = subprocess.run(["./multiclust", "-f", "SS2alleles001.str", "-s", "2", "-M", "-o", unique_filename], capture_output=True, text=True)
    
    worker_start = time.perf_counter()
    ##prints worker number and then multiclust output, but prints output after every single worker number
    ##print("This is worker [{0}]:\n".format(n) + result.stdout)

    lock.acquire()
    f = open(str(unique_filename) + "_stdout.txt", "a+")
    
    ##clock function added to each worker to determine worker speed
    worker_finish = time.perf_counter()
    f.write("[{0}]".format(n) + " " + result.stdout + "Worker [{0}] time to complete: ".format(n) + str(worker_finish-worker_start) + " seconds\n\n")
    lock.release()

def master():
    worker_num = 0
    worker_num_line = 0
    current_worker = 0
    best_logL = -12345678901234567890
    best_logL_line = 0
    
    with open(str(unique_filename) + "_stdout.txt") as f:
        lines = f.readlines()
    for x in range(0, len(lines)):
        if lines[x][0] == '[':
            current_worker = lines[x][1]
        if lines[x][0] == '-':
            if float(lines[x]) > best_logL:
                best_logL = float(lines[x])
                best_logL_line = lines[x]
                worker_num = current_worker
                worker_num_line = lines[x-1]

    f.close()

    ##prints out worker information and logL before writing over file, used for debugging purposes
    ##print(worker_num, current_worker, best_logL, worker_num_line, best_logL_line)
    
    ##master will open a new file instead writing over all old worker files for debugging purposes
    ##g = open(str(unique_filename) + "_stdout.txt", "w")
    g = open(str(unique_filename) + "_master.txt", "w")
    g.write("Worker number: " + worker_num_line + "\nBest log likelihood: " + best_logL_line)

    

processes = [ ]
for i in range(10):
    t = multiprocessing.Process(target=multiclust, args=(i,))
    processes.append(t)
    t.start()

for one_process in processes:
    one_process.join()

master()

## time.perf_counter() is a precise clock we can use to measure runtime of workers + masters down to nanoseconds, used for debug
finish = time.perf_counter()
h = open(str(unique_filename) + "_master.txt", "a")
h.write("\nTotal time to complete: " + str(finish-start) + " seconds")
print("Done! Time to complete: " + str(finish-start) + " seconds")


