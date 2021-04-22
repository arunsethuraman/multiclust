import multiprocessing
import time
import random
import subprocess

unique_filename = "str"
lock = multiprocessing.Lock()

def multiclust(n):
    time.sleep(random.randint(1,3))
    result = subprocess.run(["./multiclust", "-f", "SS2alleles001.str", "-s", "2", "-M"], capture_output=True, text=True)
    
    ##prints worker number and then multiclust output, but prints output after every single worker number
    ##print("This is worker [{0}]:\n".format(n) + result.stdout)

    lock.acquire()
    f = open(str(unique_filename) + "_stdout.txt", "a+")
    f.write("[{0}]".format(n) + " " + result.stdout + "\n")
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

    g = open(str(unique_filename) + "_stdout.txt", "w")
    g.write(worker_num_line + "\nBest log likelihood: " + best_logL_line)

processes = [ ]
for i in range(10):
    t = multiprocessing.Process(target=multiclust, args=(i,))
    processes.append(t)
    t.start()

for one_process in processes:
    one_process.join()

master()

print("Done!")
