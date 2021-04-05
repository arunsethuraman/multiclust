import multiprocessing
import time
import random
import subprocess

def multiclust(n):
    time.sleep(random.randint(1,3))
    filename = "Blandings-PhyloStat.stru"
    result = subprocess.run(["./multiclust", "-f", filename, "-a", "-s", "2", "-M", "1"], capture_output=True, text=True)
    ## prints worker number without multiclust output
    ##print("This is worker [{0}]".format(n))
    ##prints worker number and then multiclust output, but prints output after every single worker number
    print("This is worker [{0}]:\n".format(n) + result.stdout)
    f = open(str(filename) + "_stdout.txt", "a+")
    f.write("[{0}]".format(n) + ": " + result.stdout + "\n")

processes = [ ]
for i in range(100):
    t = multiprocessing.Process(target=multiclust, args=(i,))
    processes.append(t)
    t.start()

for one_process in processes:
    one_process.join()

print("Done!")

## to create individually named output files, look at write_data() in write_file.c