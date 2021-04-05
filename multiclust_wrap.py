import subprocess
import argparse

result = subprocess.run(["./multiclust", "-f", "SS2alleles001.str", "-a","-M", "1"], capture_output=True, text=True)

f = open("output_stdout.txt", "w+")
f.write(result.stdout)

print(result.stdout)
