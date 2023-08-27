import os
import subprocess 
import sys
import math

size = 100 if len(sys.argv) == 1 else int(sys.argv[1])

print("\n---------------------------------------")
print("Serial Version...")  
print("---------------------------------------\n")

os.chdir("serial")
os.makedirs("results", exist_ok=True)
subprocess.run(["make"])
subprocess.run(["./infoli_simple", str(size)])

print("\n---------------------------------------") 
print("Parallel Version...")
print("---------------------------------------\n")

os.chdir("../parallel") 
os.makedirs("results", exist_ok=True)
subprocess.run(["make"])
subprocess.run(["./infoli_simple", str(size)])

print("\n---------------------------------------")
print("Comparing results...")
print("---------------------------------------\n") 

os.chdir("..")

with open('serial/results/lastStateDump.txt') as f1, open('parallel/results/lastStateDump.txt') as f2:
  close = True
  for line1, line2 in zip(f1, f2):
    values1 = [float(x) for x in line1.split()]
    values2 = [float(x) for x in line2.split()]
    
    if len(values1) != len(values2):
      print("Line lengths differ")
      break

    if not all(math.isclose(v1, v2, rel_tol=1e-5) for v1, v2 in zip(values1, values2)):
      print("Values not close")
      close = False
      break
  if close:
    print("Checks passed: output is equal.")