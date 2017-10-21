# usage: "python3 run_mpi_tests.py"

import os, sys
import time
import datetime

start_time = time.time()


input_file = "test_StructuralMechnaicsApplication_mpi.py"

list_processors = []

for i in range(1,31):
    list_processors.append(i)

list_processors.extend([40, 50])

os.system("export OMP_NUM_THREADS=1")
print("OMP threads set to 1")

# parallel executions
for num_processors in list_processors:
    system_cmd = "mpiexec -np " + str(num_processors) + " python3 " + input_file
    os.system(system_cmd)

test_runtime = datetime.timedelta(seconds=int((time.time() - start_time)))
print("Runtime: " + str(test_runtime))
