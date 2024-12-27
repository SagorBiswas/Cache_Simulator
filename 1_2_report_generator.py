import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import sys
import copy
import math

# # Ensure enough arguments are passed
# if len(sys.argv) < 11:
#     print("Usage: python script.py <block_size> <L1_size> <L1_assoc> <L2_size> <L2_assoc> <policy> <inclusion> <trace_file> <Hit_time_L1> <Hit_time_L2>")
#     sys.exit(1)

# Convert the command line arguments
block_size = 32 #int(sys.argv[1])
#L1_size = int(sys.argv[2])
#L1_assoc = int(sys.argv[3])
L2_size = 0 #int(sys.argv[4])
L2_assoc = 0    #int(sys.argv[5])
policy = 0  #sys.argv[6]
inclusion = 0   #sys.argv[7]
trace_file = "gcc_trace.txt"    # sys.argv[8]
#Hit_time_L1 = float(sys.argv[9])  # ns
#Hit_time_L2 = float(sys.argv[10]) # ns
miss_penalty = 100  # ns (fixed value as per your description)



# load the catcl Excel file
df = pd.read_excel("cacti_table.xls")

# Preview the DataFrame (optional)
#print(df.head())
# Create a nested dictionary based on Cache Size and Associativity :: dictionary format = { (Cache Size, Associativity): Access Time }
cache_access_data = {
    (row["Cache_Size"], row["Associativity"]): row["Access_Time"]
    for _, row in df.iterrows()
}
# Now you can access the access time using Cache Size and Associativity as keys
#cache_size = 1024  # Example cache size in bytes
#associativity = 4  # Example associativity
#access_time = cache_access_data.get((cache_size, associativity))
#print(f"Access Time for {cache_size} bytes and associativity {associativity}: {access_time} ns")file
print("excel file loaded")

associativities = [1, 2, 4, 8, 32]
cache_sizes = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576]

value_dict_X = {}
value_dict_miss = {}
value_dict_AAT = {}

sim_count = 0

for assoc in associativities:
    X = []
    miss_Y = []
    AAT_Y = []
    for L1_size in cache_sizes:
        L1_assoc = assoc
        if assoc == 8 and L1_size == 1024:
            continue


        Hit_time_L1 = cache_access_data.get((L1_size, L1_assoc))
        Hit_time_L2 = 0 #cache_access_data.get((L2_size, L2_assoc))

        if assoc == 32:  #Trigger for full assoc
            L1_assoc = L1_size/block_size
            Hit_time_L1 = cache_access_data.get((L1_size, " FA"))

        sim_count += 1
        input_string = str(sim_count) + "=> " + "./sim_cache.exe 32 " +  str(L1_size) + " " + str(L1_assoc) + " " + str(L2_size) + " " + str(L2_assoc) + " " + str(policy) + " " + str(inclusion) + " " + trace_file + " " + str(Hit_time_L1) + " " + str(Hit_time_L2) # + " " + str(miss_penalty)
        print(input_string)

        # input_args = [
        #     #"./sim_cache.exe", 
        #     "32", 
        #     str(L1_size), 
        #     str(L1_assoc), 
        #     str(L2_size), 
        #     str(L2_assoc), 
        #     str(policy), 
        #     str(inclusion), 
        #     trace_file, 
        #     str(Hit_time_L1), 
        #     str(Hit_time_L2)
        # ]

        process = subprocess.Popen(
            ["./sim_cache.exe", 
            "32", 
            str(L1_size), 
            str(L1_assoc), 
            str(L2_size), 
            str(L2_assoc), 
            str(policy), 
            str(inclusion), 
            trace_file, 
            str(Hit_time_L1), 
            str(Hit_time_L2)],
            stdout=subprocess.PIPE,  # Capture standard output
            stderr=subprocess.PIPE   # Capture standard error
        )

        # Read the output
        stdout, stderr = process.communicate()

        # Decode the byte output to string
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        result = stdout
        # print(result)

        data = [float(line) for line in str(result).splitlines()]
        # print("data is = ", data)

        # # Run the C++ program
        # result = subprocess.run([input_args], capture_output=True, text=True)
        # print("result = ", result)

        # # Parse the output
        # data = [float(line) for line in result.stdout.strip().splitlines()]



        miss_rate = data[0]
        AAT = data[1]


        X.append(math.log(L1_size, 2))
        miss_Y.append(miss_rate)
        AAT_Y.append(AAT)
    
    value_dict_X[assoc] = copy.deepcopy(X)
    value_dict_miss[assoc] = copy.deepcopy(miss_Y)
    value_dict_AAT[assoc] = copy.deepcopy(AAT_Y)



print (value_dict_X)
print (value_dict_miss)
print (value_dict_AAT)



# Plot each curve
for assoc in associativities:
    # plt.plot(value_dict_X[assoc], value_dict_miss[assoc], label="associativity = "+assoc if assoc!=32 else "full associativity")
    plt.plot(value_dict_X[assoc], value_dict_miss[assoc], label="associativity = " + str(assoc) if assoc != 32 else "full associativity")


plt.xlabel("L1 Cache Size")
plt.ylabel("miss rate")
plt.title("Graph 1")
# Show legend
plt.legend()
# Show the plot
plt.show()




# Plot each curve
for assoc in associativities:
    plt.plot(value_dict_X[assoc], value_dict_AAT[assoc], label="associativity = "+ str(assoc) if assoc!=32 else "full associativity")

plt.xlabel("L1 Cache Size")
plt.ylabel("AAT")
plt.title("Graph 2")
# Show legend
plt.legend()
# Show the plot
plt.show()






