import matplotlib.pyplot as plt
import csv

brute = []
ddot = []
cblas = []
with open("speed_test_results_brute.txt", "r") as file:
    reader = csv.reader(file,delimiter=",")
    for row in reader:
        brute.append(float(row[0])*10e-9)
with open("speed_test_results_ddot.txt", "r") as file:
    reader = csv.reader(file,delimiter=",")
    for row in reader:
        ddot.append(float(row[0])*10e-9)
with open("speed_test_results_cblas_direct.txt", "r") as file:
    reader = csv.reader(file,delimiter=",")
    for row in reader:
        cblas.append(float(row[0])*10e-9)
ns = [100000 + 10000*i for i in range(1,len(ddot)+1)]

fig, axs = plt.subplots()
axs.plot(ns,brute,color='k',linestyle='--',label='brute')
axs.plot(ns,ddot,color='r',linestyle='--',label='ddot')
#axs.plot(ns,cblas,color='b',linestyle='--',label='cblas')
axs.set_yscale('log')
axs.set_ylabel(r'sec $(\log)$')
axs.set_xlabel('dimension')
axs.set_title(r'sec $(\log)$ vs. dimension')
plt.suptitle('Brute vs. DDOT Speed Test')
plt.legend()
plt.show()
