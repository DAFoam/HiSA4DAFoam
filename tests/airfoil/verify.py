f = open("log.txt", "r")
lines = f.readlines()
f.close()

CD = []
for line in lines:
    if "Cd       :" in line:
        cols = line.split()
        CD.append(float(cols[2]))
print(CD)
CD_Final = CD[-1]

if abs(CD_Final - 0.27398) / 0.27398 > 1e-7:
    print("HiSA test failed!")
    exit(1)
else:
    print("HiSA test passed!")
