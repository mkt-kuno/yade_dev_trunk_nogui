

#combine all 4 elements of the Quaternion into one structure

import numpy as np 

print("Quaternion creation started")

#store 4 Quat elements in arrays

input0 = inputs[0]  

for i in input0.PointData.keys():#KB: I just explicitely copy all the data

    output.PointData.append(input0.PointData[i], i)


dataArray = input0.PointData["Quatw"]

dataArray1 = input0.PointData["Quatx"]

dataArray2 = input0.PointData["Quaty"]

dataArray3 = input0.PointData["Quatz"]



#######

#stack the arrays together to create 4 element Quat

All_Data2=np.vstack((dataArray,dataArray1))

All_Data2=np.vstack((All_Data2,dataArray2))

All_Data2=np.vstack((All_Data2,dataArray3))

#transpose the array to create the form (number point X 4)

All_Data2=np.transpose(All_Data2)

print("Quaternion of sphere")

print(All_Data2)

output.PointData.append(All_Data2, "Quaternion")

########## 

print("done")
