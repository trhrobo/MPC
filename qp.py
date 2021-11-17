import numpy as np

B = np.array([[0], [0]])

K = []

temp_val = [B, B, B, B]

temp0 = np.array([[1], [1]])
temp1 = np.array([[2], [2]])

temp_val[0] = temp0
temp_val[1] = temp0
temp_val[2] = temp1
temp_val[3] = temp1

val = np.array([[temp_val[0], temp_val[1]]])

print(val)

K.append(val)

val = np.array([[temp_val[2], temp_val[3]]])
K.append(val)

print("----")
print(K)