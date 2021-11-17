import numpy as np

A = np.array([[1, 1], [0, 1]])
B = np.array([[0.5], [1]])

ans = np.dot(A, B)

print(ans)

ans *= 2

print(ans)