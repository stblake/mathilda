# 4.3.7 DSolve -- linear constant-coefficient equations
DSolve[y''[x] + 4 y'[x] + 5 y[x] == 0, y[x], x]
DSolve[y''[x] - y[x] == E^x, y[x], x]
DSolve[y''[x] == 7, y, x]
DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi/2] == 1}, y[x], x]
DSolve[{y''[x] + y[x] == 0, y[0] == 1, y[Pi] == 1}, y, x]
