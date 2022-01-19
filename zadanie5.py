from math import *
from numpy import *

#Punkty
A = [50.25, 20.75]
B = [50, 20.75]
C = [50.25, 21.25]
D = [50, 21.25]
E = [50.125, 21.0]
F = [50.12527, 21.00065]

fi = [A[0], B[0], C[0], D[0], E[0], F[0]]
lam = [A[1], B[1], C[1], D[1], E[1], F[1]]

#GRS 80
a = 6378137
e2 = 0.0066943800290

#elipsoida krasowskiego
ak = 6378245
e2k = 0.0066934215520
x0 = -33.4297
y0 = 146.5746
z0 = 76.2865
m = 0.8407728 / 1000000
ex = radians(-0.35867/3600)
ey = radians(-0.05283/3600)
ez = radians(0.84354/3600)

def geo_xyz(fi, lam, a, e2, H):
    fi = radians(fi)
    lam = radians(lam)
    N = a / sqrt(1 - e2 * (sin(fi) ** 2))
    x = (N + H) * cos(fi) * cos(lam)
    y = (N + H) * cos(fi) * sin(lam)
    z = (N * (1 - e2) + H) * sin(fi)

    return x, y, z

def Hirvonen(x, y, z, a, e2):
    r = (x ** 2 + y ** 2)**0.5
    fi1 = atan((z / r) / (1 - e2))
    N = a / sqrt((1 - e2 * (sin(fi1) ** 2)))
    h = (r / cos(fi1)) - N
    fi2 = atan((z / r) / (1 - e2 * (N / (N + h))))

    while abs(fi2 - fi1) >= radians(0.00005 / 3600):
        fi1 = fi2
        N = a / sqrt((1 - e2 * (sin(fi1) ** 2)))
        h = (r / cos(fi1)) - N
        fi2 = atan((z / r) / (1 - e2 * (N / (N + h))))

    N = a / sqrt((1 - e2 * (sin(fi2) ** 2)))
    h = (r / cos(fi2)) - N
    lam = atan(y/x)

    return st_min_sek(degrees(fi2)), st_min_sek(degrees(lam)), h

def st_min_sek(stopnie):
    deg = int(stopnie)
    minutes = int((stopnie - deg)*60)
    seconds = round(((stopnie - deg - minutes/60) * 3600), 5)
    deg = f"{deg:02d}"
    minutes = f"{minutes:02d}"
    int_sec = f"{int(seconds):02d}"
    float_seconds = str(round(seconds-int(seconds), 5))[1:]

    return f'{deg}°{minutes}\'{int_sec}{float_seconds}"'

def transformacja(x, y, z):
    macierz = array([[x0], [y0], [z0]])
    macierz_wsp = array([[x], [y], [z]])
    macierz_obrotu = array([[m, ez, -ey], [-ez, m, ex], [ey, -ex, m]])
    wynik = macierz_wsp + macierz_obrotu.dot(macierz_wsp) + macierz

    return transpose(wynik)


x_GRS80 = []
y_GRS80 = []
z_GRS80 = []

for i in range(len(fi)):
    x, y, z = geo_xyz(fi[i], lam[i], a, e2, 0)
    x_GRS80.append(x)
    y_GRS80.append(y)
    z_GRS80.append(z)

print("Współrzędne xyz w GRS80")
print(x_GRS80)
print(y_GRS80)
print(z_GRS80)
print("--------------------------------------------------")
print("Współrzędne xyz w Krasowskim")
x_Krasowski = []
y_Krasowski = []
z_Krasowski = []

for i in range(len(fi)):
    x = transformacja(x_GRS80[i], y_GRS80[i], z_GRS80[i])
    x_Krasowski.append(x[0][0])
    y_Krasowski.append(x[0][1])
    z_Krasowski.append(x[0][2])

print(x_Krasowski)
print(y_Krasowski)
print(z_Krasowski)

print("Współrzędne w elipsoidzie Krasowskiego")
for i in range(len(fi)):
    x = Hirvonen(x_Krasowski[i], y_Krasowski[i], z_Krasowski[i], ak, e2k)
    print(x)





