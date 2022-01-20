from math import *
from shapely.geometry import Polygon

a = 6378137
e2 = 0.00669437999013
b2 = a ** 2 * (1 - e2)
ep2 = ((a ** 2) - b2) / b2
L0 = 19 * pi / 180

A0 = 1 - (e2 / 4) - ((3 * (e2 ** 2)) / 64) - ((5 * (e2 ** 3)) / 256)
A2 = (3 / 8) * (e2 + ((e2 ** 2) / 4) + ((15 * (e2 ** 3)) / 128))
A4 = (15 / 256) * (e2 ** 2 + ((3 * (e2 ** 3)) / 4))
A6 = (35 * (e2 ** 3)) / 3072


def wsp_gk(fi, lam, uklad):
    t = tan(fi)
    mi2 = ep2 * (cos(fi) ** 2)
    nr = 0

    if uklad == "1992":
        L0 = 19 * pi / 180
    else:
        if degrees(lam) < 16.5:
            L0 = radians(15)
            nr = 5
        elif 19.5 >= degrees(lam) > 16.5:
            L0 = radians(18)
            nr = 6
        elif 22.5 >= degrees(lam) > 19.5:
            L0 = radians(21)
            nr = 7
        else:
            L0 = radians(24)
            nr = 8

    l = lam - L0
    sigma = a * (A0 * fi - A2 * sin(2 * fi) + A4 * sin(4 * fi) - A6 * sin(6 * fi))
    N = a / (sqrt(1 - e2 * (sin(fi)) ** 2))

    xgk = sigma + ((l ** 2) / 2) * N * sin(fi) * cos(fi) * (
                1 + (l ** 2 / 12) * (cos(fi) ** 2) * (5 - t ** 2 + 9 * mi2 + 4 * (mi2 ** 2)) + ((l ** 4) / 360) *
                (cos(fi) ** 4) * (61 - 58 * (t ** 2) + (t ** 4) + 270 * mi2 - 330 * mi2 * (t ** 2)))

    ygk = l * N * cos(fi) * (1 + ((l ** 2) / 6) * (cos(fi) ** 2) * (1 - (t ** 2) + mi2) + ((l ** 4) / 120) *
                             (cos(fi) ** 4) * (5 - 18 * (t ** 2) + (t ** 4) + 14 * mi2 - 58 * mi2 * (t ** 2)))

    return xgk, ygk, nr

def wsp_1992(fi, lam):
    m01992 = 0.9993
    xgk, ygk, nr = wsp_gk(fi, lam, "1992")
    x1992 = m01992 * xgk - 5300000
    y1992 = m01992 * ygk + 500000

    return x1992, y1992

def wsp_2000(fi, lam):
    m02000 = 0.999923
    xgk, ygk, nr = wsp_gk(fi, lam, "2000")
    x2000 = m02000 * xgk
    y2000 = m02000 * ygk + nr * 1000000 + 500000

    return x2000, y2000, nr

def wsp_1992_to_gk(x, y):
    m01992 = 0.9993
    xgk = (x + 5300000) / m01992
    ygk = (y - 500000) / m01992

    return xgk, ygk

def wsp_2000_to_gk(x, y, nr):
    m02000 = 0.999923
    xgk = x / m02000
    ygk = (y - (nr * 1000000) - 500000) / m02000

    return xgk, ygk

def gk_to_fi_lam(x, y):
    fi0 = x / (a * A0)
    while True:
        sigma = a * (A0 * fi0 - A2 * sin(2 * fi0) + A4 * sin(4 * fi0) - A6 * sin(6 * fi0))
        fi1 = fi0 + (x - sigma) / a * A0
        if abs(fi1 - fi0) < radians(0.000001 / 3600):
            break
        else:
            fi0 = fi1

    t = tan(fi1)
    n2 = ep2 * (cos(fi1) ** 2)
    N = a / sqrt(1 - e2 * sin(fi1) ** 2)
    M = (a * (1 - e2)) / sqrt((1 - e2 * sin(radians(fi1)) ** 2) ** 3)

    fi = fi1 - (y ** 2 * t) / (2 * M * N) * (1 - (y ** 2) / (12 * N ** 2) *
        (5 + 3 * t ** 2 + n2 - 9 * n2 * t ** 2 - 4 * n2 ** 2) + (y ** 4) / (360 * N ** 4) * (61 + 90 * t ** 2 + 45 * t ** 4))
    lam = radians(19) + y / (N * cos(fi1)) * (1 - y ** 2 / (6 * N ** 2) * (1 + 2 * t ** 2 + n2) +
        y ** 4 / (120 * N ** 4) * (5 + 28 * t ** 2 + 24 * t ** 4 + 6 * n2 + 8 * n2 * t ** 2))

    return degrees(fi), degrees(lam)

def skala_gk(ygk, fi):
    fi = radians(fi)
    M = (a * (1 - e2)) / sqrt((1 - e2 * sin(radians(fi)) ** 2) ** 3)
    N = a / (sqrt(1 - e2 * (sin(fi)) ** 2))
    q = sqrt(M * N)

    m = 1 + ((ygk**2) / (2*q**2)) + ((ygk**2)/(24*q**4))
    Z = (1 - m) * 1000
    m2 = m**2
    Zha = (1 - m2) * 10000
    return m, Z, m2, Zha

def skala_2000(ygk, fi):
    fi = radians(fi)
    m = skala_gk(ygk, fi)[0]
    m2000 = 0.999923 * m
    Z = (1 - m2000) * 1000
    m2 = 0.999923**2 * m**2
    Zha = (1 - m2) * 10000
    return m2000, Z, m2, Zha

def skala_1992(ygk, fi):
    fi = radians(fi)
    m = skala_gk(ygk, fi)[0]
    m1992 = 0.9993 * m
    Z = (1 - m1992) * 1000
    m2 = 0.9993**2 * m**2
    Zha = (1 - m2) * 10000
    return m1992, Z, m2, Zha

A = [50.25, 20.75]
B = [50, 20.75]
C = [50.25, 21.25]
D = [50, 21.25]
E = [50.125, 21.0]
F = [50.12527, 21.00065]

fi = [A[0], B[0], C[0], D[0], E[0], F[0]]
lam = [A[1], B[1], C[1], D[1], E[1], F[1]]

x1992 = []
y1992 = []
x2000 = []
y2000 = []
y2000_nr = []
xgk = []
ygk = []

#transformacja 1992
for i in range(len(fi)):
    x, y = wsp_1992(radians(fi[i]), radians(lam[i]))
    x1992.append(x)
    y1992.append(y)

    x, y, num = wsp_2000(radians(fi[i]), radians(lam[i]))
    x2000.append(x)
    y2000.append(y)
    y2000_nr.append(num)

    x, y, nr = wsp_gk(radians(fi[i]), radians(lam[i]), '1992')
    xgk.append(x)
    ygk.append(y)

print("Współrzędne w 1992:")
print(x1992)
print(y1992)
print("Współrzędne w 2000:")
print(x2000)
print(y2000)
print("Współrzędne gk:")
print(xgk)
print(ygk)
print('\n')

#Pola
cords_1992 = [(x1992[0], y1992[0]), (x1992[1], y1992[1]), (x1992[3], y1992[3]), (x1992[2], y1992[2]), (x1992[0], y1992[0])]
cords_2000 = [(x2000[0], y2000[0]), (x2000[1], y2000[1]), (x2000[3], y2000[3]), (x2000[2], y2000[2]), (x2000[0], y2000[0])]
cords_gk = [(xgk[0], ygk[0]), (xgk[1], ygk[1]), (xgk[3], ygk[3]), (xgk[2], ygk[2]), (xgk[0], ygk[0])]

poly_1992 = Polygon(cords_1992)
poly_2000 = Polygon(cords_2000)
poly_gk = Polygon(cords_gk)

print("Pole w 1992", poly_1992.area/1000000)
print("Pole w 2000", poly_2000.area/1000000)
print("Pole gr", poly_gk.area/1000000)

new_x1992_gk = []
new_y1992_gk = []
new_x2000_gk = []
new_y2000_gk = []
new_xgk_filam = []

for i in range(len(fi)):
    x, y = wsp_1992_to_gk(x1992[i], y1992[i])
    new_x1992_gk.append(x)
    new_y1992_gk.append(y)

    x, y = wsp_2000_to_gk(x2000[i], y2000[i], y2000_nr[i])
    new_x2000_gk.append(x)
    new_y2000_gk.append(y)

    x, y = gk_to_fi_lam(xgk[i], ygk[i])
    new_xgk_filam.append(x)

fiGK92 = []
fiGK20 = []

for i in range(len(fi)):
    fii = gk_to_fi_lam(new_x1992_gk[i], new_y1992_gk[i])[0]
    fiGK92.append(fii)
    fii = gk_to_fi_lam(new_x2000_gk[i], new_y2000_gk[i])[0]
    fiGK20.append(fii)

print("Skala zniekształceń:")
s_GK = []
s_2000 = []
s_1992 = []

for i in range(len(fi)):
    s = skala_gk(ygk[i], new_xgk_filam[i])
    s_GK.append(s)

    s = skala_1992(new_y1992_gk[i], fiGK92[i])
    s_1992.append(s)

    s = skala_2000(new_y2000_gk[i], fiGK20[i])
    s_2000.append(s)

print(s_GK)
print(s_1992)
print(s_2000)





