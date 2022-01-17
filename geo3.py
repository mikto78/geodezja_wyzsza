from math import *

#dane fi i lambda
fi_a = 50.25000
lam_a = 20.75000
fi_d = 50.00000
lam_d = 21.25000
a = 6378137
e2 = 0.00669437999013

#wsp punktu średniego położenia
fi_srednie = (fi_a + fi_d) / 2
lam_srednie = (lam_a + lam_d) / 2

#algorytm Vincenta
def vincent(fi_a, fi_b, lam_a, lam_b, a, e2):
    b = a * sqrt(1-e2)
    f = 1 - (b/a)
    delta_lam = radians(lam_b - lam_a)
    Ua = atan((1-f) * tan(fi_a))
    Ub = atan((1-f) * tan(fi_b))
    L = delta_lam
    L0 = 0
    sin_sigma, cos_sigma, sigma, sin_alfa, cos2_alfa, cos2_sigma = 0, 0, 0, 0, 0, 0

    while abs(L - L0) > 0.0000000000027:
        L0 = L
        sin_sigma = ((cos(Ub) * sin(L)) ** 2 + (cos(Ua) * sin(Ub) - sin(Ua) * cos(Ub) * cos(L)) ** 2) ** 0.5
        cos_sigma = sin(Ua) * sin(Ub) + cos(Ua) * cos(Ub) * cos(L)
        sigma = atan(sin_sigma/cos_sigma)

        sin_alfa = (cos(Ua) * cos(Ub) * sin(L)) / sin_sigma
        cos2_alfa = 1 - sin_alfa ** 2
        cos2_sigma = cos_sigma - ((2 * sin(Ua) * sin(Ub)) / cos2_alfa)

        C = (f / 16) * cos2_alfa * (4 + f * (4 - 3 * cos2_alfa))
        L = delta_lam + (1 - C) * f * sin_alfa * (sigma + C * sin_sigma * (cos2_sigma + C * cos_sigma * ((-1) + 2 * cos2_sigma ** 2)))

    u2 = ((a ** 2 - b ** 2) / b ** 2) * cos2_alfa
    A = 1 + (u2/16384) * (4096 + u2 * ((-768) + u2 * (320 - 175 * u2)))
    B = (u2 / 1024) * (256 + u2 * ((-128) + u2 * (74 - 47 * u2)))

    delta_sigma = B * sin_sigma * (cos2_sigma + (1 / 4) * B * (cos_sigma * ((-1) + 2 * cos2_sigma ** 2) - (1 / 6) * B * cos2_sigma * ((-3) + 4 * sin_sigma ** 2) * ((-3) + 4 * cos2_sigma ** 2)))
    s_ab = b * A * (sigma - delta_sigma)
    y_ab = cos(Ub) * sin(L)
    x_ab = cos(Ua) * sin(Ub) - sin(Ua) * cos(Ub) * cos(L)
    y_ba = cos(Ua) * sin(L)
    x_ba = -sin(Ua) * cos(Ub) - cos(Ua) * sin(Ub) * cos(L)

#wyrówanie azymutów
    if y_ab > 0 and x_ab > 0:
        A_ab = atan(y_ab / x_ab)
    elif y_ab > 0 and x_ab < 0:
        A_ab = atan(y_ab / x_ab) + pi
    elif y_ab < 0 and x_ab < 0:
        A_ab = atan(y_ab / x_ab) + pi
    elif y_ab < 0 and x_ab > 0:
        A_ab = atan(y_ab / x_ab) + 2 * pi

    if y_ba > 0 and x_ba > 0:
        A_ba = atan(y_ba / x_ba)
    elif y_ba > 0 and x_ba < 0:
        A_ba = atan(y_ba / x_ba) + 2 * pi
    elif y_ba < 0 and x_ba < 0:
        A_ba = atan(y_ba / x_ba) + 2 * pi
    elif y_ba < 0 and x_ba > 0:
        A_ba = atan(y_ba / x_ba) + 3 * pi

    return s_ab, degrees(A_ab), degrees(A_ba)


#algorytm Kivoji
def kivoji(fi, lam, s, az):
    n = int(s / 1000)
    ds = s / n
    for i in range(n):
        M = (a*(1 - e2)) / sqrt((1 - e2*sin(fi)**2)**3)
        N = a / sqrt(1 - e2*sin(fi)**2)

        dfi = (ds * cos(az)) / M
        daz = (sin(az) * tan(fi) * ds) / N

        fisr = fi + (dfi / 2)
        azsr = az + (daz / 2)

        M = (a * (1 - e2)) / sqrt((1 - e2 * sin(fisr) ** 2) ** 3)
        N = a / sqrt(1 - e2 * sin(fisr) ** 2)

        dfi = (ds * cos(azsr)) / M
        dlambda = (ds * sin(azsr)) / (N * cos(fisr))
        daz = (sin(azsr) * tan(fisr) * ds) / N

        fi = fi + dfi
        lam = lam + dlambda
        az = az + daz

    return degrees(fi), degrees(lam), degrees(az)

#obliczenia
s, A_ad, A_da = vincent(radians(fi_a), radians(fi_d), lam_a, lam_d, a, e2)
fi_sr, lam_sr, Az_ad = kivoji(radians(fi_a), radians(lam_a), s/2, radians(A_ad))
delta_s, Az_ad, Az_da = vincent(radians(fi_srednie), radians(fi_sr), lam_srednie, lam_sr, a, e2)
delta_s, Az_da, Az_aaa = vincent(radians(fi_sr), radians(fi_srednie), lam_sr, lam_srednie, a, e2)

#funkcja na obliczenie pola
def pole(f_a, l_a, f_b, l_b):
    e = e2 ** 0.5
    b = (a * sqrt(1 - e2)) ** 2
    A = sin(f_a)/(1 - e2*(sin(f_a)**2)) + log((1+e*sin(f_a))/(1-e*sin(f_a)))/(2*e)
    B = sin(f_b)/(1 - e2*(sin(f_b)**2)) + log((1+e*sin(f_b))/(1-e*sin(f_b)))/(2*e)
    P = b*(l_b - l_a)/2*(A - B)
    return P

#zamiana stopni na stopnie-minuty-sekundy
def st_min_sek(stopnie):
    deg = int(stopnie)
    minutes = int((stopnie - deg)*60)
    seconds = round(((stopnie - deg - minutes/60) * 3600), 5)
    deg = f"{deg:02d}"
    minutes = f"{minutes:02d}"
    int_sec = f"{int(seconds):02d}"
    float_seconds = str(round(seconds-int(seconds), 5))[1:]

    return f'{deg}°{minutes}\'{int_sec}{float_seconds}"'


#wypisanie obliczonych wartości
print("1. Punkt średniej szerokości:" + str(st_min_sek(fi_srednie)) + " " + str(st_min_sek(lam_srednie)))
print("2. Azymut A -> D:   " + str(st_min_sek(A_ad)) + '\n' + "Punkt środkowy: " + str(st_min_sek(fi_sr)) + " " + str(st_min_sek(lam_sr)) )
print("3. Rożnica odległości: " + str(round(delta_s, 3)) + " metrów")
print("4. Azymut w tych punktach: " + str(st_min_sek(Az_ad)) + " " + str(st_min_sek(Az_da)))
print("5. Pole powierzchni czworokąta: " + str(round(pole(radians(fi_a), radians(lam_a), radians(fi_d), radians(lam_d)), 6)) + " m^2")


