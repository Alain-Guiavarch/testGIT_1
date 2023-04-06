'''
Module permettant de calculer les sollicitations météorologiques : 
Température extérieur et rayonnement incident dans le plan des capteurs.

-----Paramètres-----
L_loc :  Longitude
Phi : Latitude
Beta : Inclinaison du capteur
Gamma : Azimut du capteur

-----Entrées-----
Ib : Irradiance directe horizontale
Id = Irradiance diffuse horizontale

-----Sorties-----
Text : température extérieure
It : Irradiance dans le plan du capteur


'''

from pathlib import Path
import math
from math import sin
from math import cos

wd = Path.cwd()

meteoFile = wd/"FRA_AC_Agen-La.Garenne.AP.075240_TMYx.2004-2018/FRA_AC_Agen-La.Garenne.AP.075240_TMYx.2004-2018.epw"
drivingForces = wd/"drivingForces.csv"

en_tete = 7
with open(meteoFile, "r") as f:
    data = f.read().splitlines()
    del data[0:en_tete+1]

#le format avec identification des colonnes est donné ici https://designbuilder.co.uk/cahelp/Content/EnergyPlusWeatherFileFormat.htm
Text, Ib, Id =[], [], []
col_Text, col_Ib, col_Id = 6, 13, 15
for line in data:
    Text.append(line.split(",")[col_Text])
    Ib.append(line.split(",")[col_Ib])
    Id.append(line.split(",")[col_Id])


with open(drivingForces, "w") as f:
    for i, Temp in enumerate(Text):
        ligne = Temp + ';' +  Ib[i] + ';' + Id[i] + '\n'
        f.write(ligne)

# Paramètres
L_loc =  - 5 #site étudié
L_st = -15 # longitude standard - 15° pour la France
phi = 45*math.pi/180
beta = 45*math.pi/180
gamma = 0*math.pi/180


Nh = 8760
temp = wd/"temp.csv"
It=[]
with open(temp, "w") as f:
    f.write("temps;jour;Ete;heure;E;ST;TSV;delta;omega;cosTeta;cosTetaZ;Ib;It\n")
    for h in range(Nh):
        n=h // 24 +1 # numéro du jour
        B=(n-1)*360/365*math.pi/180
        E = 229.2*(0.000075+0.001868*cos(B)-0.032077*sin(B)-0.014615*cos(2*B)-0.04089*sin(2*B))# correction du temps en minute
        heure = h % 24
        ST = heure # temps standard
        Ete=0
        if n == 31+28+25+1 and heure >=3:
            Ete = 1
        elif n>(31+28+25)+1 and n<(31+28+31+30+31+30+31+31+30+29):
            Ete = 1
        Ete = 0
        TSV = ST-Ete+(4*(L_st-L_loc)+E)/60 #temps solaire
        if heure == 0:
            TSV = 24+TSV
        delta = 23.45*sin(360*(284+n)/365*math.pi/180)
        omega = 15*(12-TSV)*math.pi/180
        cosTeta = (  sin(delta)*sin(phi)*cos(beta)
                    -sin(delta)*cos(phi)*sin(beta)*cos(gamma) 
                    +cos(delta)*cos(phi)*cos(beta)*cos(omega) 
                    +cos(delta)*sin(phi)*sin(beta)*cos(gamma)*cos(omega)
                    +cos(delta)*sin(beta)*sin(gamma)*sin(omega)           )
        cosTetaZ=cos(delta)*cos(phi)*cos(omega)+sin(delta)*sin(phi)
        Rb = cosTeta / cosTetaZ
        It.append(str(float(Ib[h])*Rb))
        f.write(f"{h};{n};{Ete};{heure};{E};{ST};{TSV};{delta};{omega*180/math.pi};{cosTeta};{cosTetaZ};{Ib[h]};{It[h]}\n")


'''
    for i, Temp in enumerate(Text):
        ligne = Temp + ',' +  Ib[i] + ',' + Id[i] + '\n'
        f.write(ligne)
        '''