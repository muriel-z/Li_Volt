# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 15:26:26 2022

@author: Muri
"""


import numpy as np 
import matplotlib.pyplot as plt


# Inicialización de variables

dt = 1e-6
n0 = 77     #comienzo con este numero de li0 depositados sobre la sup
n0max = 600 #numero maximo de particulas li0
rli0 = 1.67e-10  #radio atomico del li0
rlim = 1.2e-10   #radio atomico del li+
nm = 50          #Este numero debe mantenerse a lo largo de la ev temporal
D = 1.4e-14      #coef de dif del Li+ en el electrolito
tita = 0
gx = 0
gy = 0
x=0
y=0
Long = 16.7e-9
datt = 1.3*rli0/Long     #que es rli0/(pi/4)
qu = np.sqrt(2*D*dt)/Long   #desplazamiento medio debido a la difusion

alfa = 0
mod_li0 = 0

Rli_mx = 0
Rli_my = 0
Rli_0x = 0
Rli_0y = 0
modd = 0

u = 5.6e-13

#Variables de la malla de voltaje y campo
N_iter = 1000
nmalla = 100
dx = Long/(nmalla*Long)
dy = Long/(nmalla*Long)

#Defino y calculo el voltaje V0
Vm = 0.85 #Volts
V0 = np.zeros([nmalla+1,nmalla+1])
V =  np.zeros([nmalla+1,nmalla+1])

V0[:,0] = Vm #Condición de borde

for k in range(N_iter): 
    for i in range(1,nmalla,1):
        for j in range(1,nmalla,1):
            V0[i,j] = (V0[i+1,j]+V0[i-1,j]+V0[i,j+1]+V0[i,j-1])/4
            
#Defino y Calculo el campo E0
E0_x = np.zeros([nmalla+1,nmalla+1])
E0_y = np.zeros([nmalla+1,nmalla+1])
E_x = np.zeros([nmalla+1,nmalla+1])
E_y = np.zeros([nmalla+1,nmalla+1])

for j in range(0,nmalla+1,1):
    for i in range(1,nmalla,1):
        E0_x[i,j] = -(V0[i+1,j]-V0[i-1,j])/dx
     
for j in range(0,nmalla+1,1): #Defino el campo en x respetando las PBC 
    E0_x[0,j] = -(V0[1,j]-V0[nmalla,j])/dx
    E0_x[nmalla,j] = - (V0[0,j]-V0[nmalla-1,j])/dx
    
for j in range(1,nmalla,1):
    for i in range(0,nmalla+1,1):
        E0_y[i,0] = -V0[i,1]/dy
        E0_y[i,j] = -(V0[i,j+1]-V0[i,j-1])/dy
        E0_y[i,nmalla] = V0[i,nmalla-1]/dy
        

#Inicialización de vectores 
ex = np.zeros(nm)
ey = np.zeros(nm)
ex_0 = np.zeros(nm)
ey_0 = np.zeros(nm)
lix_d = np.zeros(n0)
liy_d = np.zeros(n0)
lix_aux = np.zeros(n0max)
liy_aux = np.zeros(n0max)

#Creación de la monocapa de Li0 depositado

for i in range(n0):
    lix_d[i] = 1.3*rli0*i/Long
    liy_d[i] = 0/Long

#Creación de las posiciones iniciales de los iones

dx = rli0/(2*Long)
dy = rli0/(2*Long)

for i in range(nm):
    ex_0[i] = np.random.rand()
    ey_0[i] = np.random.rand()
    ex_0[i] = ex_0[i] % 1
    if (ey_0[i]<0):
        ey_0[i] = abs(ey_0[i])
    elif (ey_0[i]>1):
        ey_0[i] = 1-(ey_0[i]-1)
    
#Guardo las coordenadas del sistema inicial
ex_init = ex_0
ey_init = ey_0
 
#Calculo los vectores rx y ry

rx = np.zeros(nm)
ry = np.zeros(nm)

indices = [] #Me guardo los indices que entran

#Bueno empiezo a escribir a ver kionda
Rcut = 0.007 #1/100 (escala de las separaciones)/2*raiz(2)
Rcut2 = Rcut*Rcut

for i in range(nm):
    for j in range(nmalla+1):
        for k in range(nmalla+1):
            d1 = ex_0[i] - j*0.01
            ka = k*0.01 + (1-2*k*0.01)
            d2 = ey_0[i] - ka
            dist2 = d1*d1 + d2*d2
            div = 1
            if dist2<=Rcut2:
                indices.append((j,k))
                div = div+1
                rx[i] = u*E0_x[j,k]*dt/(Long*div)
                ry[i] = -u*E0_y[j,k]*dt/(Long*div)
                

#mu+*E*dt/Long es el desplazamiento debido al campo electrico

        
#Comienza el loop de evolución temporal


m = n0      #Número inicial de Li0
nt = 50  #Número de pasos temporales


#Inicialmente los vectores de Li0 son iguales al litio depositado
lix_0 = lix_d
liy_0 = liy_d


lix_aux = np.zeros(n0max)
liy_aux = np.zeros(n0max)

porcent = 9 #Porcentaje inicial para el contador

exys = []
pqs = []
list_p = []
list_q = []

for i in range(nt):
    break_out_of_i = False
    for j in range(nm):
        break_out_of_j = False
        tita = 2*np.pi*np.random.rand()
        #Definición del vector unitario g
        gx = np.cos(tita)
        gy = np.sin(tita)
        #Evolución de la posición de los iones
        ex[j] = ex_0[j] + qu*gx + rx[j]
        ey[j] = ey_0[j] + qu*gy + ry[j]
        #Meto las PBC en la caja de tamaño 1x1 (normalizada)
        ex[j] = ex[j] % 1
        if (ey[j]<0):
            ey[j] = abs(ey[j])
        elif (ey[j]>1):
            ey[j] = 1-(ey[j]-1)
    
        #Definición de la condición Li+-->Li0
        for k in range(m):
            if (m<=n0):
                distx = ex[j] - lix_d[k]
                disty = ey[j] - liy_d[k]
                dist = np.sqrt(distx*distx + disty*disty)
            elif (m>n0):
                distx = ex[j] - lix_0[k]
                disty = ey[j] - liy_0[k]
                dist = np.sqrt(distx*distx + disty*disty)
                
            if (dist<datt):
                if (m<=n0):
                    modd = np.sqrt((ex[j] - lix_d[k])*(ex[j] - lix_d[k]) + (ey[j] - liy_d[k])*(ey[j] - liy_d[k]))
                    exs = (ex[j]-lix_d[k])*datt/modd + lix_d[k]
                    eys = (ey[j]-liy_d[k])*datt/modd + liy_d[k]
                elif (m>n0):
                    modd = np.sqrt((ex[j] - lix_0[k])*(ex[j] - lix_0[k]) + (ey[j] - liy_0[k])*(ey[j] - liy_0[k]))
                    exs = (ex[j]-lix_0[k])*datt/modd + lix_0[k]
                    eys = (ey[j]-liy_0[k])*datt/modd + liy_0[k]
                
                exys.append([exs,eys])
                m = m+1
                if (m==n0max):
                    break_out_of_i = True
                    break_out_of_j = True
                    print('Se alcanzó la cantidad máxima de Li0')
                    break
                lix_0 = np.zeros(m)
                liy_0 = np.zeros(m)
                for l in range(n0):
                    lix_aux[l] = lix_d[l]
                    liy_aux[l] = liy_d[l]
                    
                lix_aux[m] = exs
                liy_aux[m] = eys
                
                for l in range(m):
                    lix_0[l] = lix_aux[l]
                    liy_0[l] = liy_aux[l]
                    
                #PBC
                for l in range(m):
                    lix_0[l] = lix_0[l] % 1
                    if (liy_0[l]<0):
                        liy_0[l] = abs(liy_0[l])
                    elif (liy_0[l]>1):
                        liy_0[l] = 1-(liy_0[l]-1)
                #Calculo el nuevo voltaje
                #Primero paso las coordenadas a indices de la malla
                p = int(exs*100)
                q0 = int(eys*100)
                q = q0 + (100-2*q0)
                list_p.append(p)
                list_q.append(q)
                pqs.append((p,q))
                V[:,:] = V0[:,:] #inicialmente
                V[:,0] = Vm
                V[:,100] = 0
                V[0,:] = 0
                V[100,:] = 0
                V[p,q] = 0 #Lugar del nuevo Li0
                for kk in range(N_iter): 
                    for ii in range(1,nmalla,1):
                        for jj in range(1,nmalla,1):
                            for ll in range(np.size(list_p)):
                                for mm in range(np.size(list_q)):
                                    site_p = list_p[ll]
                                    site_q = list_q[mm]
                                    V[site_p,site_q] = 0
                                    V[ii,jj] = (V[ii+1,jj]+V[ii-1,jj]+V[ii,jj+1]+V[ii,jj-1])/4
                        for ll in range(np.size(list_p)):
                            for mm in range(np.size(list_q)):
                                site_p = list_p[ll]
                                site_q = list_q[mm]
                                V[site_p,site_q] = 0
                            
                V0[:,:] = V[:,:] #Para la prox iteracion
                #Calculamos el campo electrico a partir de este potencial
                for jj in range(0,nmalla+1,1):
                    for ii in range(1,nmalla,1):
                        E_x[ii,jj] = -(V0[ii+1,jj]-V0[ii-1,jj])/dx
                     
                for jj in range(0,nmalla+1,1): #Defino el campo en x respetando las PBC 
                    E_x[0,j] = -(V0[1,j]-V0[nmalla,jj])/dx
                    E_x[nmalla,jj] = - (V0[0,jj]-V0[nmalla-1,jj])/dx
                    
                for jj in range(1,nmalla,1):
                    for ii in range(0,nmalla+1,1):
                        E_y[ii,0] = -V0[ii,1]/dy
                        E_y[ii,jj] = -(V0[ii,jj+1]-V0[ii,jj-1])/dy
                        E_y[ii,nmalla] = V0[ii,nmalla-1]/dy
                #Calculo los vectores rx y ry
                for ii in range(nm):
                    for jj in range(nmalla+1):
                        for kk in range(nmalla+1):
                            d1 = ex_0[ii] - jj*0.01
                            ka = kk*0.01 + (1-2*kk*0.01)
                            d2 = ey_0[ii] - ka
                            dist2 = d1*d1 + d2*d2
                            div = 1
                            if dist2<=Rcut2:
                                indices.append((jj,kk))
                                div = div+1
                                rx[ii] = u*E_x[jj,kk]*dt/(Long*div)
                                ry[ii] = -u*E_y[jj,kk]*dt/(Long*div)
                #Repongo el ion
                ex[j] = np.random.rand()
                ey[j] = np.random.rand()
                #Las PBC
                ex[j] = ex[j] % 1
                if (ey[j]<0):
                    ey[j] = abs(ey[j])
                elif (ey[j]>1):
                    ey[j] = 1-(ey[j]-1)
        if break_out_of_j: break
                
    t = (i+1)*dt
    ex_0 = ex
    ey_0 = ey
    #Agrego el pocentaje de corrida del programa
    ip = int((i/nt)*100)
    if (ip > porcent):
        if (ip % 10 == 0):
            print(ip,'%')
        porcent = ip
    if break_out_of_i: break
    
print('Cantidad de Li0 =', m)
              
#%%  

#En este bloque voy a guardar los sitios de la malla que corresponden a equipotenciales.

equip_sites_x1 = []
equip_sites_y1= []
equip_sites_x2 = []
equip_sites_y2= []
equip_sites_x3 = []
equip_sites_y3= []
equip_sites_x4 = []
equip_sites_y4= []
equip_sites_x5 = []
equip_sites_y5= []
volt = 0.84
vol = [0.84,0.67,0.4,0.2,0]


for i in range(nmalla+1):
    for j in range(nmalla+1):
        if 0.85>=V0[i,j]>=0.84:
            equip_sites_x1.append(i/100)
            equip_sites_y1.append((j+(100-2*j))/100)
            
for i in range(nmalla+1):
    for j in range(nmalla+1):
        if 0.67>=V0[i,j]>=0.66:
            equip_sites_x2.append(i/100)
            equip_sites_y2.append((j+(100-2*j))/100)
            
for i in range(nmalla+1):
    for j in range(nmalla+1):
        if 0.4>=V0[i,j]>=0.39:
            equip_sites_x3.append(i/100)
            equip_sites_y3.append((j+(100-2*j))/100)
            
for i in range(nmalla+1):
    for j in range(nmalla+1):
        if 0.2>=V0[i,j]>=0.19:
            equip_sites_x4.append(i/100)
            equip_sites_y4.append((j+(100-2*j))/100)


for i in range(nmalla+1):
    for j in range(nmalla+1):
        if 0.001>=V0[i,j]>=0.0008:
            equip_sites_x5.append(i/100)
            equip_sites_y5.append((j+(100-2*j))/100)


plt.figure(40)
ax = plt.subplot(111)
ax.plot(equip_sites_x1,equip_sites_y1,'-',linewidth=2, label=' Equipotencial en 0.85 V')
ax.plot(equip_sites_x2,equip_sites_y2,'-',linewidth=2, label=' Equipotencial en 0.67 V')
ax.plot(equip_sites_x3,equip_sites_y3,'-',linewidth=2, label=' Equipotencial en 0.4 V')
ax.plot(equip_sites_x4,equip_sites_y4,'-',linewidth=2, label=' Equipotencial en 0.2 V')
ax.plot(equip_sites_x5,equip_sites_y5,'-',linewidth=2, label=' Equipotencial en 0.0 V')
plt.xlabel('sitios de la malla')
plt.ylabel(' sitios de la malla ')
ax.legend()
ax.set_aspect('equal', adjustable='box')







#%%
#Gráfico del sistema inicial

plt.figure(10)
x = ex_init
f = ey_init
y = lix_d
g = liy_d
ax = plt.subplot(111)
ax.plot(x,f,linestyle='', marker='.', markersize=7,color='#1f77b4', label=' Li$^+$')
ax.plot(y,g,linestyle='', marker='.', markersize=7,color='#d62728', label=' Li$^0$')
plt.ylim([-0.05,1.05])
plt.xlabel('x')
plt.ylabel('y')
plt.title(' init ')
ax.legend()      
ax.set_aspect('equal', adjustable='box')

#Gráfico del sistema final
        
plt.figure(20)
x = ex
f = ey
y = lix_0
g = liy_0
ax = plt.subplot(111)
ax.plot(x,f,linestyle='', marker='.', markersize=7,color='#1f77b4', label=' Li$^+$')
ax.plot(y,g,linestyle='', marker='.', markersize=7,color='#d62728', label=' Li$^0$')
plt.ylim([-0.05,1.05])
plt.xlabel('x')
plt.ylabel('y')
plt.title(' ')
ax.legend()
ax.set_aspect('equal', adjustable='box')


#Intento graficar el sistema final mas los euipotenciales
plt.figure(30)
plt.pcolormesh(V0[:,:])
plt.xlabel('')
plt.ylabel('')


#plt.figure(40)
#ax.plot_surface(V0)

plt.show()
