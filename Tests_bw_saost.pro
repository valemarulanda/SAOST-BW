;----------------------
;Pupilles tx et rx uniformément éclairées
;----------------------
lambda =  1.55e-6
dsat =  1000e3
thetad =  10e-6
x = sqrt(8)
w0 =  lambda/(!pi*thetad)
dtx = w0*x

d1 = 26
polaire2, rt=d1/2., largeur=d1, masque=pup1
;tvwin,  pup1,  1,  zoom =  10

dx = (dtx/d1)/lambda

tt = total(pup1*dx^2)
print,  tt

d2 = 1.5*d1/dtx
polaire2, rt=d2/2., largeur=round(d2), masque=pup2
;tvwin,  pup2,  2

dx2 = (dtx/d1)/dsat

tt2 = total(pup2*dx2^2)
print,  tt2

print,  tt2*tt
;-------------------
;Pupille tx uniformément éclairée -> tache d'Airy en rx
;-------------------
lambda =  1.55e-6
dsat =  1000e3
 
thetad =  10e-6
x = sqrt(8)
w0 =  lambda/(!pi*thetad)
dtx = w0*x

lar =  1024

d1 = 16
polaire2, rt=d1/2., largeur=lar, masque=pup1
tvwin,  pup1,  1

dx = (dtx/d1)/lambda

tt = total(pup1*dx^2)
print,  tt

rx =  fftshift(pup1)
rx = rx/max(rx)
tvwin,  abs(rx),   2

p =  plot(abs(rx[lar/2, *]), xrange =  [0, lar])

mm = min(abs(rx[lar/2,lar/2:600]),ii)
print,  ii
rairy =  1.22*lambda*dsat/dtx

d2 = 1.5*ii/rairy
polaire2, rt=d2/2., largeur=1024, masque=pup2
tvwin,  pup2,  3

rec =  pup2*rx
tvwin,  abs(rec),  4

dx2 = (1.5/d2)/dsat

tt2 = total(rec*dx2^2)
print,  tt2

print,  tt2*tt

;------------------
;Couplage complexe pupille rx + mode SMF
;------------------

Drx =  1.5
dsat =  400e3
thetad =  10e-6
wp = thetad*dsat
t02 =  1 - exp(-2*((Drx/2)/wp)^2)
print,  t02

dim =  80
polaire2, rt=dim/2., largeur=dim, masque=a

ss = size(a)
npix = float(ss[1])

fwo = 1./(!pi*0.71)
wo = fwo*npix
R = distc(npix,cx=npix/2,cy=npix/2)
Mo = sqrt(2/(!pi*wo^2))*exp(-R^2/wo^2)    
    
Po = un(npix, npix, /double)

temp = scalar_prod(a, Mo, pup = Po)
omega = temp/sqrt(scalar_prod(A, A, pup = Po)) 

tflux = (abs2(omega))

print,  tflux

print,  tflux*t02

;------------------

lambda =  1.55e-6
dsat =  1000e3
 
thetad =  10e-6
x = sqrt(8)
w0 =  lambda/(!pi*thetad)
dtx = w0*x
drx = 1.5

lar =  6000.

d1 = 10
polaire2, rt=d1/2., largeur=lar, masque=p1
tvwin,  p1,  1, zoom = -1

du = (dtx/d1)/lambda
rx =  fftshift(p1)*du^2*lar^2
print,  max(rx)
;rx = rx*(drx/d2)^2;/max(rx)
tvwin,  abs(rx),   3, zoom = -1

rairy =  1.22*lambda*dsat/dtx

p =  plot(abs(rx[lar/2, *]),  thick =  2,  xrange =  [0, lar])

mm =  min(abs(rx[lar/2, lar/2:lar/2+lar/6]),jj)
print, jj

d2 = round(1.5*jj/rairy)
print,  d2

llim =  (lar-d2)/2
ulim =  llim+d2-1
rx2 = rx[llim:ulim, llim:ulim]
;tvwin, rx2, 4,  zoom =  -1

polaire2, rt=d2/2., largeur=d2, masque=p2
;tvwin,  p2,  5, zoom = -1

;tvwin,  p2*rx2, 6, zoom = -1

fwo = 1./(!pi*0.71)
wo = fwo*d2
R = distc(d2,cx=d2/2,cy=d2/2)
Mo = sqrt(2/(!pi*wo^2))*exp(-R^2/wo^2)    

;------------------

b = total(abs2(p1*du))
print,  b

da = (drx/d2)/dsat
c = total(abs2(p2*da))
print,  c

print, b*c ;pertes géométriques, approx: produit des surfaces
;------------------

d =  abs2(total(p2*Mo*da))
print,  d/c     ; couplage pupille + mode SMF
;------------------
print,  d*b ; pertes géométriques + couplage SMF (approx)
;------------------
k =  abs2(total(rx2*p2*Mo*da))
print,  k/b ; pertes géométriques + couplage SMF 
