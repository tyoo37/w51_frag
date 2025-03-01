# python forGETSF0.py 
import sys
import os
sys.path.insert(0, '/Users/flouvet/python')
from mymodules_python3 import *
from astropy.io import fits
import math
from scipy import special


######################################
# data from .dat file

freq = 'B6'

au = 1.4959784260968E+11
msun = 1.9891E30 #kg
arcsec = 4.8481368110954E-06
h = 6.62607004e-34 
k = 1.3806488E-23
beta = 1.5 # indice spectral
T = 20.
maxsource = 3. # les sources doivent etre plus petites que maxsource fois le beam
cormax = 4.    # correction maximale acceptee entre la taille a 1mm et la taille a 3mm
temperror = 0.25 # pourcentage d'erreur sur la temperature
printtable     = True
tempcorrection = True	

# liste de base (nom, distance [kpc], frequence de ref B3 [GHz], frequence de ref B6 [GHz], coeurs freefree)
sources =     [['G327',2.5,101.776,229.507,[0]],                            ['G328',2.5,101.500,227.575,[0]],   ['G337',2.7,101.602,227.503,[0]],        ['G338',3.9,100.882,229.226,[0]],     ['W43-MM1',5.5,99.759,229.680,[0]],        ['W43-MM2',5.5,101.017,227.597,[0]],                       ['W43-MM3',5.5,100.911,228.931,[0]],        ['W51-E',5.4,101.426,228.918,[0]],                                                                       ['G351',2.0,100.228,227.991,[0]],                     ['G353',2.0,100.547,229.431,[0]],                                                                                                          ['G008',3.4,100.526,228.732,[0]],          ['G010',4.95,100.704,229.268,[0]],                                                               ['G012',2.4,100.680,229.080,[0]],                                                                                                             ['G333',4.2,100.756,229.062,[0]],                                                                                                                                                                                                                     ['W51-IRS2',5.4,101.263,228.530,[0]]                                                                               ]
#
# j'ai enleve toutes les sources toutes les sources vertes (alpha<2) et toutes les sources dans du halpha 
sourcesoeil = [['G327',2.5,101.776,229.507,[0]],                            ['G328',2.5,101.500,227.575,[6]],   ['G337',2.7,101.602,227.503,[0]],        ['G338',3.9,100.882,229.226,[0]],     ['W43-MM1',5.5,99.759,229.680,[0]],        ['W43-MM2',5.5,101.017,227.597,[0]],                       ['W43-MM3',5.5,100.911,228.931,[33]],        ['W51-E',5.4,101.426,228.918,[2, 16, 18, 20, 21, 22, 26, 29, 30, 31, 33, 34, 38]],                                                                       ['G351',2.0,100.228,227.991,[0]],                     ['G353',2.0,100.547,229.431,[7, 11, 14, 15, 21, 22, 26, 29, 35, 45, 46, 47]],                                                                                                          ['G008',3.4,100.526,228.732,[10,15]],          ['G010',4.95,100.704,229.268,[8,16,39,50]],                                                               ['G012',2.4,100.680,229.080,[11, 16, 27, 30, 33, 35, 36, 43, 46, 62]],                                                                                                             ['G333',4.2,100.756,229.062,[4, 7, 11, 26, 29, 49, 50, 62, 64, 67, 70, 90, 102, 105, 112, 115]],                                                                                                                                                                                                                     ['W51-IRS2',5.4,101.263,228.530,[7, 22, 25, 34, 41, 54, 68, 75, 97, 98, 108, 113, 118, 122, 123, 124]]                                                                               ]
#      																																																																																									  																													                                                       																																													      																																																																     																																																														   																		
# methode ratio des flux log(F6/F3)/log(v6/v3)<2, avec une correction lineaire pour les flux en B3
sources =     [['G327', 2.5, 101.776, 229.507, [0], 1.08, 1.53], ['G328', 2.5, 101.5, 227.575, [0], 1.08, 1.53], ['G337', 2.7, 101.602, 227.503, [0], 1.0, 1.13], ['G338', 3.9, 100.882, 229.226, [0], 0.69, 1.41], ['W43-MM1', 5.5, 99.759, 229.68, [0, 19], 0.49, 1.54], ['W43-MM2', 5.5, 101.017, 227.597, [0], 0.49, 1.60], ['W43-MM3', 5.5, 100.911, 228.931, [0, 7], 0.49, 1.33], ['W51-E', 5.4, 101.426, 228.918, [0, 1, 3, 5, 6, 9, 10, 27], 0.5, 3.86], ['G351', 2.0, 100.228, 227.991, [0, 3], 1.35, 0.80], ['G353', 2.0, 100.547, 229.431, [0, 1], 1.35, 1.18], ['G008', 3.4, 100.526, 228.732, [0, 1], 0.79, 1.53], ['G010', 4.95, 100.704, 229.268, [0, 1, 10, 12, 13, 15, 23, 28, 41], 0.55, 0.96], ['G012', 2.4, 100.68, 229.08, [0, 2, 7, 17, 21, 28, 37, 49, 50, 56], 1.13, 1.30], ['G333', 4.2, 100.756, 229.062, [0, 6, 9, 10, 12, 13, 16, 17, 20, 24, 27, 36, 40, 42, 48, 51, 53, 55, 56, 57, 58, 59, 60, 69, 82, 83, 86, 93, 95, 96, 97, 106, 107, 114, 118], 0.64, 1.28], ['W51-IRS2', 5.4, 101.263, 228.53, [0, 3, 5, 18, 36, 56, 58], 0.5, 1.64]]

## en confrontant les deux (je concatene, puis j'exclu les sources mesurees a 3mm qui ne montrent pas de freefree, et celle non mesurees qui auraient eu une mesure a 3mm avec alpha = 2)
sortie =[]
nboeil=0

#for source in sources:
#	cat = source[0]+'/catalog/'+source[0]+'-getsf.cat'
#	gets = np.loadtxt(cat, skiprows=122)
#	newlist  = source[4]
#	compt = 0
#	for sourceoeil in sourcesoeil[sources.index(source)][4]:
#		if sourceoeil not in newlist:
#			if gets[sourceoeil-1][46] < 1. and gets[sourceoeil-1][31]*(source[2]/source[3])**1 < gets[sourceoeil-1][49]*2:
#				newlist.append(sourceoeil)
#				compt+=1
#	newlist.sort()
#	field =[source[0],source[1],source[2],source[3],newlist,source[5]]
#	sortie.append(field)
#	nboeil=nboeil+compt
##	print(source[0]+' '+str(compt))
#print(nboeil)

for source in sources:
	cat = source[0]+'/catalog/'+source[0]+'-getsf.cat'
	gets = np.loadtxt(cat, skiprows=122)
	newlist  = source[4]
	compt = 0
	for sourceoeil in sourcesoeil[sources.index(source)][4]:
		if sourceoeil not in newlist:
			if math.log(gets[sourceoeil-1][33]/abs(gets[sourceoeil-1][50]))/math.log(source[3]/source[2])<1.:
				newlist.append(sourceoeil)
				compt+=1
	newlist.sort()
	field =[source[0],source[1],source[2],source[3],newlist,source[5],source[6]]
	sortie.append(field)
	nboeil=nboeil+compt
#	print(source[0]+' '+str(compt))
#print(nboeil)

sources = sortie

## liste des sources "propres" dans les regions intermediaires et evoluees
cleansources =  [['G010', 4.95, 100.704, 229.268, [7, 11, 17, 29, 33, 37, 40, 48, 58, 62]], ['G012', 2.4, 100.68, 229.08, [1, 4, 5, 8, 12, 13, 23, 26, 32, 34, 40, 47, 51, 53, 61, 66]], ['G333', 4.2, 100.756, 229.062,[3, 5, 8, 14, 18, 19, 22, 25, 35, 43, 44, 66, 77, 81, 84, 85, 108, 116]], ['W51-IRS2', 5.4, 101.263, 228.53, [4, 8, 9, 10, 12, 13, 14, 15, 16, 20, 24, 27, 28, 29, 31, 32, 35, 37, 43, 46, 47, 49, 51, 55, 57, 59, 60, 61, 62, 64, 71, 74, 76, 78, 79, 80, 81, 85, 86, 87, 88, 91, 95, 100, 104, 107, 109, 111, 114, 117, 119]]]


interevo = []
i=0
while i<len(cleansources):
	interevo.append(cleansources[i][0])
	i+=1
	pass

def Planck(temp, frequ):
	""" Retourne la fonction de Planck"""
	h = 6.62607004e-34 
	c = 299792458.
	k = 1.3806488E-23
	return 2*h*frequ**3/(c**2)*1./(math.exp(h*frequ/(k*temp))-1)

def computemass(flux, peack, omega, distance, beta = 1.5, temp = 20., frequ = 232.45):
	""" Permet de convertir un flux en Jy en mass (kg), distance en pc, frequ en Hz """
	frequ = frequ*1e9
	msun = 1.9891E30 #kg
	pc = 3.085677E16 #m
	kappa = 0.1*(frequ/1.0e12)**beta*0.1 # x0.1 pour passer de cm^2/g a des m^2/kg
	omega = ((omega/3600.*np.pi/180.)/2.)**2*np.pi
	mass = -omega*(distance*pc)**2/kappa*flux/peack*np.log(1-peack*1.e-26/(omega*Planck(temp,frequ)))
#	mass = flux*1.e-26*(distance*pc)**2/kappa*1./Planck(temp,frequ)
	return mass #en kg

def computemass1(flux, peack, omega, distance, beta = 1.5, temp = 20., frequ = 232.45):
	""" Permet de convertir un flux en Jy en mass (kg), distance en pc, frequ en Hz """
	frequ = frequ*1e9
	msun = 1.9891E30 #kg
	pc = 3.085677E16 #m
	kappa = 0.1*(frequ/1.0e12)**beta*0.1 # x0.1 pour passer de cm^2/g a des m^2/kg
	omega = ((omega/3600.*np.pi/180.)/2.)**2*np.pi
#	print(peack*1.e-26/(omega*Planck(temp,frequ)))
#	input()
	mass = -omega*(distance*pc)**2/kappa*flux/peack*np.log(1-peack*1.e-26/(omega*Planck(temp,frequ)))
#	masso = flux*1.e-26*(distance*pc)**2/kappa*1./Planck(temp,frequ)
#	print(mass/msun)
	return mass #en kg

def radius(major, minor, distance):
	""" Retourne le rayon moyen d une sphere de meme volume que l'ellipsoide rentree """
	au = 1.4959784260968E+11 #m
	a = major*distance*au/2
	b = minor*distance*au/2
	r = (a*b**2)**(1/3.)
	return r # en m

def alphavir(R,M,T=20):
	""" Retourne le parametre du viriel alpha_vir """
	msun = 1.9891E30 # kg
	G  = 6.6738480E-11 #m^3 kg^{-1} s^{-2}
	kB = 1.3806488E-23 # J/K
	mu = 2.4 # mean molecular weight per free particle  
	mp = 1.672627E-27  # kg
	gamma = 7/5. #index adiabatique
	sigma = (gamma*kB*T/(mu*mp))**0.5
	alpha = sigma**2*R/(3/5.*G*M)
	return alpha

def bonnorebert(R,T=20.):
	"""Retourne la masse de Bonnor Ebert"""
	msun = 1.9891E30 # kg
	G  = 6.6738480E-11 #m^3 kg^{-1} s^{-2}
	mp = 1.672627E-27  # kg
	mu = 2.8  
	kB = 1.3806488E-23 # J/K
	gamma = 1. # #index adiabatique
	sigma = (gamma*kB*T/(mu*mp))**0.5
	MBEc = 2.4*R*sigma**2/G	
	return MBEc#/msun

def compfreefree(mu1,mu3,unmm,troismm,alphadust,alphafree):
	old1mm = 0.
	i=0
	while (abs(unmm-old1mm))/unmm*100>0.1 and (i < 100):
		old1mm = unmm
		theo3mm = unmm*(mu3/mu1)**alphadust
		free3mm = troismm - theo3mm
		troismm = troismm - free3mm
		free1mm = free3mm*(mu1/mu3)**alphafree  
		unmm    = unmm-free1mm
		i+=1
	return unmm

def expomodinv(x,C=1.,mu=0.2,sigma=0.2,lammbda=0.2):
	return C*lammbda/2.*np.exp(lammbda/2.*(2*x-2*mu-lammbda*sigma**2))*(special.erf((mu+lammbda*sigma**2-x)/(2**0.5*sigma))+1.)


def estimateflux(field,sourcename,flux1mm,flux3mm,freq1mm,freq3mm,nbtries = 500000):
	alphadusts    = np.random.normal(3,1,nbtries)
	alphafreefree = np.random.normal(1,1,nbtries)
	thermique = []
	i=0 
	while i < len(alphadusts):
		thermique.append(compfreefree(freq1mm,freq3mm,flux1mm,flux3mm,alphadusts[i],alphafreefree[i]))
		i+=1
	thermique = np.array(thermique)
	thermique = thermique[np.where(thermique>0.)]
	thermique = thermique[np.where(thermique<2*flux1mm)]
	thermique.sort()
	i,k = 0,0
	nbcoup, valflu = [],[]
	step = (thermique[-1]-thermique[0])/100.
	palier = step
	while i < len(thermique):
		if thermique[i]>palier:
			nbcoup.append(k)
			valflu.append(palier-step/2.)
			palier=palier+step
			k=0
		i+=1
		k+=1
	valflu = np.array(valflu)
	nbcoup = np.array(nbcoup)
#	popt,pcov = curve_fit(expomodinv,valflu,nbcoup,p0=[9.5e+04,5.7e+01,4.9,3.6e-02],maxfev=100000)
	popt,pcov = curve_fit(expomodinv,valflu,nbcoup,maxfev=100000)
	result = expomodinv(valflu,popt[0],popt[1],popt[2],popt[3])
	mostlikely=valflu[np.where(result==max(result))[0][0]]
	i=0
	while result[i] < max(result)/2.:
		i+=1
	errorg = mostlikely-valflu[i]
	while result[i] > max(result)/2. and i<len(result)-1:
		i+=1
	errord = valflu[i]-mostlikely
	plt.figure(figsize=(20,15))
	plt.hist(thermique,bins=100)
	plt.plot(valflu,expomodinv(valflu,popt[0],popt[1],popt[2],popt[3]),color='grey',ls='--',lw=4.)
	plt.axvline(mostlikely,color='black',lw=2.)
	plt.axvline(mostlikely-errorg,color='green',lw=2.)
	plt.axvline(mostlikely+errord,color='green',lw=2.)
	plt.title('Flux estimate of source '+str(int(sourcename))+' in field '+str(field))
	plt.xlabel('Corrected flux at 1mm [Jy] (measured at '+str(round(flux1mm*1000.,4))+' mJy)')
	plt.ylabel('Number of counts')
	plt.savefig('freefree-estimate-'+str(field)+'-'+str(int(sourcename)),format='png')
	plt.show()
	return mostlikely,errorg,errord

#sources =     [['G327',2.5,101.776,229.507,[0]],['G328',2.5,101.500,227.575,[0]],['G337',2.7,101.602,227.503,[0]],        ['G338',3.9,100.882,229.226,[0]],     ['W43-MM1',5.5,99.759,229.680,[0]],        ['W43-MM2',5.5,101.017,227.597,[0]],                       ['W43-MM3',5.5,100.911,228.931,[0]],        ['W51-E',5.4,101.426,228.918,[0]],                                                                       ['G351',2.0,100.228,227.991,[0]],                     ['G353',2.0,100.547,229.431,[0]],                                                                                                          ['G008',3.4,100.526,228.732,[0]],          ['G010',4.95,100.704,229.268,[0]],                                                               ['G012',2.4,100.680,229.080,[0]],                                                                                                             ['G333',4.2,100.756,229.062,[0]],                                                                                                                                                                                                                     ['W51-IRS2',5.4,101.263,228.530,[0]]                                                                               ]
#print(computemass(2.562E-02,7e-3,0.1,2500.,frequ=229.507))
#print(computemass(1.245E-03,7e-3,0.1,2500.,frequ=229.507))

newfreefree = []
opacitycor = []

badmeas3mm=0

maxalphaBE=[]

os.system('rm quiet.txt')
os.system('rm active.txt')
fichier4 = open('quiet.txt', "w")
fichier5 = open('active.txt', "w")
fichier4.write( 'Number' + '          RA(User)' + '          DEC(User)' +'         GOOD' +'        SIG-MONO_B6'+'       Flux-tot'+ '       Error'+'           Grd-axe'+'         Pt-axe'+'        Angle'+'      Mass'+'      Err-Mass'+'      Mass BE'+'     alpha_BE'+'     SIG-MONO_B3'+'       '+'Flux_c 3mm '+ '     Err_c 3mm'+'     ratio B3/B6'+'     size[au] '+'      Temperature    \n')
fichier4.write( '  0   ' + '             1    ' + '              2    ' +'           3 ' +'             4     '+'           5   '+ '         6  '+'              7   '+'            8  '+'          9  '+'       10 '+'         11   '+'         12  '+'        13   '+'          14    '+'       '+'    15     '+ '        16    '+'          17    '+'        18    '+'         [19]        \n')
fichier5.write( 'Number' + '          RA(User)' + '          DEC(User)' +'         GOOD' +'        SIG-MONO_B6'+'       Flux-tot'+ '       Error'+'           Grd-axe'+'         Pt-axe'+'        Angle'+'      Mass'+'      Err-Mass'+'      Mass BE'+'     alpha_BE'+'     SIG-MONO_B3'+'       '+'Flux_c 3mm '+ '     Err_c 3mm'+'     ratio B3/B6'+'     size[au] '+'      Temperature    \n')
fichier5.write( '  0   ' + '             1    ' + '              2    ' +'           3 ' +'             4     '+'           5   '+ '         6  '+'              7   '+'            8  '+'          9  '+'       10 '+'         11   '+'         12  '+'        13   '+'          14    '+'       '+'    15     '+ '        16    '+'          17    '+'        18    '+'         [19]        \n')

os.system('rm doublemass.txt')
fichier7 = open('doublemass.txt', "w")

temppmap = []

toutout = 0

for source in sources:
	allsources=0
	freefreecore = []
	nofre        = []
	notbound     = []

	freefree = [0]
	#nom du catalogue de sortie et du fichier mapping de sortie
	fileout = source[0]
	######################################
	# data from .dat file
	cat = source[0]+'/catalog/'+source[0]+'-getsf.cat'
	gets = np.loadtxt(cat, skiprows=122)
	if tempcorrection:
		temperatures = 'PPMAP/'+source[0]+'_core_temperature_smooth_catalog.dat'
		temps = np.loadtxt(temperatures, skiprows=4)
		labelcore = []
		tempcore = []
		errtempcore =[]
		for val in temps:
			labelcore.append(int(val[0]))
			tempcore.append(val[3])
			errtempcore.append(int(round(temperror*val[3],0)))
#			tempcore.append(20)
#			errtempcore.append(5)
		temperatures = 'CH3OCHO_mom0_maps/'+source[0]+'.txt'
		tempsb = np.loadtxt(temperatures, skiprows=0)
		for val in tempsb:
			if val[0]!=0.0:
#				print(tempcore[labelcore.index(val[0])],val[1])
				tempcore[labelcore.index(val[0])]=val[1]
				errtempcore[labelcore.index(val[0])]=val[2]
			if val[0]>98.:
				nbhot +=1
		temps=[labelcore,tempcore,errtempcore]
		for val in tempcore:
			temppmap.append(val)
	# Ecriture des catalogues
	os.system('rm '+fileout+'-all.txt')
	os.system('rm '+fileout+'-filter.txt')
	os.system('rm '+fileout+'-nofreefree.txt')
	os.system('rm '+fileout+'-BE.txt')

	fichier  = open(fileout+'-all.txt', "w")
	fichier1 = open(fileout+'-filter.txt', "w")
	fichier2 = open(fileout+'-nofreefree.txt', "w")
	fichier3 = open(fileout+'-BE.txt', "w")
#	fichier6 = open(fileout+'-freefreecor.txt', "w")

	fichier.write( 'Number' + '          RA(User)' + '          DEC(User)' +'         GOOD' +'        SIG-MONO_B6'+'       Flux-tot'+ '       Error'+'           Grd-axe'+'         Pt-axe'+'        Angle'+'      Mass'+'      Err-Mass'+'      Mass BE'+'     alpha_BE'+'     SIG-MONO_B3'+'       '+'Flux_c 3mm '+ '     Err_c 3mm'+'     ratio B3/B6'+'     size[au] '+'      Temperature    \n')
	fichier.write( '  0   ' + '             1    ' + '              2    ' +'           3 ' +'             4     '+'           5   '+ '         6  '+'              7   '+'            8  '+'          9  '+'       10 '+'         11   '+'         12  '+'        13   '+'          14    '+'       '+'    15     '+ '        16    '+'          17    '+'        18    '+'          19        \n')
	fichier1.write('Number' + '          RA(User)' + '          DEC(User)' +'         GOOD' +'        SIG-MONO_B6'+'       Flux-tot'+ '       Error'+'           Grd-axe'+'         Pt-axe'+'        Angle'+'      Mass'+'      Err-Mass'+'      Mass BE'+'     alpha_BE'+'     SIG-MONO_B3'+'       '+'Flux_c 3mm '+ '     Err_c 3mm'+'     ratio B3/B6'+'     size[au] '+'      Temperature    \n')
	fichier1.write( '  0   ' + '             1    ' + '              2    ' +'           3 ' +'             4     '+'           5   '+ '         6  '+'              7   '+'            8  '+'          9  '+'       10 '+'         11   '+'         12  '+'        13   '+'          14    '+'       '+'    15     '+ '        16    '+'          17    '+'        18    '+'         19        \n')
	fichier2.write('Number' + '          RA(User)' + '          DEC(User)' +'         GOOD' +'        SIG-MONO_B6'+'       Flux-tot'+ '       Error'+'           Grd-axe'+'         Pt-axe'+'        Angle'+'      Mass'+'      Err-Mass'+'      Mass BE'+'     alpha_BE'+'     SIG-MONO_B3'+'       '+'Flux_c 3mm '+ '     Err_c 3mm'+'     ratio B3/B6'+'     size[au] '+'      Temperature  '+'    Peak_1mm '+'    Err_peak-1mm \n')
	fichier2.write( '  0   ' + '             1    ' + '              2    ' +'           3 ' +'             4     '+'           5   '+ '         6  '+'              7   '+'            8  '+'          9  '+'       10 '+'         11   '+'         12  '+'        13   '+'          14    '+'       '+'    15     '+ '        16    '+'          17    '+'        18    '+'          19       '+'       20    '+'          21     \n')
	fichier3.write('Number' + '          RA(User)' + '          DEC(User)' +'         GOOD' +'        SIG-MONO_B6'+'       Flux-tot'+ '       Error'+'           Grd-axe'+'         Pt-axe'+'        Angle'+'      Mass'+'      Err-Mass'+'      Mass BE'+'     alpha_BE'+'     SIG-MONO_B3'+'       '+'Flux_c 3mm '+ '     Err_c 3mm'+'     ratio B3/B6'+'     size[au] '+'      Temperature    \n')
	fichier3.write( '  0   ' + '             1    ' + '              2    ' +'           3 ' +'             4     '+'           5   '+ '         6  '+'              7   '+'            8  '+'          9  '+'       10 '+'         11   '+'         12  '+'        13   '+'          14    '+'       '+'    15     '+ '        16    '+'          17    '+'        18    '+'         19        \n')
#	fichier6.write('Number' + '          RA(User)' + '          DEC(User)' +'         GOOD' +'        SIG-MONO_B6'+'       Flux-tot'+ '       Error'+'           Grd-axe'+'         Pt-axe'+'        Angle'+'      Mass'+'      Err-Mass'+'      Mass BE'+'     alpha_BE'+'     SIG-MONO_B3'+'       '+'Flux_c 3mm '+ '     Err_c 3mm'+'     ratio B3/B6'+'     size[au] '+'      Temperature    \n')
#	fichier6.write( '  0   ' + '             1    ' + '              2    ' +'           3 ' +'             4     '+'           5   '+ '         6  '+'              7   '+'            8  '+'          9  '+'       10 '+'         11   '+'         12  '+'        13   '+'          14    '+'       '+'    15     '+ '        16    '+'          17    '+'        18    '+'         [19]        \n')
	i=0
	while i<len(gets):
		temperature = 20.0
		errtemp     = 5.0
		if tempcorrection:
			if gets[i][0] in temps[0]:
				temperature = temps[1][temps[0].index(gets[i][0])] 
				errtemp     = temps[2][temps[0].index(gets[i][0])]
		cor = gets[i][37]*gets[i][38]/(gets[i][54]*gets[i][55]) #correction de taille 1mm versus 3mm
		#ecriture des tables pour le papier
		if printtable:
			if round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(abs(gets[i][33]), abs(gets[i][31]), source[5],source[1]*1000.,frequ=source[3],temp=temperature),2) == 0.00:
				alphaBEprint = '<0.01'
			else:
				alphaBEprint = imposelength(imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(abs(gets[i][33]), abs(gets[i][31]), source[5],source[1]*1000.,frequ=source[3],temp=temperature),2)),4,True,'0'),5,False,' ')
			if round(computemass(abs(gets[i][34]),abs(gets[i][32]),source[5],source[1]*1000.,frequ=source[3],temp=temperature)/msun,1) == 0.0:
				errormass = '0.1'
			else:
				errormass = round(computemass(abs(gets[i][34]),abs(gets[i][32]),source[5],source[1]*1000.,frequ=source[3],temp=temperature)/msun,1)
			if round(gets[i][49]*1000,1) == 0.0:
				errorpeak3mm = '0.1'
			else:
				errorpeak3mm = round(gets[i][49]*1000,1)
			if gets[i][46]<0.1 or gets[i][48]<0.:
				badmeas3mm +=1
			# filtres de sasha + source freefree
			if (gets[i][9]>1. and gets[i][10]>1. and gets[i][31]/gets[i][32]>2. and gets[i][33]/gets[i][34]>2. and gets[i][37]/gets[i][38]<2. and gets[i][43]/gets[i][37]>1.15 and (gets[i][37]*gets[i][38])**0.5<maxsource*source[5]) and (int(gets[i][0]) in source[4]):
				# source detected at 3mm
				if gets[i][46]>0.1 and gets[i][48]>0. and gets[i][50]>0. and (1./cormax <= cor <= cormax):
					core = [int(gets[i][0]), degtoRAstr(gets[i][6])+degtoDECstr(gets[i][7]), degtoRA(gets[i][6]),degtoDEC(gets[i][7]),round(gets[i][37],2),round(gets[i][38],2),int(round(gets[i][41],0)),round(gets[i][31]*1000,1),round(gets[i][32]*1000,1),  round(gets[i][33]*1000,1),round(gets[i][34]*1000,1)   ,round(gets[i][54],2),round(gets[i][55],2),int(round(gets[i][58],0)),round(gets[i][48]*1000,1),errorpeak3mm,  round(gets[i][50]*cor*1000,1),round(gets[i][51]*cor*1000,1)      ,round(computemass(abs(gets[i][33]), abs(gets[i][31]), source[5],source[1]*1000.,frequ=source[3],temp=temperature)/msun,1),errormass,2*int(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,0)),alphaBEprint,round(math.log(abs(gets[i][33])/abs((gets[i][50]*cor)))/math.log(source[3]/source[2]),1), temperature, errtemp]
				else:
					core = [int(gets[i][0]), degtoRAstr(gets[i][6])+degtoDECstr(gets[i][7]), degtoRA(gets[i][6]),degtoDEC(gets[i][7]),round(gets[i][37],2),round(gets[i][38],2),int(round(gets[i][41],0)),round(gets[i][31]*1000,1),round(gets[i][32]*1000,1),  round(gets[i][33]*1000,1),round(gets[i][34]*1000,1)   ,'- ','-','-','-','-','-','-',round(computemass(abs(gets[i][33]), abs(gets[i][31]), source[5],source[1]*1000.,frequ=source[3],temp=temperature)/msun,1),errormass,2*int(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,0)),alphaBEprint,round(math.log(abs(gets[i][33])/abs((gets[i][50]*cor)))/math.log(source[3]/source[2]),1), temperature, errtemp]
				freefreecore.append(core)
				maxalphaBE.append(alphaBEprint)
			# filtres de sasha + source non contaminee
			if (gets[i][9]>1. and gets[i][10]>1. and  gets[i][31]/gets[i][32]>2. and gets[i][33]/gets[i][34]>2. and gets[i][37]/gets[i][38]<2. and gets[i][43]/gets[i][37]>1.15 and (gets[i][37]*gets[i][38])**0.5<maxsource*source[5]) and (int(gets[i][0]) not in source[4]):
				# source virialise
				if (bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(abs(gets[i][33]), abs(gets[i][31]), source[5],source[1]*1000.,frequ=source[3],temp=temperature)<2):
					if gets[i][46]>0.1 and gets[i][48]>0. and gets[i][50]>0. and (1./cormax <= cor <= cormax):
						core = [int(gets[i][0]), degtoRAstr(gets[i][6])+degtoDECstr(gets[i][7]), degtoRA(gets[i][6]),degtoDEC(gets[i][7]),round(gets[i][37],2),round(gets[i][38],2),int(round(gets[i][41],0)),round(gets[i][31]*1000,1),round(gets[i][32]*1000,1),  round(gets[i][33]*1000,1),round(gets[i][34]*1000,1)   ,round(gets[i][54],2),round(gets[i][55],2),int(round(gets[i][58],0)),round(gets[i][48]*1000,1),errorpeak3mm,  round(gets[i][50]*cor*1000,1),round(gets[i][51]*cor*1000,1)      ,round(computemass(abs(gets[i][33]), abs(gets[i][31]), source[5],source[1]*1000.,frequ=source[3],temp=temperature)/msun,1),errormass,2*int(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,0)),alphaBEprint,round(math.log(abs(gets[i][33])/abs((gets[i][50]*cor)))/math.log(source[3]/source[2]),1), temperature, errtemp]
					else:
						core = [int(gets[i][0]), degtoRAstr(gets[i][6])+degtoDECstr(gets[i][7]), degtoRA(gets[i][6]),degtoDEC(gets[i][7]),round(gets[i][37],2),round(gets[i][38],2),int(round(gets[i][41],0)),round(gets[i][31]*1000,1),round(gets[i][32]*1000,1),  round(gets[i][33]*1000,1),round(gets[i][34]*1000,1)   ,'- ','-','-','-','-','-','-',round(computemass(abs(gets[i][33]), abs(gets[i][31]), source[5],source[1]*1000.,frequ=source[3],temp=temperature)/msun,1),errormass,2*int(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,0)),alphaBEprint,round(math.log(abs(gets[i][33])/abs((gets[i][50]*cor)))/math.log(source[3]/source[2]),1), temperature, errtemp]
					nofre.append(core)
				# source non virialise
				else:
					if gets[i][46]>0.1 and gets[i][48]>0. and gets[i][50]>0. and (1./cormax <= cor <= cormax):
						core = [int(gets[i][0]), degtoRAstr(gets[i][6])+degtoDECstr(gets[i][7]), degtoRA(gets[i][6]),degtoDEC(gets[i][7]),round(gets[i][37],2),round(gets[i][38],2),int(round(gets[i][41],0)),round(gets[i][31]*1000,1),round(gets[i][32]*1000,1),  round(gets[i][33]*1000,1),round(gets[i][34]*1000,1)   ,round(gets[i][54],2),round(gets[i][55],2),int(round(gets[i][58],0)),round(gets[i][48]*1000,1),errorpeak3mm,  round(gets[i][50]*cor*1000,1),round(gets[i][51]*cor*1000,1)      ,round(computemass(abs(gets[i][33]), abs(gets[i][31]), source[5],source[1]*1000.,frequ=source[3],temp=temperature)/msun,1),errormass,2*int(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,0)),alphaBEprint,round(math.log(abs(gets[i][33])/abs((gets[i][50]*cor)))/math.log(source[3]/source[2]),1), temperature, errtemp]
					else:
						core = [int(gets[i][0]), degtoRAstr(gets[i][6])+degtoDECstr(gets[i][7]), degtoRA(gets[i][6]),degtoDEC(gets[i][7]),round(gets[i][37],2),round(gets[i][38],2),int(round(gets[i][41],0)),round(gets[i][31]*1000,1),round(gets[i][32]*1000,1),  round(gets[i][33]*1000,1),round(gets[i][34]*1000,1)   ,'- ','-','-','-','-','-','-',round(computemass(abs(gets[i][33]), abs(gets[i][31]), source[5],source[1]*1000.,frequ=source[3],temp=temperature)/msun,1),errormass,2*int(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,0)),alphaBEprint,round(math.log(abs(gets[i][33])/abs((gets[i][50]*cor)))/math.log(source[3]/source[2]),1), temperature,  errtemp]
					notbound.append(core)
		fichier.write(imposelength(str(int(gets[i][0])),4)+'         '+imposelength(str(gets[i][6]),12)+'       '+imposelength(str(gets[i][7]),12)+'        '+ imposelength(str(gets[i][10]),8) +'      '+ imposelength(str(gets[i][29]),10) +'     '+ imposelength(str("%.2e"%(abs(gets[i][33]))),9,False) +'     '+ imposelength(str("%0.2e"%(gets[i][34])),9,False) +'        '+ imposelength(str("%0.2e"%(gets[i][37])),9,False) +'      '+ imposelength(str("%0.2e"%(gets[i][38])),9,False) +'    '+ imposelength(str("%0.2e"%(gets[i][41])),9,False) + '    '+imposelength(str(round(computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+ '       '+imposelength(str(round(computemass(gets[i][34], gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+'        '+  imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.), T= temperature)/msun,2)),5) +'       '+ imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature),4)),6)+'      '+imposelength(str(gets[i][46]),10)+'       '+imposelength(str("%.2e"%(abs(gets[i][50]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][51]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][50]*cor/gets[i][33]))),9,False)+'       '+imposelength(str(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,2)),8,True)+'          '+imposelength(str(temperature),5,False,' ')+ '\n')
		# Filtering selon les recommendations de Sasha + filtre sur la taille des sources
		if gets[i][9]>1. and gets[i][10]>1. and  gets[i][31]/gets[i][32]>2. and gets[i][33]/gets[i][34]>2. and gets[i][37]/gets[i][38]<2. and gets[i][43]/gets[i][37]>1.15 and (gets[i][37]*gets[i][38])**0.5<maxsource*source[5]:
			fichier1.write(imposelength(str(int(gets[i][0])),4)+'         '+imposelength(str(gets[i][6]),12)+'       '+imposelength(str(gets[i][7]),12)+'        '+ imposelength(str(gets[i][10]),8) +'      '+ imposelength(str(gets[i][29]),10) +'     '+ imposelength(str("%.2e"%(abs(gets[i][33]))),9,False) +'     '+ imposelength(str("%0.2e"%(gets[i][34])),9,False) +'        '+ imposelength(str("%0.2e"%(gets[i][37])),9,False) +'      '+ imposelength(str("%0.2e"%(gets[i][38])),9,False) +'    '+ imposelength(str("%0.2e"%(gets[i][41])),9,False) + '    '+imposelength(str(round(computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+ '       '+imposelength(str(round(computemass(gets[i][34], gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+'        '+  imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.), T= temperature)/msun,2)),5) +'       '+ imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature),4)),6)+'      '+imposelength(str(gets[i][46]),10)+'       '+imposelength(str("%.2e"%(abs(gets[i][50]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][51]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][50]*cor/gets[i][33]))),9,False)+'       '+imposelength(str(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,2)),8,True)+'          '+imposelength(str(temperature),5,False,' ')+ '\n')
			# Si les sources sont detectees a 3mm
			if gets[i][46]>0.1 and gets[i][48]>0. and gets[i][50]>0. and (1./cormax <= cor <= cormax) :
				# Calcul de l'indice spectral, flag s'il est inferieur a 2
				if math.log(gets[i][33]/(gets[i][50]*cor))/math.log(source[3]/source[2])<2. :
					freefree.append(int(gets[i][0]))
					# si la source n'est pas identifiee commme etant contaminee par le freefree
			if int(gets[i][0]) not in source[4]:
				opacitycorrected = computemass( gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun
				notcorrected     = computemass1(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun
				percentage       = (opacitycorrected-notcorrected)/notcorrected*100
				if percentage > 20.:
					print(source[0], gets[i][0], round(notcorrected,1), round(opacitycorrected,1),round(percentage), temperature)
				fichier7.write(imposelength(str(round(opacitycorrected,2)),6,False," ")+ '       '+imposelength(str(round(notcorrected,2)),6, False," ")+'          '+imposelength(str(round(abs(notcorrected-opacitycorrected)/notcorrected*100.,1)),5,False," ")+    '\n')
				fichier2.write(imposelength(str(int(gets[i][0])),4)+'         '+imposelength(str(gets[i][6]),12)+'       '+imposelength(str(gets[i][7]),12)+'        '+ imposelength(str(gets[i][10]),8) +'      '+ imposelength(str(gets[i][29]),10) +'     '+ imposelength(str("%.2e"%(abs(gets[i][33]))),9,False) +'     '+ imposelength(str("%0.2e"%(gets[i][34])),9,False) +'        '+ imposelength(str("%0.2e"%(gets[i][37])),9,False) +'      '+ imposelength(str("%0.2e"%(gets[i][38])),9,False) +'    '+ imposelength(str("%0.2e"%(gets[i][41])),9,False) + '    '+imposelength(str(round(computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+ '       '+imposelength(str(round(computemass(gets[i][34], gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+'        '+  imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.), T= temperature)/msun,2)),5) +'       '+ imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature),4)),6)+'      '+imposelength(str(gets[i][46]),10)+'       '+imposelength(str("%.2e"%(abs(gets[i][50]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][51]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][50]*cor/gets[i][33]))),9,False)+'       '+imposelength(str(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,2)),8,True)+'          '+imposelength(str(temperature),5,False,' ')+'        '+imposelength(str("%.2e"%(abs(gets[i][31]))),9,False)+'        '+imposelength(str("%.2e"%(gets[i][32])),10)+   '\n')
				if computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun > source[6]:
					allsources+=1
#					if source[0]=='W51-IRS2':
#						print(gets[i][0])
#				fichier6.write(imposelength(str(int(gets[i][0])),4)+'         '+imposelength(str(gets[i][6]),12)+'       '+imposelength(str(gets[i][7]),12)+'        '+ imposelength(str(gets[i][10]),8) +'      '+ imposelength(str(gets[i][29]),10) +'     '+ imposelength(str("%.2e"%(abs(gets[i][33]))),9,False) +'     '+ imposelength(str("%0.2e"%(gets[i][34])),9,False) +'        '+ imposelength(str("%0.2e"%(gets[i][37])),9,False) +'      '+ imposelength(str("%0.2e"%(gets[i][38])),9,False) +'    '+ imposelength(str("%0.2e"%(gets[i][41])),9,False) + '    '+imposelength(str(round(computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+ '       '+imposelength(str(round(computemass(gets[i][34], gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+'        '+  imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.), T= temperature)/msun,2)),5) +'       '+ imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature),4)),6)+'      '+imposelength(str(gets[i][46]),10)+'       '+imposelength(str("%.2e"%(abs(gets[i][50]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][51]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][50]*cor/gets[i][33]))),9,False)+'       '+imposelength(str(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,2)),8,True)+'          '+imposelength(str(temperature),5,False,' ')+'        '+imposelength(str("%.2e"%(abs(gets[i][31]))),9,False)+  '\n')
#			if int(gets[i][0]) in source[4]:
#				if gets[i][46]>0.1 and gets[i][48]>0. and gets[i][50]>0. and (1./cormax <= cor <= cormax):
					#on a pu mesurer son indice spectral, elle est freefree et peut etre corrigee
#					print(source[0],gets[i][0],gets[i][33],gets[i][50]*cor,source[3],source[2])	
#					fluxunmm,errg,errd = estimateflux(source[0],gets[i][0],gets[i][33],gets[i][50]*cor,source[3],source[2])		
#					fichier6.write(imposelength(str(int(gets[i][0])),4)+'         '+imposelength(str(gets[i][6]),12)+'       '+imposelength(str(gets[i][7]),12)+'        '+ imposelength(str(gets[i][10]),8) +'      '+ imposelength(str(gets[i][29]),10) +'     '+ imposelength(str("%.2e"%(abs(fluxunmm))),9,False) +'     '+ imposelength(str("%0.2e"%(errg)),9,False) +'        '+ imposelength(str("%0.2e"%(gets[i][37])),9,False) +'      '+ imposelength(str("%0.2e"%(gets[i][38])),9,False) +'    '+ imposelength(str("%0.2e"%(gets[i][41])),9,False) + '    '+imposelength(str(round(computemass(fluxunmm,gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+ '       '+imposelength(str(round(computemass(errg, gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+'        '+  imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.), T= temperature)/msun,2)),5) +'       '+ imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(fluxunmm,gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature),4)),6)+'      '+imposelength(str(gets[i][46]),10)+'       '+imposelength(str("%.2e"%(abs(gets[i][50]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][51]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][50]*cor/fluxunmm))),9,False)+'       '+imposelength(str(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,2)),8,True)+'          '+imposelength(str(temperature),5,False,' ')+'        '+imposelength(str("%.2e"%(abs(gets[i][31]))),9,False)+  '\n')
				if source[0] in interevo:
					if (int(gets[i][0]) in cleansources[interevo.index(source[0])][4]):
						fichier4.write(imposelength(str(int(gets[i][0])),4)+'         '+imposelength(str(gets[i][6]),12)+'       '+imposelength(str(gets[i][7]),12)+'        '+ imposelength(str(gets[i][10]),8) +'      '+ imposelength(str(gets[i][29]),10) +'     '+ imposelength(str("%.2e"%(abs(gets[i][33]))),9,False) +'     '+ imposelength(str("%0.2e"%(gets[i][34])),9,False) +'        '+ imposelength(str("%0.2e"%(gets[i][37])),9,False) +'      '+ imposelength(str("%0.2e"%(gets[i][38])),9,False) +'    '+ imposelength(str("%0.2e"%(gets[i][41])),9,False) + '    '+imposelength(str(round(computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+ '       '+imposelength(str(round(computemass(gets[i][34], gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+'        '+  imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.), T= temperature)/msun,2)),5) +'       '+ imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature),4)),6)+'      '+imposelength(str(gets[i][46]),10)+'       '+imposelength(str("%.2e"%(abs(gets[i][50]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][51]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][50]*cor/gets[i][33]))),9,False)+'       '+imposelength(str(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,2)),8,True)+'          '+imposelength(str(temperature),5,False,' ')+'        '+imposelength(str("%.2e"%(abs(gets[i][31]))),9,False)+  '\n')
					else:
						fichier5.write(imposelength(str(int(gets[i][0])),4)+'         '+imposelength(str(gets[i][6]),12)+'       '+imposelength(str(gets[i][7]),12)+'        '+ imposelength(str(gets[i][10]),8) +'      '+ imposelength(str(gets[i][29]),10) +'     '+ imposelength(str("%.2e"%(abs(gets[i][33]))),9,False) +'     '+ imposelength(str("%0.2e"%(gets[i][34])),9,False) +'        '+ imposelength(str("%0.2e"%(gets[i][37])),9,False) +'      '+ imposelength(str("%0.2e"%(gets[i][38])),9,False) +'    '+ imposelength(str("%0.2e"%(gets[i][41])),9,False) + '    '+imposelength(str(round(computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+ '       '+imposelength(str(round(computemass(gets[i][34], gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+'        '+  imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.), T= temperature)/msun,2)),5) +'       '+ imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature),4)),6)+'      '+imposelength(str(gets[i][46]),10)+'       '+imposelength(str("%.2e"%(abs(gets[i][50]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][51]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][50]*cor/gets[i][33]))),9,False)+'       '+imposelength(str(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,2)),8,True)+'          '+imposelength(str(temperature),5,False,' ')+'        '+imposelength(str("%.2e"%(abs(gets[i][31]))),9,False)+  '\n')
				# Filtre des sources non virialise
				if bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)<2 :
					fichier3.write(imposelength(str(int(gets[i][0])),4)+'         '+imposelength(str(gets[i][6]),12)+'       '+imposelength(str(gets[i][7]),12)+'        '+ imposelength(str(gets[i][10]),8) +'      '+ imposelength(str(gets[i][29]),10) +'     '+ imposelength(str("%.2e"%(abs(gets[i][33]))),9,False) +'     '+ imposelength(str("%0.2e"%(gets[i][34])),9,False) +'        '+ imposelength(str("%0.2e"%(gets[i][37])),9,False) +'      '+ imposelength(str("%0.2e"%(gets[i][38])),9,False) +'    '+ imposelength(str("%0.2e"%(gets[i][41])),9,False) + '    '+imposelength(str(round(computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+ '       '+imposelength(str(round(computemass(gets[i][34], gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature)/msun,2)),6)+'        '+  imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.), T= temperature)/msun,2)),5) +'       '+ imposelength(str(round(bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.),T=temperature)/computemass(gets[i][33],gets[i][31], source[5], source[1]*1000.,frequ=source[3], temp = temperature),4)),6)+'      '+imposelength(str(gets[i][46]),10)+'       '+imposelength(str("%.2e"%(abs(gets[i][50]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][51]*cor))),9,False)+'      '+imposelength(str("%.2e"%(abs(gets[i][50]*cor/gets[i][33]))),9,False)+'       '+imposelength(str(round(radius(gets[i][37],gets[i][38],source[1]*1000.)/au,2)),8,True)+'          '+imposelength(str(temperature),5,False,' ')+ '\n')

		i+=1  
	coch='tets'
	if printtable:
#		print(source[0],notbound)
		table  = open(source[0]+'-catpap.tex', "w")
		table.write('\\begin{landscape}\n')
		table.write('\\begin{table}\n')
		table.write('\\caption{Compact sources extracted by \\textsl{getsf} in '+source[0]+', after smoothing.}\n')
		table.write('\\label{t:core-extraction-smoothed-'+source[0]+'}\n')
		table.write('\\centerline{ \n')    
		table.write('\\addtolength{\\tabcolsep}{-2pt}\n')  
		table.write('\\begin{tabular}{cccc|cccc|cccc|ccccccc}\n') 
		table.write('\\hline\\hline \n')
		table.write('\\noalign{\smallskip}\n')
		table.write('  &           &             &              &      \multicolumn{4}{c}{1.3\,mm}                                                               &  \multicolumn{4}{c}{3\,mm}                            \\\\ \n') 
		table.write('n & Core name & RA          & DEC          & $a  \\times b $                                  & PA    &$S^{\\rm Peak}$                   &$S^{\\rm int}$                 & $a \\times b$               & PA    &$S^{\\rm Peak}$                 &$S^{\\rm int}$            & T    & Mass$^{\,(1)}$& $Size^{\,(1)}$ & $\\alpha_{\\rm BE}^{\,(1)}$ & $\\alpha$ & $\\textsl{GExt2D}^{\,(2)}$\\\\ \n') 
		table.write('& & [IRCS]                & [ICRS]                 & [$\\arcsec\\times\\arcsec$]            & [deg]              & [mJy.Beam$^{-1}$]           & [mJy]      & [$\\arcsec\\times\\arcsec$]            & [deg]              & [mJy.Beam$^{-1}$]           & [mJy]          & [K]            & [$M_{\odot}$] & [au]   &                   &                      &          &                          \\\\ \n')
		table.write('\\noalign{\smallskip}\n')
		table.write('\hline \n')
		table.write('\\noalign{\smallskip} \n')
		tick = chargtab(source[0]+'-tick.txt')
		nbcore = 0
		for core in nofre:
			if float(core[18])>2.05:
				nbcore+=1
#		print(source[0],nbcore)
		nbcore=0
		for core in nofre:
			if core[21] != '$<10^{-2}$' and float(core[21])>2.	:
				print('changement!')
				print(source[0])
				print(core)
			nbcore+=1
			if nbcore==42:
				table.write('\end{tabular}}\n')
				table.write('\end{table}\n')
				table.write('\end{landscape}\n\n\n')
				table.write('\\begin{landscape}\n')
				table.write('\\begin{table}\n')
				table.write('\caption{Continuation of Table~\\ref{t:core-extraction-smoothed-'+source[0]+'}}\n')
				table.write('\label{t:core-extraction-smoothed-w51irs2-suite2}\n')    
				table.write('\centerline{\n')
				table.write('\\addtolength{\\tabcolsep}{-2pt}\n')  
				table.write('\\begin{tabular}{cccc|cccc|cccc|ccccccc}\n')
				table.write('  &           &             &              &      \multicolumn{4}{c}{1.3\,mm}                                                               &  \multicolumn{4}{c}{3\,mm}                            \\\\ \n') 
				table.write('n & Core name & RA          & DEC          & $a  \\times b $                                  & PA    &$S^{\\rm Peak}$                   &$S^{\\rm int}$                 & $a \\times b$               & PA    &$S^{\\rm Peak}$                 &$S^{\\rm int}$            & T    & Mass$^{\,(1)}$& $Size^{\,(1)}$ & $\\alpha_{\\rm BE}^{\,(1)}$ & $\\alpha$ & $\\textsl{GExt2D}^{\,(2)}$\\\\ \n') 
				table.write('\\noalign{\smallskip}\n')
				table.write('\hline \n')
				table.write('\\noalign{\smallskip} \n')
				nbcore=0 
			for valeur in tick:
				if int(valeur[0])==core[0]:
					coch = valeur[1]
			if str(core[12])=='-':
				table.write(imposelength(str(core[0]),3,False)+' & '+str(core[1])+' & '+str(core[2])+' & '+str(core[3])+' & '+ imposelength(str(core[4]),4,True,'0')+'$\\times$'+imposelength(str(core[5]),4,True,'0')+' & '+imposelength(str(core[6]),3,False,'\,\,\,')+ ' & '+ imposelength(str(core[7]),6,False,'\,\,\,')+'$\\pm$'+str(core[8]) +' & '+ imposelength(str(core[9]),5,False,'\,\,\,')+'$\\pm$'+imposelength(str(core[10]),3,True,'0') + ' &        -         &      -      &        \,\,\,-      &       \,\,\,\,\,\,-        & '+imposelength(str(int(core[23])),3,False,'\,\,\,')+'$\\pm$'+imposelength(str(int(core[24])),2,False,'0')+'      &    ' + imposelength(str(core[18]),5,False,'\,\,\,')+'$\\pm$'+imposelength(str(core[19]),4,True,'\,\,\,') + ' & '+ str(core[20])+ ' & '+ imposelength(str(core[21]),4,True,'0')+ ' &  -  & ' +coch+' \\\\ \n')
			else:
				table.write(imposelength(str(core[0]),3,False)+' & '+str(core[1])+' & '+str(core[2])+' & '+str(core[3])+' & '+ imposelength(str(core[4]),4,True,'0')+'$\\times$'+imposelength(str(core[5]),4,True,'0')+' & '+imposelength(str(core[6]),3,False,'\,\,\,')+ ' & '+ imposelength(str(core[7]),6,False,'\,\,\,')+'$\\pm$'+str(core[8]) +' & '+ imposelength(str(core[9]),5,False,'\,\,\,')+'$\\pm$'+imposelength(str(core[10]),3,True,'0') + ' & '+   imposelength(str(core[11]),4,True,'0')+'$\\times$'+imposelength(str(core[12]),4,True,'0')+' & '+imposelength(str(core[13]),3,False,'\,\,\,')+ ' & '+ imposelength(str(core[14]),5,False,'\,\,\,')+'$\\pm$'+imposelength(str(core[15]),4,True,'\,\,\,') +' & '+ imposelength(str(core[16]),5,False,'\,\,\,')+'$\\pm$'+str(core[17])+' & '+imposelength(str(int(core[23])),3,False,'\,\,\,')+'$\\pm$'+imposelength(str(int(core[24])),2,False,'0') +'      &    ' + imposelength(str(core[18]),5,False,'\,\,\,')+'$\\pm$'+imposelength(str(core[19]),4,True,'\,\,\,') + ' & '+ str(core[20])+ ' & '+ imposelength(str(core[21]),4,True,'0')+ ' & '+ str(core[22])+ ' & ' +coch+' \\\\ \n')
		if freefreecore != []:
			table.write('\\noalign{\smallskip} \n')
			table.write('\hline \n')
			table.write('\\noalign{\smallskip} \n')   
			for core in freefreecore:
				nbcore+=1
				if nbcore==42:
					table.write('\end{tabular}}\n')
					table.write('\end{table}\n')
					table.write('\end{landscape}\n\n\n')
					table.write('\\begin{landscape}\n')
					table.write('\\begin{table}\n')
					table.write('\caption{Continuation of Table~\\ref{t:core-extraction-smoothed-'+source[0]+'}}\n')
					table.write('\label{t:core-extraction-smoothed-w51irs2-suite2}\n')    
					table.write('\centerline{\n')
					table.write('\\addtolength{\\tabcolsep}{-2pt}\n')  
					table.write('\\begin{tabular}{cccc|cccc|cccc|cccccc}\n')
					table.write('  &           &             &              &      \multicolumn{4}{c}{1.3\,mm}                                                               &  \multicolumn{4}{c}{3\,mm}                            \\\\ \n') 
					table.write('n & Core name & RA          & DEC          & $a  \\times b $                                  & PA    &$S^{\\rm Peak}$                   &$S^{\\rm int}$                 & $a \\times b$               & PA    &$S^{\\rm Peak}$                 &$S^{\\rm int}$            & T    & Mass$^{\,(1)}$& $Size^{\,(1)}$ & $\\alpha_{\\rm BE}^{\,(1)}$ & $\\alpha$ & $\\textsl{GExt2D}^{\,(2)}$\\\\ \n') 
					table.write('\\noalign{\smallskip}\n')
					table.write('\hline \n')
					table.write('\\noalign{\smallskip} \n')
					nbcore=0 
				for valeur in tick:
					if int(valeur[0])==core[0]:
						coch = valeur[1]
				if str(core[12])=='-':
					table.write(imposelength(str(core[0]),3,False)+' & '+str(core[1])+' & '+str(core[2])+' & '+str(core[3])+' & '+ imposelength(str(core[4]),4,True,'0')+'$\\times$'+imposelength(str(core[5]),4,True,'0')+' & '+imposelength(str(core[6]),3,False,'\,\,\,')+ ' & '+ imposelength(str(core[7]),6,False,'\,\,\,')+'$\\pm$'+str(core[8]) +' & '+ imposelength(str(core[9]),5,False,'\,\,\,')+'$\\pm$'+imposelength(str(core[10]),3,True,'0') + ' &        -         &      -      &        \,\,\,-      &       \,\,\,\,\,\,-        &      -      &    -  & '+ str(core[20])+ ' & -  &  -  & ' +coch+' \\\\ \n')
				else:
					table.write(imposelength(str(core[0]),3,False)+' & '+str(core[1])+' & '+str(core[2])+' & '+str(core[3])+' & '+ imposelength(str(core[4]),4,True,'0')+'$\\times$'+imposelength(str(core[5]),4,True,'0')+' & '+imposelength(str(core[6]),3,False,'\,\,\,')+ ' & '+ imposelength(str(core[7]),6,False,'\,\,\,')+'$\\pm$'+str(core[8]) +' & '+ imposelength(str(core[9]),5,False,'\,\,\,')+'$\\pm$'+imposelength(str(core[10]),3,True,'0') + ' & '+   imposelength(str(core[11]),4,True,'0')+'$\\times$'+imposelength(str(core[12]),4,True,'0')+' & '+imposelength(str(core[13]),3,False,'\,\,\,')+ ' & '+ imposelength(str(core[14]),5,False,'\,\,\,')+'$\\pm$'+imposelength(str(core[15]),4,True,'\,\,\,') +' & '+ imposelength(str(core[16]),5,False,'\,\,\,')+'$\\pm$'+str(core[17])+' &  -      &    -   & '+ str(core[20])+ ' & -  & '+ str(core[22])+ ' & ' +coch+' \\\\ \n')
		table.write('\\noalign{\smallskip}\n')
		table.write('\hline\n')
		table.write('\\noalign{\smallskip}\n')
		table.write('\end{tabular}}\n')
		table.write('\\tablefoot\n')
		table.write('{\n')
		table.write('The coordinates are given at the J2000 Epoch. \n')
		table.write('(1) Estimated from the measures at 1.3\,mm. \n')
		table.write('(2) Indicates if the source found by \\textsl{getsf} was also found by \\textsl{GExt2D} (True or False). \n')
		table.write('The $^*$ next to the uncertainty indicates that the true value is $<$0.05. \n')
		table.write('}\n')
		table.write('\end{table}\n')	
		table.write('\\end{landscape}')
		table.close()

		tableascii  = open(source[0]+'-catalogue-smoothed.txt', "w")
		tableascii.write('  n     Core name            RA             DEC         a x b     PA   S_Peak-1.3mm   S_Int-1.3mm     a x b    PA     S_Peak-3mm   S_Int-3mm         T        Mass           Size  alpha_BE    alpha   GExt2D \n') 
		tableascii.write('                           [J2000]        [J2000]       ["x"]   [deg]   [mJy/beam]       [mJy]        ["x"]   [deg]   [mJy/beam]     [mJy]          [K]      [Msun]          [au]   \n') 
		tick = chargtab(source[0]+'-tick.txt')
#		if source[0]=='G012':
#			print(source[0],nofre[10])
		for core in nofre:
			for valeur in tick:
				if int(valeur[0])==core[0]:
					coch = valeur[1]
			if str(core[12])=='-':
				tableascii.write(imposelength(str(core[0]),3,False)+'  '+str('J'+core[1])+'  '+str(core[2])+'  '+str(core[3])+'  '+ imposelength(str(core[4]),4,True,'0')+'x'+imposelength(str(core[5]),4,True,'0')+'  '+imposelength(str(core[6]),3,False,' ')+ '   '+ imposelength(str(core[7]),6,False,' ')+'±'+imposelength(str(core[8]),4) +'    '+ imposelength(str(core[9]),6,False,' ')+'±'+imposelength(str(core[10]),4,True,' ') + '       -        -         -           -          '+imposelength(str(int(round(core[23]))),3,False)+'±'+imposelength(str(int(round(core[24]))),3)+'  '+ imposelength(str(core[18]),5,False)+'±'+imposelength(str(core[19]),7) +'     '+ str(core[20])+'    '+ core[21]+ '        -       ' +coch+'  \n')
			else:
				tableascii.write(imposelength(str(core[0]),3,False)+'  '+str('J'+core[1])+'  '+str(core[2])+'  '+str(core[3])+'  '+ imposelength(str(core[4]),4,True,'0')+'x'+imposelength(str(core[5]),4,True,'0')+'  '+imposelength(str(core[6]),3,False,' ')+ '   '+ imposelength(str(core[7]),6,False,' ')+'±'+imposelength(str(core[8]),4) +'    '+ imposelength(str(core[9]),6,False,' ')+'±'+imposelength(str(core[10]),4,True,' ') + '   '+   imposelength(str(core[11]),4,True,'0')+'x'+imposelength(str(core[12]),4,True,'0')+'  '+imposelength(str(core[13]),3,False,' ')+ '  '+ imposelength(str(core[14]),7,False,' ')+'±'+imposelength(str(core[15]),4,True,' ') +' '+ imposelength(str(core[16]),6,False,' ')+'±'+imposelength(str(core[17]),4)+'      '+imposelength(str(int(round(core[23]))),3,False)+'±'+imposelength(str(int(round(core[24]))),3)+'  '+ imposelength(str(core[18]),5,False)+'±'+imposelength(str(core[19]),3) +'         '+str(core[20]) +'    '+ core[21]+'      '+ imposelength(str(core[22]),4,False,' ')+ '      ' +coch+'  \n')
		if freefreecore != []:
			tableascii.write('\n')
			for core in freefreecore:
				for valeur in tick:
					if int(valeur[0])==core[0]:
						coch = valeur[1]
				if str(core[12])=='-':
					tableascii.write(imposelength(str(core[0]),3,False)+'  '+str('J'+core[1])+'  '+str(core[2])+'  '+str(core[3])+'  '+ imposelength(str(core[4]),4,True,'0')+'x'+imposelength(str(core[5]),4,True,'0')+'  '+imposelength(str(core[6]),3,False,' ')+ '   '+ imposelength(str(core[7]),6,False,' ')+'±'+imposelength(str(core[8]),4) +'    '+ imposelength(str(core[9]),6,False,' ')+'±'+imposelength(str(core[10]),4,True,' ') + '       -        -         -            -            -          -            '+ str(core[20])+ '      -          -       ' +coch+' \n')
				else:
					tableascii.write(imposelength(str(core[0]),3,False)+'  '+str('J'+core[1])+'  '+str(core[2])+'  '+str(core[3])+'  '+ imposelength(str(core[4]),4,True,'0')+'x'+imposelength(str(core[5]),4,True,'0')+'  '+imposelength(str(core[6]),3,False,' ')+ '   '+ imposelength(str(core[7]),6,False,' ')+'±'+imposelength(str(core[8]),4) +'    '+ imposelength(str(core[9]),6,False,' ')+'±'+imposelength(str(core[10]),4,True,' ') + '   '+   imposelength(str(core[11]),4,True,'0')+'x'+imposelength(str(core[12]),4,True,'0')+'  '+imposelength(str(core[13]),3,False,' ')+ '  '+ imposelength(str(core[14]),7,False,' ')+'±'+imposelength(str(core[15]),4,True,' ') +'  '+ imposelength(str(core[16]),6,False,' ')+'±'+imposelength(str(core[17]),4)+'        -          -            '+ str(core[20]) +'      -        '+ imposelength(str(core[22]),4,False,' ')+ '      ' +coch+' \n')
		tableascii.close()

	newfreefree.append([source[0],source[1],source[2],source[3],freefree,source[5]])
	fichier.close()
	fichier1.close()
	fichier2.close()
	fichier3.close()
	
	# Ecriture des fichiers .map
	os.system('rm '+fileout+'-all.map')
	os.system('rm '+fileout+'-filter.map')
	os.system('rm '+fileout+'-nofreefree.map')
	os.system('rm '+fileout+'-BE.map')

	fichier  = open(fileout+'-all.map', "w")
	fichier1 = open(fileout+'-filter.map', "w")
	fichier2 = open(fileout+'-nofreefree.map', "w")
	fichier3 = open(fileout+'-BE.map', "w")
	
	fichier2.write('set char 0.15 \n')
#	fichier1.write('set char 0.15 \n')
#	fichier2.write('@pen 0 2 \n')
	i=0
	while i<len(gets):
		fichier.write('ellipse '+str(gets[i][37]*arcsec/2)+' '+str(gets[i][38]*arcsec/2)+' '+str(90-gets[i][41])+' /user '+degtoRAmap(gets[i][6])+' '+degtoDECmap(gets[i][7])+' abs \n')
		# filtre de Sasha + taille max des sources
		if gets[i][9]>1. and gets[i][10]>1. and  gets[i][31]/gets[i][32]>2. and gets[i][33]/gets[i][34]>2. and gets[i][37]/gets[i][38]<2. and gets[i][43]/gets[i][37]>1.15 and (gets[i][37]*gets[i][38])**0.5<maxsource*source[5]:
			fichier1.write('ellipse '+str(gets[i][37]*arcsec/2)+' '+str(gets[i][38]*arcsec/2)+' '+str(90-gets[i][41])+' /user '+degtoRAmap(gets[i][6])+' '+degtoDECmap(gets[i][7])+' abs \n')
#			fichier1.write('dr text '+degtoRA(gets[i][6]+0*arcsec*180./np.pi)+' '+degtoDEC(gets[i][7]+0*arcsec*180./np.pi)+' "'+str(int(gets[i][0]))+'" /user abs \n')
			if int(gets[i][0]) in sourcesoeil[sources.index(source)][4]:
				fichier2.write('@pen 5 3 \n')
				fichier2.write('ellipse '+str(gets[i][37]*arcsec/2)+' '+str(gets[i][38]*arcsec/2)+' '+str(90-gets[i][41])+' /user '+degtoRAmap(gets[i][6])+' '+degtoDECmap(gets[i][7])+' abs \n')
				fichier2.write('@pen 0 2 \n')
			if int(gets[i][0]) in source[4]:
				fichier2.write('@pen 7 1 \n')
				fichier2.write('ellipse '+str(gets[i][37]*arcsec/2)+' '+str(gets[i][38]*arcsec/2)+' '+str(90-gets[i][41])+' /user '+degtoRAmap(gets[i][6])+' '+degtoDECmap(gets[i][7])+' abs \n')
				fichier2.write('@pen 6 1 \n')
				fichier2.write('ellipse '+str(gets[i][37]*arcsec/2)+' '+str(gets[i][38]*arcsec/2)+' '+str(90-gets[i][41])+' /user '+degtoRA(gets[i][6])+' '+degtoDEC(gets[i][7])+' abs \n')
#				fichier2.write('ellipse '+str(gets[i][54]*arcsec/2)+' '+str(gets[i][55]*arcsec/2)+' '+str(90-gets[i][58])+' /user '+degtoRAmap(gets[i][6])+' '+degtoDECmap(gets[i][7])+' abs \n')
#				fichier2.write('@pen 0 2 \n')
#				fichier2.write('dr text '+degtoRA(gets[i][6]+0*arcsec*180./np.pi)+' '+degtoDEC(gets[i][7]+0*arcsec*180./np.pi)+' "'+str(int(gets[i][0]))+'" /user abs \n')
#				fichier2.write('@pen 7 1 \n')
#				fichier2.write('dr text '+degtoRA(gets[i][6]+0*arcsec*180./np.pi)+' '+degtoDEC(gets[i][7]+0*arcsec*180./np.pi)+' "'+str(int(gets[i][0]))+'" /user abs \n')
#				fichier2.write('@pen 0 2 \n')					
			else:
				fichier2.write('@pen 7 2 \n')
				fichier2.write('ellipse '+str(gets[i][37]*arcsec/2)+' '+str(gets[i][38]*arcsec/2)+' '+str(90-gets[i][41])+' /user '+degtoRAmap(gets[i][6])+' '+degtoDECmap(gets[i][7])+' abs \n')
#				fichier2.write('dr text '+degtoRA(gets[i][6]+0*arcsec*180./np.pi)+' '+degtoDEC(gets[i][7]+0*arcsec*180./np.pi)+' "'+str(int(gets[i][0]))+'" /user abs \n')
				fichier2.write('@pen 0 1 \n')
				fichier2.write('ellipse '+str(gets[i][37]*arcsec/2)+' '+str(gets[i][38]*arcsec/2)+' '+str(90-gets[i][41])+' /user '+degtoRAmap(gets[i][6])+' '+degtoDECmap(gets[i][7])+' abs \n')
#				fichier2.write('dr text '+degtoRA(gets[i][6]+0*arcsec*180./np.pi)+' '+degtoDEC(gets[i][7]+0*arcsec*180./np.pi)+' "'+str(int(gets[i][0]))+'" /user abs \n')
#				if bonnorebert(radius(gets[i][37],gets[i][38],source[1]*1000.), T=temperature)/computemass(gets[i][33], gets[i][31], source[5], source[1]*1000.,frequ=source[3],temp=temperature)<2.:
#					fichier3.write('ellipse '+str(gets[i][37]*arcsec/2)+' '+str(gets[i][38]*arcsec/2)+' '+str(90-gets[i][41])+' /user '+degtoRAmap(gets[i][6])+' '+degtoDECmap(gets[i][7])+' abs \n')
		i+=1  
#	print(str(source[0])+', all:'+str(allsources))
#	toutout = toutout+allsources
	fichier2.write('set char 0.6 \n')
#	fichier1.write('set char 0.6 \n')

	fichier.close()
	fichier1.close()
	fichier2.close()
	fichier3.close()
#	fichier6.close()
#	fichier5.write('@ visualize.map '+source[0]+' -all.map \n')
#	fichier5.write('sic wait 1 \n')

#print(maxalphaBE)
#print(newfreefree)
#print('mapping @ visualize.map G008 BE')
fichier4.close()
fichier5.close()
fichier7.close()

#print(toutout)
#print(badmeas3mm)