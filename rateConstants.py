import os
import sys, math, time
from array import array
import csv
import numpy
import math
import scipy.optimize as optimization

def read_TS_corrections(reaction_name, T, pathOutput):
    file_name = pathOutput + "thermo-data-TS/" + reaction_name + "/" + reaction_name + ".data"
    f = open(file_name, "r")
    lines = f.readlines()
    T = str(T) + ' '
    for line in lines:
        if line.startswith(T):
            Hcorr = float(line.strip().split()[1])
            Gcorr = float(line.strip().split()[2])
            S = float(line.strip().split()[3])
    
    return Hcorr, Gcorr, S
        
def read_minima_corrections(species_name, T, pathOutput):
    file_name = pathOutput + "thermo-data-minima/" + species_name + "/" + species_name + ".data"
    f = open(file_name, "r")
    lines = f.readlines()
    
    T = str(T) + ' '
    for line in lines:
        if line.startswith(T):
            Hcorr = float(line.strip().split()[1])
            Gcorr = float(line.strip().split()[2])
            S = float(line.strip().split()[3])
    
    return Hcorr, Gcorr, S

def cTST(T,reaction_name, species, conv_f, conv_b, SCF, pathOutput):
    h = 6.62606957e-34
    kB = 1.3806488e-23
    R = 8.3144598
    
    SCF_A = SCF[0]
    SCF_B = SCF[1]
    SCF_TS = SCF[2]
    SCF_C = SCF[3]
    SCF_D = SCF[4]
    SCF_E = SCF[5]
    SCF_F = SCF[6]
    
    A = species[0]
    B = species[1]
    C = species[2]
    D = species[3]
    E = species[4]
    F = species[5]
        
    Hcorr_TS, Gcorr_TS, S_TS = read_TS_corrections(reaction_name,T,pathOutput)
    Hcorr_A, Gcorr_A, S_A = read_minima_corrections(A,T,pathOutput)
    Hcorr_C, Gcorr_C, S_C = read_minima_corrections(C,T,pathOutput) 
    
    G_TS = SCF_TS + Gcorr_TS
    G_A = SCF_A + Gcorr_A
    G_C = SCF_C + Gcorr_C
    
    H_TS = SCF_TS + Hcorr_TS
    H_A = SCF_A + Hcorr_A
    H_C = SCF_C + Hcorr_C
    
    
    if conv_f == 1:
        H_B = 0
        G_B = 0
        S_B = 0
    else:
        Hcorr_B, Gcorr_B, S_B = read_minima_corrections(B,T,pathOutput)
        H_B = SCF_B + Hcorr_B
        G_B = SCF_B + Gcorr_B
        
    if conv_b == 1:
        H_D = 0
        G_D = 0
        S_D = 0
        
        H_E = 0
        G_E = 0
        S_E = 0
        
        H_F = 0
        G_F = 0
        S_F = 0
    elif conv_b == 1000:
        Hcorr_D, Gcorr_D, S_D = read_minima_corrections(D,T,pathOutput)
        H_D = SCF_D + Hcorr_D
        G_D = SCF_D + Gcorr_D
        
        H_E = 0
        G_E = 0
        S_E = 0
        
        H_F = 0
        G_F = 0
        S_F = 0
    elif conv_b == 1000000:
        Hcorr_D, Gcorr_D, S_D = read_minima_corrections(D,T,pathOutput)
        H_D = SCF_D + Hcorr_D
        G_D = SCF_D + Gcorr_D
        
        Hcorr_E, Gcorr_E, S_E = read_minima_corrections(E,T,pathOutput)
        H_E = SCF_E + Hcorr_E
        G_E = SCF_E + Gcorr_E
        
        H_F = 0
        G_F = 0
        S_F = 0
    else:
        Hcorr_D, Gcorr_D, S_D = read_minima_corrections(D,T,pathOutput)
        H_D = SCF_D + Hcorr_D
        G_D = SCF_D + Gcorr_D
        
        Hcorr_E, Gcorr_E, S_E = read_minima_corrections(E,T,pathOutput)
        H_E = SCF_E + Hcorr_E
        G_E = SCF_E + Gcorr_E
        
        Hcorr_F, Gcorr_F, S_F = read_minima_corrections(F,T,pathOutput)
        H_F = SCF_F + Hcorr_F
        G_F = SCF_F + Gcorr_F    
    
    dGf = (G_TS - G_A - G_B)*627.509*4184
    dGb = (G_TS - G_C - G_D - G_E - G_F)*627.509*4184
    
    if dGf < 0:
        dGf = 0
    
    if dGb < 0:
        dGb = 0  
    
    dHf = (H_TS - H_A - H_B)*627.509*4184
    dHb = (H_TS - H_C - H_D - H_E - H_F)*627.509*4184
    
    dSf = (S_TS - S_A - S_B)
    dSb = (S_TS - S_C - S_D - S_E - S_F)
    
    kf = (kB*T/h)*(conv_f)*math.exp(-dGf/(R*T))
    kb = (kB*T/h)*(conv_b)*math.exp(-dGb/(R*T))
    
    # modified for gas phase: C = P/RT (Standard state concentration with P = 1 atm reference)
    if conv_f == 1:
        n = 1
    elif conv_f == 1000:
        n = 2
    if conv_b == 1:
        m = 1
    elif conv_b == 1000:
        m = 2
    elif conv_b == 1000000:
        m = 3
    elif conv_b == 1000000000:
        m = 4
    
    kf = (kB*T/h)*((T*82.057)**(n-1))*math.exp(-dGf/(R*T))
    kb = (kB*T/h)*((T*82.057)**(m-1))*math.exp(-dGb/(R*T))
    
    return kf,kb, dGf/4184, dGb/4184, dHf/4184, dHb/4184, dSf, dSb        

def Wigner_tunneling(img_freq,T):
    h = 6.62606957e-34
    kB = 1.3806488e-23
    img_freq = abs(img_freq)* 29979245800.0
    tau_wig = 1.0 + (1.0/24.0)*((h*img_freq)/(kB*T))**2
    return tau_wig
    
def Eckart_tunneling(TEMPERATURE, SCF, ZPE, mu, fc):
    # PHYSICAL CONSTANTS
    GAS_CONSTANT = 8.3144621
    PLANCK_CONSTANT = 6.62606957e-34
    BOLTZMANN_CONSTANT = 1.3806488e-23
    SPEED_OF_LIGHT = 2.99792458e10
    AVOGADRO_CONSTANT = 6.0221415e23
    AMU_to_KG = 1.66053886E-27
    autokcal = 627.509541
    kjtokcal = 4.184
    atmos = 101.325
    PI = 3.14159265359
    k = 3.1668114E-6 #Boltzmann Constant in atomic units 

    #10-point Gauss-Legendre Quadrature abscissa and weight (exact solution for up to 21st order polynomial)
    x = array('d',[-0.9739065285,-0.8650633667,-0.6794095683,-0.4333953941,-0.1488743390,0.1488743390,0.4333953941,0.6794095683,0.8650633667,0.9739065285])
    w = array('d',[0.0666713443,0.1494513492,0.2190863625,0.2692667193,0.2955242247,0.2955242247,0.2692667193,0.2190863625,0.1494513492,0.0666713443])

    #Parameters B, ALPHA, a, b, d of Eckart Potential
    def Bee(V_max,V_r, V_p):
	    bee = (V_max ** 0.5 + ((V_max - (V_p - V_r))) ** 0.5) ** 2
	    return bee

    def ALPHA(B,F_s,V_max,V_r, V_p):
	    alpha = (B * F_s / (2 * V_max * (V_max - (V_p - V_r)))) ** 0.5
	    return alpha

    def A(E,mu,ALPHA):
	    a =  2 * PI * (2 * mu * E)**0.5 / ALPHA
	    return a

    def B(E,mu,V_p,V_r,ALPHA):
	    b =  2 * PI * (2 * mu * ((E - (V_p - V_r))))**0.5 / ALPHA
	    return b

    def D(bee,mu,ALPHA):
	    d =  2 * PI * abs((2 * mu * bee - (ALPHA/2)**2))**0.5 / ALPHA
	    return d

    #Calculation of Transmission Probabilty of Eckart Potential
    def T(a,b,d):
        if d > 700:
            return 0
        else:
            T = (math.cosh(a+b) - math.cosh(a-b))/(math.cosh(a+b) + math.cosh(d))
            return T
	
    #Calculation of SINH function of Kappa
    def S(V_max,E):
	    S = math.sinh(((V_max-E)) / (TEMPERATURE*k))
	    return S
    
    SCF_A = SCF[0]
    SCF_B = SCF[1]
    SCF_TS = SCF[2]
    SCF_C = SCF[3]
    SCF_D = SCF[4]
    SCF_E = SCF[5]
    SCF_F = SCF[6]
    
    ZPE_A = ZPE[0]
    ZPE_B = ZPE[1]
    ZPE_TS = ZPE[2]
    ZPE_C = ZPE[3]
    ZPE_D = ZPE[4]
    ZPE_E = ZPE[5]
    ZPE_F = ZPE[6]
    
    E_r = SCF_A + SCF_B
    E_p = SCF_C + SCF_D + SCF_E + SCF_F
    
    ZPE_r = ZPE_A + ZPE_B
    ZPE_p = ZPE_C + ZPE_D + ZPE_E + ZPE_F

    V_r = (E_r + ZPE_r)
    V_p = (E_p + ZPE_p)
    V_max = (SCF_TS + ZPE_TS)
    
    if V_r > V_p:
	    E_o = V_r
    else:
	    E_o = V_p

    #Scaling of Energies(define V_r == 0)
    V_max = V_max - V_r
    V_p = V_p - V_r
    E_o = E_o - V_r
    V_r = V_r - V_r
    
    if V_max < 0:
        return 1
    if V_max - (V_p - V_r) < 0:
        return 1
       
    y = (V_max - E_o)/2.0 
    z = (V_max + E_o)/2.0 
    
    # Specifing Parameters for the Eckart Potential
    mu = mu*1836
    F_s = fc/15.569141
    bee = Bee(V_max,V_r,V_p)
    alpha = ALPHA(bee,F_s,V_max,V_r,V_p)
    d = D(bee,mu,alpha)
    
    #Calculation of Eckart tunneling correction using 10-point Gauss-Legendre Quadrature
    kappa = 1
    for i in range(0,10):
	    a = A((x[i] * y + z),mu,alpha)
	    b = B((x[i] * y + z),mu,V_p,V_r,alpha)
	    kappa = (2 * y  / (TEMPERATURE*k) * w[i] * S((V_max),(x[i] * y + z)) * T(a,b,d)) + kappa
    
    return kappa

def calculate_rate_constants(pathTS,pathOutput):
    def read_TS_data(reaction_name,pathOutput):
        file_name = pathOutput+"TS-SCF.txt"
        f = open(file_name, "r")
        lines = f.readlines()
        for line in lines:
            if line.startswith(reaction_name):
                SCF = float(line.strip().split()[1])
                ZPE = float(line.strip().split()[2])
                mu = float(line.strip().split()[3])
                fc = float(line.strip().split()[4])
                img_freq = float(line.strip().split()[5])
        
        return SCF, ZPE, mu, fc, img_freq
        
    def read_minima_data(species_name, pathOutput):
        species_name = species_name + ' '
        file_name = pathOutput + "minima-SCF.txt"
        f = open(file_name, "r")
        lines = f.readlines()
        for line in lines:
            if line.startswith(species_name):
                SCF = float(line.strip().split()[1])
                ZPE = float(line.strip().split()[2])
        
        return SCF, ZPE  
       
    for file_name in os.listdir(pathTS):
        if file_name.endswith('.log'):
            reaction_name = file_name.split('.log')[0]

            lines_to_write = []
            temperatures = [300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
                      
            SCF_TS, ZPE_TS, mu, fc, img_freq = read_TS_data(reaction_name, pathOutput)
            
            reactants,products = reaction_name.split('-TS-')
        
            if len(reactants.split('-p-'))==1:
                A = reactants
                B = ''
                conv_f = 1
                SCF_A, ZPE_A = read_minima_data(A,pathOutput)
                SCF_B, ZPE_B = 0, 0
                
            if len(reactants.split('-p-'))==2:
                A,B = reactants.split('-p-')
                conv_f = 1000
                SCF_A, ZPE_A = read_minima_data(A,pathOutput)
                SCF_B, ZPE_B = read_minima_data(B,pathOutput)
                        
            if len(products.split('-p-'))==1:
                C = products
                D = ''
                E = ''
                F = ''
                conv_b = 1
                SCF_C, ZPE_C = read_minima_data(C,pathOutput)
                SCF_D, ZPE_D = 0, 0
                SCF_E, ZPE_E = 0, 0
                SCF_F, ZPE_F = 0, 0
                              
            if len(products.split('-p-'))==2:
                C,D = products.split('-p-')
                E = ''
                F = ''
                conv_b = 1000
                SCF_C, ZPE_C = read_minima_data(C,pathOutput)
                SCF_D, ZPE_D = read_minima_data(D,pathOutput)
                SCF_E, ZPE_E = 0, 0
                SCF_F, ZPE_F = 0, 0
                                
            if len(products.split('-p-'))==3:
                C,D,E = products.split('-p-')
                F = ''
                conv_b = 1000000
                SCF_C, ZPE_C = read_minima_data(C,pathOutput)
                SCF_D, ZPE_D = read_minima_data(D,pathOutput)
                SCF_E, ZPE_E = read_minima_data(E,pathOutput)
                SCF_F, ZPE_F = 0, 0
                
            if len(products.split('-p-'))==4:
                C,D,E,F = products.split('-p-')
                conv_b = 1000000000
                SCF_C, ZPE_C = read_minima_data(C,pathOutput)
                SCF_D, ZPE_D = read_minima_data(D,pathOutput)
                SCF_E, ZPE_E = read_minima_data(E,pathOutput)
                SCF_F, ZPE_F = read_minima_data(F,pathOutput)
            
            SCF = [SCF_A, SCF_B, SCF_TS, SCF_C, SCF_D, SCF_E, SCF_F]
            ZPE = [ZPE_A, ZPE_B, ZPE_TS, ZPE_C, ZPE_D, ZPE_E, ZPE_F]
            species = [A, B, C, D, E, F]
                
            for T in temperatures:
                kf_cTST, kb_cTST, dGf, dGb, dHf, dHb, dSf, dSb = cTST(T,reaction_name, species, conv_f, conv_b, SCF, pathOutput)
                
                tau_wig = Wigner_tunneling(img_freq,T)
                kf_cTST_Wig = kf_cTST * tau_wig
                kb_cTST_Wig = kb_cTST * tau_wig
                
                tau_Eck_f = Eckart_tunneling(T, SCF, ZPE, mu, fc)
                kf_cTST_Eck = kf_cTST * tau_Eck_f
                kb_cTST_Eck = kb_cTST * tau_Eck_f
                               
                contents = "%-20d %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.6E %15.6E %15.6E %15.6E %15.6E %15.6E" %(T, dGf, dGb, dHf, dHb, dSf, dSb, kf_cTST, kf_cTST_Wig, kf_cTST_Eck, kb_cTST, kb_cTST_Wig, kb_cTST_Eck)
                lines_to_write.append(contents)
            
            lines_to_write = sorted(list(set(lines_to_write)),key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item))     
            text = "\n".join(lines_to_write) + "\n"
            
            header1 = "%-20s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s" %('Temperature', 'dGf (kcal/mol)', 'dGb (kcal/mol)', 'dHf (kcal/mol)', 'dHb (kcal/mol)', 'dSf (cal/mol)', 'dSb (cal/mol)','kf(cTST)', 'kf(cTST/Wig)', 'kf(cTST/Eck)','kb(cTST)', 'kb(cTST/Wig)', 'kb(cTST/Eck)') + '\n'
            header2 = "%-20s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s" %('-----------', '--------------', '--------------', '--------------', '--------------', '---------------', '--------------', '---------------', '----------------', '------------','---------------', '----------------', '------------') + '\n' 
            
            rate_constants_dir = pathOutput+"rate-constants/"
                  
            if not os.path.exists(rate_constants_dir):
                os.makedirs(rate_constants_dir)
            
            with open(os.path.join(rate_constants_dir,reaction_name + '.txt'), "w") as File:
                File.writelines(header1+header2)
                File.writelines(text)

def get_symmetry_factors(pathCSV,fileCSV):
    # Get symmetry factors

    ifile  = open(os.path.join(pathCSV,fileCSV), "r")
    reader = csv.reader(ifile)

    rownum = 0
    species_list = []
    TS_list = []
    reactions_list = []
    symmetry = []
     
    for row in reader:
        
        if rownum < 3:
            pass
        else:
            reactants,products = row[0].split('=')
            species = reactants.split('+') + products.split('+')
            species_list.extend(species)
            TS = "-p-".join(reactants.split('+')) + '-TS-' + "-p-".join(products.split('+'))
            TS_list.append(TS)
            reactions_list.append(row[0])
            symmetry.append([row[0],row[1],row[2]])        
               
        rownum = rownum + 1
    
    return symmetry
    
def curvefit_arrhenius(pathOutput,pathCSV,fileCSV):
    R = 8.3144621
    Path = pathOutput + "rate-constants/"
    Path2final = pathOutput

    def kf(T, A, n, Ea):
                return A*(T**n)*numpy.exp(-Ea/R/T)
                
    def modified_arrhenius(T, lnA, n, Ea):
                return -Ea/R/T + n*numpy.log(T) + lnA
                
    def arrhenius(T, lnA, Ea):
                return -Ea/R/T + lnA

    lines_to_write = []
    for file_name in os.listdir(Path):
        if file_name.endswith('.txt'):
            reaction_name = file_name.split('.txt')[0]
            
            reactants, products = reaction_name.split('-TS-')
            reactants = '+'.join(reactants.split('-p-'))
            products = '+'.join(products.split('-p-'))
            
            reaction_name = '='.join([reactants,products])
            
            data = numpy.genfromtxt(Path + file_name,skip_header=2)
            xdata = data[3:18,0]
            ydata = data[3:18,9]
            ydata = numpy.log(data[3:18,9])
            
            dHf = numpy.mean(data[3:18,3])
            
            x0 = numpy.array([5.0,0.0,dHf])
            x0_form2 = numpy.array([5.0,dHf])
            sigma = numpy.ones((ydata.size,))
            
            form2 = False    
            #coeff = optimization.curve_fit(arrhenius, xdata, ydata, x0_form2, sigma)
            #coeff = optimization.curve_fit(arrhenius, xdata, ydata, x0_form2, sigma, bounds = ([-numpy.inf, -numpy.inf],[numpy.inf, numpy.inf]))
            
            form1 = True
            coeff = optimization.curve_fit(modified_arrhenius, xdata, ydata, x0, sigma)
            #coeff = optimization.curve_fit(modified_arrhenius, xdata, ydata, x0, sigma, bounds = ([-numpy.inf, -1.0, -numpy.inf],[numpy.inf, 1.0, numpy.inf]))
            
            a = coeff[0]
            perr = numpy.sqrt(numpy.diag(coeff[1]))
            
            if form1:             
                A = math.exp(a[0])
                n = a[1]
                Ea = a[2]/4.184
                
            if form2:             
                A = math.exp(a[0])
                n = 0.0
                Ea = a[1]/4.184
                
            for sym in get_symmetry_factors(pathCSV,fileCSV):
                if sym[0]== reaction_name:
                    for_sym = float(sym[1])
                    back_sym = float(sym[2])
                
            A = for_sym*A
            
            contents = "%-50s %15.3E %15.2f %15.3f" %(reaction_name, A, n, Ea)

            lines_to_write.append(contents)
            
            k = A*(xdata**n)*numpy.exp(-Ea*4.184/R/xdata)
            
            kf = for_sym*numpy.exp(ydata)
            
            SSE = numpy.sum(numpy.square(kf-k))
            avg = numpy.sum(kf)
            SST = numpy.sum(numpy.square(kf-avg))
            
            R_square = 1 - (SSE/SST)
            print R_square, reaction_name
            
            f = open(Path+file_name.split('.txt')[0]+'.data', 'w')
            numpy.savetxt(f,numpy.c_[xdata, kf, k, kf/k], fmt='%5d %15.3E %15.3E %10.2f', header='%4s %14s %15s %10s'%('T(K)', 'Data', 'Fit', 'Data/fit'), footer='\n')
            numpy.savetxt(f,numpy.c_[A, n, Ea], fmt='%15.3E %15.3f %15.3f',header='%10s %15s %22s'%('A', 'n', 'Ea (cal/mol)'))
            numpy.savetxt(f,numpy.c_[perr[0], perr[1], perr[2]], fmt='%15.3E %15.3f %15.3f')
            numpy.savetxt(f,numpy.c_[R_square], fmt='%15.4f', header='')
            f.close()
              
    lines_to_write = sorted(list(set(lines_to_write)),key=str.lower)     
    text = "\n".join(lines_to_write) + "\n"

    header1 = "%-50s %15s %15s %15s" %('Reaction', 'A', 'n', 'Ea(cal/mol)') + '\n'
    header2 = "%-50s %15s %15s %15s" %('--------------------------------------------------', '---------------', '---------------', '---------------') + '\n' 

    with open(os.path.join(Path2final,'kinetics.txt'), "w") as File:
        File.writelines(header1+header2)
        File.writelines(text)
        
###############################################################################
# Change input output paths and set output file names here
###############################################################################
base_directory = '/gpfs/group/umt/default/HMX/nitramine_gas_phase_mechanism/improvements-to-caltech-mechanism/final-gas-phase-mechanism/'

pathCSV = base_directory
pathB3LYP = base_directory + 'M062X-minima'
pathTS = base_directory + 'M062X-TS'
pathOutput = base_directory + 'Output-files/'

fileCSV = 'PSU-gas-phase-mechanism.csv'                  
###############################################################################
calculate_rate_constants(pathTS,pathOutput)

curvefit_arrhenius(pathOutput,pathCSV, fileCSV)


