import os
import sys, math, time
from array import array
import csv
import numpy
import math
import scipy.optimize as optimization

def get_barrierless_reactions(pathCSV, fileCSV, pathTS):
    ifile  = open(os.path.join(pathCSV,fileCSV), "r")
    reader = csv.reader(ifile)

    rownum = 0
    TS_list = []
    reactions_list = []
    symmetry = []
     
    for row in reader:
        
        if rownum < 3:
            pass
        else:
            reactants,products = row[0].split('=')
            species = reactants.split('+') + products.split('+')
            TS = "-p-".join(reactants.split('+')) + '-TS-' + "-p-".join(products.split('+'))
            TS_list.append(TS)
            reactions_list.append(row[0])
            symmetry.append([row[0],row[1],row[2]])        
               
        rownum = rownum + 1

    file_list = []
    file_list_TS = []    

    for file_name in os.listdir(pathTS):
        if file_name.endswith('.log'):
            file_list_TS.append(file_name.split('.log')[0])

    NO_TS_list = []
    print '\nFiles missing from the TS folder'
    for item in TS_list:
        if item not in file_list_TS:
            print item
            NO_TS_list.append(item)

    print '\n'
              
    ifile.close()
    return NO_TS_list
       
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
    SCF_C = SCF[2]
    SCF_D = SCF[3]
    SCF_E = SCF[4]
    SCF_F = SCF[5]
    
    A = species[0]
    B = species[1]
    C = species[2]
    D = species[3]
    E = species[4]
    F = species[5]
        
    Hcorr_A, Gcorr_A, S_A = read_minima_corrections(A,T,pathOutput)
    Hcorr_C, Gcorr_C, S_C = read_minima_corrections(C,T,pathOutput) 
    
    G_A = SCF_A + Gcorr_A
    G_C = SCF_C + Gcorr_C
    
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
    
    dHr = (H_C + H_D + H_E + H_F - H_A - H_B)*627.509*4184   
    dGr = (G_C + G_D + G_E + G_F - G_A - G_B)*627.509*4184
    dSr = (S_C + S_D + S_E + S_F - S_A - S_B)
    
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
    
    if dHr > 0:
        dGf = dHr
        dGb = 0
    else:
        dGf = 0
        dGb = -dHr
    
    kf = (kB*T/h)*((T*82.057)**(n-1))*math.exp(-dGf/(R*T))
    kb = (kB*T/h)*((T*82.057)**(m-1))*math.exp(-dGb/(R*T))
    
    return kf,kb, dGr/4184, dHr/4184, dSr        

def calculate_rate_constants(pathTS,pathOutput,listBarrierless):       
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
       
    for reaction_name in listBarrierless:
        lines_to_write = []
        temperatures = [300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
                         
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
        
        SCF = [SCF_A, SCF_B, SCF_C, SCF_D, SCF_E, SCF_F]
        ZPE = [ZPE_A, ZPE_B, ZPE_C, ZPE_D, ZPE_E, ZPE_F]
        species = [A, B, C, D, E, F]
            
        for T in temperatures:
            kf_cTST, kb_cTST, dGr, dHr, dSr = cTST(T,reaction_name, species, conv_f, conv_b, SCF, pathOutput)
                                       
            contents = "%-20d %15.3f %15.3f %15.3f %15.6E %15.6E" %(T, dGr, dHr, dSr, kf_cTST, kb_cTST)
            lines_to_write.append(contents)
    
        lines_to_write = sorted(list(set(lines_to_write)),key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item))     
        text = "\n".join(lines_to_write) + "\n"
        
        header1 = "%-20s %15s %15s %15s %15s %15s" %('Temperature', 'dGr (kcal/mol)', 'dHr (kcal/mol)', 'dSr (cal/mol)', 'kf(cTST)', 'kb(cTST)') + '\n'
        header2 = "%-20s %15s %15s %15s %15s %15s" %('-----------', '--------------', '--------------', '--------------', '--------------', '---------------') + '\n' 
        
        rate_constants_dir = pathOutput+"rate-constants-barrierless/"
              
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
    Path = pathOutput + "rate-constants-barrierless/"
    Path2final = pathOutput
    
    if not os.path.exists(Path):
        os.mkdir(Path)

    def kf(T, A, n, Ea):
                return A*(T**n)*numpy.exp(-Ea/R/T)
                
    def logkf(T, lnA, n, Ea):
                return -Ea/R/T + n*numpy.log(T) + lnA

    lines_to_write = []
    for file_name in os.listdir(Path):
        if file_name.endswith('.txt'):
            reaction_name = file_name.split('.txt')[0]
            
            reactants, products = reaction_name.split('-TS-')
            reactants = '+'.join(reactants.split('-p-'))
            products = '+'.join(products.split('-p-'))
            
            reaction_name = '='.join([reactants,products])
            
            data = numpy.genfromtxt(Path + file_name,skip_header=2)
            xdata = data[0:18,0]
            ydata = data[0:18,4]
            ydata = numpy.log(data[0:18,4])
            
            dHf = numpy.mean(data[0:18,2])
            
            x0 = numpy.array([1.0,1.0,dHf])
            sigma = numpy.ones((ydata.size,))
                
            coeff = optimization.curve_fit(logkf, xdata, ydata)
            #coeff = optimization.curve_fit(logkf, xdata, ydata, x0, sigma, bounds = ([0, -10, -numpy.inf],[numpy.inf, 10, numpy.inf]))
            
            a = coeff[0]
            perr = numpy.sqrt(numpy.diag(coeff[1]))
            
            A = a[0]
            n = a[1]
            Ea = a[2]/4.184
            
            A = math.exp(a[0])
            n = a[1]
            Ea = a[2]/4.184
            
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

    with open(os.path.join(Path2final,'kinetics-barrierless.txt'), "w") as File:
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
list_of_barrierless_reactions = get_barrierless_reactions(pathCSV, fileCSV, pathTS)

calculate_rate_constants(pathTS,pathOutput,list_of_barrierless_reactions)

curvefit_arrhenius(pathOutput,pathCSV,fileCSV)


