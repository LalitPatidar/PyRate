import os
import shutil
import glob
import numpy
import math
import scipy.optimize as optimization

def extractG4MP2(pathOutput,fileOutput,pathInput):
    """
    Extract G4MP2 log file data for all minima and calculates Heat of formation at 298 K
    Keeps updating the text file as new species are added
    """
    
    if not os.path.exists(pathOutput):
        os.mkdir(pathOutput)
    
    if not os.path.exists(pathOutput+fileOutput):
        header1 = "%-24s %-30s %-30s" %('Species', 'Enthalpy of formation (kcal/mol)', 'S_298 (cal/mol)') + '\n'
        header2 = "%-24s %-30s %-30s" %('-------', '--------------------------------', '---------------') + '\n' 
        with open(pathOutput+fileOutput, "w") as File:
            File.writelines(header1+header2)
    
    fdata = open(os.path.join(pathOutput,fileOutput), "r").read()
    
    lines_to_write = []
    flag_for_entropy = 0
    counter_for_entropy = 0
    flag_for_ZPE = 0

    E_o_carb=-37.794203              #Free Energy=-37.808748 triplet
    E_o_hydr= -0.502094              #Free Energy= -0.512748 doublet
    E_o_nitr=-54.532825              #Free Energy=-54.547860 quartet
    E_o_oxyg=-75.002483              #Free Energy=-75.017435 triplet

    # Sum of electronic and zero-point energy
    #E_o_carb= -37.747164              # Singlet 
    #E_o_hydr= -0.502094              # Doublet 
    #E_o_nitr= -54.439864              # Doublet
    #E_o_oxyg= -74.928184             # Singlet
           
    for file_name in os.listdir(pathInput):
        if file_name.endswith('.log') and (file_name.split('.log')[0]+ ' ') not in fdata:
            species_name = file_name.split('.log')[0]
            f = open(os.path.join(pathInput,file_name), "r")
            lines = f.readlines()
            
            N_carb = 0
            N_hydr = 0
            N_nitr = 0
            N_oxyg = 0
            
            count_begin = 0
            count_end = 0
            
            for line in lines:
                # Read elemental composition
                if line.startswith(' Redundant internal coordinates found in file.'):
                    count_begin = 1
                if line.startswith(' Recover connectivity data from disk.'):
                    count_end = 1
                if count_begin==1 and count_end==0 and line.startswith(' C'):
                    N_carb = N_carb + 1
                if count_begin==1 and count_end==0 and line.startswith(' H'):
                    N_hydr = N_hydr + 1
                if count_begin==1 and count_end==0 and line.startswith(' N'):
                    N_nitr = N_nitr + 1
                if count_begin==1 and count_end==0 and line.startswith(' O'):
                    N_oxyg = N_oxyg + 1
                             
                if line.startswith(' E(ZPE)='):
                    ZPE = line.split(' E(ZPE)=')[1].split('E(Thermal)=')[0]
                    ZPE = float(ZPE)

                    H_corr = line.split('E(ZPE)=')[1].split('E(Thermal)=')[1]
                    H_corr = float(H_corr)
                    
                if line.startswith(' G4MP2 Enthalpy='):
                    H = line.split('G4MP2 Enthalpy=')[1].split('G4MP2 Free Energy=')[0]
                    H = float(H)

                    G = line.split('G4MP2 Enthalpy=')[1].split('G4MP2 Free Energy=')[1]
                    G = float(G)
                                   
            H_f = (N_carb*(169.98-0.25) + N_hydr*(51.63-1.01) + N_nitr*(112.53-1.04) + N_oxyg*(58.99-1.04)) - 627.509*(N_carb*E_o_carb + N_hydr*E_o_hydr + N_nitr*E_o_nitr + N_oxyg*E_o_oxyg- H) 
            S = (H - G)*627.509*1000/298.15
            contents = "%-24s %20.3f %20.3f" %(file_name.split('.log')[0], H_f, S)
            lines_to_write.append(contents)
            
        flag_for_entropy = 0
        counter_for_entropy = 0

    #lines_to_write = sorted(list(set(lines_to_write)),key=str.lower)     
    text = "\n".join(lines_to_write) #+ "\n"

    with open(pathOutput+fileOutput, "a") as File:
        File.writelines(text)
        
def performThermochemistry(PathOutput,Path2B3LYP,Path2TS):
    #-----------------------------------------------------------------------------------------
    # Thermochemistry calculations minima
    for file_name in os.listdir(Path2B3LYP):
        if file_name.endswith('.chk'):
            species_name = file_name.split('.chk')[0]
            species_dir = PathOutput+"thermo-data-minima/" + file_name.split('.chk')[0]
                  
            if not os.path.exists(species_dir):
                os.makedirs(species_dir)
              
            shutil.copy(Path2B3LYP+'/'+file_name, species_dir)
            
            os.chdir(PathOutput+"thermo-data-minima")    
            os.chdir(species_name)        
            os.system('module load gaussian/g09d01')

            temperatures = [200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 850, 900, 950, 1000, \
            1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000]
            
            for temp in temperatures:
                output_file = str(temp) + '.txt'
                with open(output_file, "w") as File:
                    os.system('cd ..')
                    os.system('freqchk '+ file_name + ' -o ' + output_file + ' N ' + str(temp) +' 0 1 Y N') 

            lines_to_write = []
            flag_for_entropy = 0
            counter_for_entropy = 0
            
            for temperature_file in os.listdir("."):
                if temperature_file.endswith('.txt'):
                    f = open(temperature_file, "r")
                    lines = f.readlines()
                    for line in lines:
                        if line.startswith(' Thermal correction to Enthalpy= '):
                            H_corr = line.split('=')[1]
                            H_corr = float(H_corr)
                            
                        if line.startswith(' Thermal correction to Gibbs Free Energy='):
                            G_corr = line.split('=')[1]
                            G_corr = float(G_corr)
                            flag_for_entropy = 1

                        if counter_for_entropy==3:
                            Cv = float(line[40:50])
                            entropy = line[50:70]
                            S = (float(entropy))    
                        if flag_for_entropy == 1:
                            counter_for_entropy = counter_for_entropy + 1
                    
                    contents = "%-20d %15.6f %15.6f %15.3f %15.3f" %(int(temperature_file.split('.txt')[0]), H_corr, G_corr, S, Cv)
                    lines_to_write.append(contents)
                    
                flag_for_entropy = 0
                counter_for_entropy = 0

            lines_to_write = sorted(list(set(lines_to_write)),key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item))     
            text = "\n".join(lines_to_write) + "\n"
            
            header1 = "%-20s %15s %15s %15s %20s" %('Temperature', 'H_corr(Hartres) ', ' G_corr(Hartrees)', 'S(cal/mol-K)', 'Cv (Cal/Mol-K)') + '\n'
            header2 = "%-20s %15s %15s %15s %20s" %('-----------', '---------------', '----------------', '------------', '--------------') + '\n' 
            with open(species_name + '.data', "w") as File:
                File.writelines(header1+header2)
                File.writelines(text)                
    #-----------------------------------------------------------------------------------------
    # Thermochemistry calculations TS
    for file_name in os.listdir(Path2TS):
        if file_name.endswith('.chk'):
            species_name = file_name.split('.chk')[0]
            species_dir = PathOutput+"thermo-data-TS/" + file_name.split('.chk')[0]
                  
            if not os.path.exists(species_dir):
                os.makedirs(species_dir)
              
            shutil.copy(Path2TS+'/'+file_name, species_dir)
            
            os.chdir(PathOutput+"thermo-data-TS")    
            os.chdir(species_name)        
            os.system('module load gaussian/g09d01')

            temperatures = [200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 850, 900, 950, 1000, \
            1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000]
            
            for temp in temperatures:
                output_file = str(temp) + '.txt'
                with open(output_file, "w") as File:
                    os.system('cd ..')
                    os.system('freqchk '+ file_name + ' -o ' + output_file + ' N ' + str(temp) +' 0 1 Y N') 

            lines_to_write = []
            flag_for_entropy = 0
            counter_for_entropy = 0
            
            for temperature_file in os.listdir("."):
                if temperature_file.endswith('.txt'):
                    f = open(temperature_file, "r")
                    lines = f.readlines()
                    for line in lines:
                        if line.startswith(' Thermal correction to Enthalpy= '):
                            H_corr = line.split('=')[1]
                            H_corr = float(H_corr)
                            
                        if line.startswith(' Thermal correction to Gibbs Free Energy='):
                            G_corr = line.split('=')[1]
                            G_corr = float(G_corr)
                            flag_for_entropy = 1

                        if counter_for_entropy==3:
                            Cv = float(line[40:50])
                            entropy = line[50:70]
                            S = (float(entropy))    
                        if flag_for_entropy == 1:
                            counter_for_entropy = counter_for_entropy + 1
                    
                    contents = "%-20d %15.6f %15.6f %15.3f %15.3f" %(int(temperature_file.split('.txt')[0]), H_corr, G_corr, S, Cv)
                    lines_to_write.append(contents)
                    
                flag_for_entropy = 0
                counter_for_entropy = 0

            lines_to_write = sorted(list(set(lines_to_write)),key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item))     
            text = "\n".join(lines_to_write) + "\n"
            
            header1 = "%-20s %15s %15s %15s %20s" %('Temperature', 'H_corr(Hartres) ', ' G_corr(Hartrees)', 'S(cal/mol-K)', 'Cv (Cal/Mol-K)') + '\n'
            header2 = "%-20s %15s %15s %15s %20s" %('-----------', '---------------', '----------------', '------------', '--------------') + '\n' 
            with open(species_name + '.data', "w") as File:
                File.writelines(header1+header2)
                File.writelines(text)
                
def curvefit(PathOutput, fileHfData, Path2B3LYP):
    Ru = 1.9872
    T_ref = 298.15
    T_low = 200
    T_high = 5000
    T_int = 1000

    Hf_data = numpy.genfromtxt(PathOutput+fileHfData, dtype=None, skip_header=2,usecols=(0,1,2))

    def Cp_eval(a,x):
        return Ru*(a[0] + a[1]*x + a[2]*(x**2) + a[3]*(x**3) + a[4]*(x**4))
        
    def H_eval(a,x):
        return Ru*x*(a[0] + a[1]*x/2.0 + a[2]*(x**2)/3.0 + a[3]*(x**3)/4.0 + a[4]*(x**4)/5.0 + a[5]/x)
        
    def S_eval(a,x):
        return Ru*(a[0]*numpy.log(x) + a[1]*x + a[2]*(x**2)/2.0 + a[3]*(x**3)/3.0 + a[4]*(x**4)/4.0 + a[6])

    for file_name in os.listdir(Path2B3LYP):
        if file_name.endswith('.chk'):
            species_name = file_name.split('.chk')[0]
            
            species_dir = PathOutput+"thermo-data-minima/" + file_name.split('.chk')[0]
            os.chdir(PathOutput+"thermo-data-minima")    
            os.chdir(species_name) 
            
            data = numpy.genfromtxt(species_name + ".data",skip_header=2)
            xdata_low = data[0:29,0]
            ydata_low = (data[0:29,4] + Ru )/Ru
            S_B3LYP_low = (data[0:29,3])/Ru
            x0_low = numpy.array([1.0,1.0,1.0,1.0,1.0])
            sigma_low = numpy.ones((ydata_low.size,))
            sigma_low[-1] = 0.000001
           
            def Cp(x, a1, a2, a3, a4, a5):
                return a1 + a2*x + a3*(x**2) + a4*(x**3) + a5*(x**4)
                
            coeff_low = optimization.curve_fit(Cp, xdata_low, ydata_low, x0_low, sigma_low)
            
           
            a = coeff_low[0]                     
            
            for i in range(0,Hf_data.size):
                if Hf_data[i][0] == species_name:
                    H_f = float(Hf_data[i][1])
                    S_298 = float(Hf_data[i][2])
            
            a6 = 1000*H_f/Ru-a[0]*T_ref-a[1]*(T_ref**2)/2.0-a[2]*(T_ref**3)/3.0 -a[3]*(T_ref**4)/4.0-a[4]*(T_ref**5)/5.0
            
            a7 = S_298/Ru-a[0]*math.log(T_ref)-a[1]*T_ref-a[2]*(T_ref**2)/2.0-a[3]*(T_ref**3)/3.0-a[4]*(T_ref**4)/4.0
            
            a = numpy.append(a,[a6,a7])
            
            # High temperature limit fit
            
            xdata_high = data[28:,0]
            ydata_high = (data[28:,4] + Ru )/Ru
            S_B3LYP_high = (data[28:,3])/Ru
            x0_high = numpy.array([1.0,1.0,1.0,1.0,1.0])
            sigma_high = numpy.ones((ydata_high.size,))
            sigma_high[0] = 0.000001

            coeff_high = optimization.curve_fit(Cp, xdata_high, ydata_high, x0_high, sigma_high)

            b = coeff_high[0]                
                                       
            b6 = a[0]*(T_int) + a[1]*(T_int**2)/2.0 + a[2]*(T_int**3)/3.0 + a[3]*(T_int**4)/4.0 + a[4]*(T_int**5)/5.0 + a[5] \
                -b[0]*(T_int) - b[1]*(T_int**2)/2.0 - b[2]*(T_int**3)/3.0 - b[3]*(T_int**4)/4.0 - b[4]*(T_int**5)/5.0
                 
            b7 = a[0]*math.log(T_int) + a[1]*(T_int) + a[2]*(T_int**2)/2.0 + a[3]*(T_int**3)/3.0 + a[4]*(T_int**4)/4.0 + a[6] \
                -b[0]*math.log(T_int) - b[1]*(T_int) - b[2]*(T_int**2)/2.0 - b[3]*(T_int**3)/3.0 - b[4]*(T_int**4)/4.0
                
            b = numpy.append(b,[b6,b7])
            
            Cp_low = Cp_eval(a,xdata_low)
            Cp_high = Cp_eval(b,xdata_high)
            H_low = H_eval(a,xdata_low)
            H_high = H_eval(b,xdata_high)
            S_low = S_eval(a,xdata_low)
            S_high = S_eval(b,xdata_high)
                 
            with open(species_name + '.thermo', "w") as File:
                f = open(PathOutput+'log-file-data-minima.txt', "r")
                for line in f:
                    if line.startswith(species_name+' '):
                        species_info = line.rstrip()[0:44] + 'G' + '%10.2f%10.2f%8.2f' %(T_low,T_high,T_int)
                        File.writelines(species_info)
                        pos = 7
                        File.writelines('%*s' % (pos,'1'))
                        File.writelines('\n')
                        File.writelines('%15.8E%15.8E%15.8E%15.8E%15.8E%5s' % (b[0],b[1],b[2],b[3],b[4],'2'))
                        File.writelines('\n')
                        File.writelines('%15.8E%15.8E%15.8E%15.8E%15.8E%5s' % (b[5],b[6],a[0],a[1],a[2],'3'))
                        File.writelines('\n')
                        File.writelines('%15.8E%15.8E%15.8E%15.8E%20s' % (a[3],a[4],a[5],a[6],'4'))
                        File.writelines('\n')
                        File.writelines('\n')
                        File.writelines('Thermodynamic properties of '+ species_name + ' at 298.15K using the G4(MP2) compound method: \n')
                        File.writelines('H_f   = ' + str(H_f) +   ' kcal/mol\n')
                        File.writelines('S_298 = ' + str(S_298) + ' cal/mol-K\n')
                        File.writelines('Cv and S at other temperatures were computed using B3LYP/6311++G(d,p)\n')
                        File.writelines('\n')
                        File.writelines('%8s %15s %15s %15s %15s %15s' %('Temp.','Cp','H','S','Cp(B3LYP)','S(B3LYP)\n'))
                        File.writelines('%8s %15s %15s %15s %15s %15s' %('K','cal/mol-K','cal/mol','cal/mol-K','cal/mol-K','cal/mol-K\n'))
                        
            f2 = open(species_name + '.thermo', "a")        
            numpy.savetxt(f2,numpy.c_[xdata_low, Cp_low, H_low, S_low, ydata_low*Ru, S_B3LYP_low*Ru], fmt = "%8.1f %15.3E %15.3E %15.3E %15.3E %15.3E")
            numpy.savetxt(f2,numpy.c_[xdata_high, Cp_high, H_high, S_high, ydata_high*Ru, S_B3LYP_high*Ru], fmt = "%8.1f %15.3E %15.3E %15.3E %15.3E %15.3E")
                    
def extractSCF(pathOutput, pathB3LYP,pathTS):
    # Minima
    header1 = "%-40s %30s %30s" %('Species', 'SCF-Energy(Hartres)', 'Zero-point correction (Hartrees)') + '\n'
    header2 = "%-40s %30s %30s" %('-------', '-------------------', '--------------------------------') + '\n' 
    with open(os.path.join(pathOutput,'minima-SCF.txt'), "w") as File:
        File.writelines(header1+header2)

    lines_to_write = []
           
    for file_name in os.listdir(pathB3LYP):
        if file_name.endswith('.log'):
            f = open(os.path.join(pathB3LYP,file_name), "r")
            lines = f.readlines()
            for line in lines:
                if line.startswith(' SCF Done:'):
                    scf = float((line.strip().split()[4]))
                
                if line.startswith(' Zero-point correction='):
                    zpe = float((line.strip().split()[2]))
                    
            contents = "%-40s %30.6f %30.6f" %(file_name.split('.log')[0], scf, zpe)
            lines_to_write.append(contents)
              
    lines_to_write = sorted(list(set(lines_to_write)),key=str.lower)     
    text = "\n".join(lines_to_write) + "\n"
    with open(os.path.join(pathOutput,'minima-SCF.txt'), "a") as File:
        File.writelines(text)

    # TS
    header1 = "%-40s %30s %30s %30s %30s %30s" %('Species', 'SCF-Energy(Hartres)', 'Zero-point correction (Hartrees)','Reduced mass (AMU)','Force constant (mDyne/A)','Imaginary freq. (cm-1)') + '\n'
    header2 = "%-40s %30s %30s %30s %30s %30s" %('-------', '-------------------', '--------------------------------','------------------','------------------------', '---------------------') + '\n' 
    with open(os.path.join(pathOutput,'TS-SCF.txt'), "w") as File:
        File.writelines(header1+header2)

    lines_to_write = []
           
    for file_name in os.listdir(pathTS):
        if file_name.endswith('.log'):
            f = open(os.path.join(pathTS,file_name), "r")
            lines = f.readlines()
            
            flag = 0
            count = 0
            
            for line in lines:
                if line.startswith(' SCF Done:'):
                    scf = float((line.strip().split()[4]))
                    
                if line.startswith(' Frequencies --'):
                    if float((line.strip().split()[2])) < 0:
                        img_freq = float((line.strip().split()[2]))
                        flag = 1
                
                if flag == 1:
                    count = count + 1
                    
                if count == 2:
                    mu = float((line.strip().split()[3]))
                
                if count == 3:
                    fc = float((line.strip().split()[3]))     
                
                if line.startswith(' Zero-point correction='):
                    zpe = float((line.strip().split()[2]))    
                    
            contents = "%-40s %30.6f %30.6f %30.6f %30.6f %30.6f" %(file_name.split('.log')[0], scf, zpe, mu, fc, img_freq)
            lines_to_write.append(contents)
              
    lines_to_write = sorted(list(set(lines_to_write)),key=str.lower)     
    text = "\n".join(lines_to_write) + "\n"
    with open(os.path.join(pathOutput,'TS-SCF.txt'), "a") as File:
        File.writelines(text)    
            
###############################################################################
# Change input output paths and set output file names here
###############################################################################

base_directory = '/gpfs/group/umt/default/HMX/nitramine_gas_phase_mechanism/improvements-to-caltech-mechanism/final-gas-phase-mechanism/'

pathB3LYP = base_directory + 'M062X-minima'
pathTS = base_directory + 'M062X-TS'
pathG4MP2 = base_directory + 'G4MP2-data'
pathOutput = base_directory + 'Output-files/'

# Prefer not to change this. Otherwise subsequent programs, which uses these files, will need to be changed and give the correct file names                   
fileHfData = 'G4MP2-Hf-data.txt'
###############################################################################
extractG4MP2(pathOutput,fileHfData, pathG4MP2)

performThermochemistry(pathOutput, pathB3LYP,pathTS)

curvefit(pathOutput, fileHfData, pathB3LYP)

extractSCF(pathOutput, pathB3LYP,pathTS)


