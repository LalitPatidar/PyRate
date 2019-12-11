import csv
import os

###############################################################################
speciesList = []
transitionStatesList = []
reactionsList = []
symmetryList = []
###############################################################################

def readInputCSV(path):
    """
    Read input csv file where all reactions and their corresponding forward and backward symmetry factors are listed. See 'RDX-all-reactions.csv' for example and format
    Creates a list of all species, a list of all reactions and a list of names of all transition state files. Also creates a list of symmetry factors.
    """
    global speciesList, transitionStatesList, reactionsList, symmetryList
    
    with open(path,'r') as ifile:
        reader = csv.reader(ifile)
        
        rownum = 0
        for row in reader:
            if rownum < 3:
                pass
            else:
                reactants,products = row[0].split('=')
                species = reactants.split('+') + products.split('+')
                for item in species:
                    if item not in speciesList:
                        speciesList.append(item)
                TS = "-p-".join(reactants.split('+')) + '-TS-' + "-p-".join(products.split('+'))
                transitionStatesList.append(TS)
                reactionsList.append(row[0])
                symmetryList.append([row[0],row[1],row[2]])        
                   
            rownum = rownum + 1
        
        longNames = []
        for species in speciesList:
            if len(species) >=15:
                longNames.append(species)
                
        if len(longNames)>0:
            print 'The names of following species is greater than 15 characters'
            print 'Consider reducing them since chemkin format does not allow long species names' 
            print '\n'.join(longNames)

def LogFilesList(path):
    """
    Read number of log files present in a folder so that missing files or extra files can be identified and removed for consistency
    """
    filesList = []
    
    for file_name in os.listdir(path):
        if file_name.endswith('.log'):
            filesList.append(file_name.split('.log')[0])    
    
    return filesList

def showStats(pathMinima, pathTS):
    """
    Displays the number of species from the reactions list and from files; displays the missing or extra files in species folder
    Displays the number of reactions from reactions list and from files; displays the missing of extra files in the transition-states folder
    """
    global speciesList, transitionStatesList, reactionsList, symmetryList
    
    minimaLogFilesList = LogFilesList(pathMinima)
    TSLogFilesList = LogFilesList(pathTS)
    
    print 'Number of species from reactions  = %4d' % len(speciesList)        
    print 'Number of species from files      = %4d' % len(minimaLogFilesList)
    print 'Number of TS from reactions       = %4d' % len(transitionStatesList)
    print 'Number of TS from files           = %4d' % len(TSLogFilesList)
    print 'Number of unique reactions        = %4d' % len(list(set(reactionsList)))

    print '\nNumber of species files missing   = %4d' % (len(speciesList) - len(minimaLogFilesList))
    print 'Number of TS files missing        = %4d' % (len(reactionsList) - len(TSLogFilesList))

    print '\nFiles missing from the minima folder'
    for item in list(set(speciesList)):
        if item not in minimaLogFilesList:
            print item

    print '\nFiles missing from the TS folder'
    for item in transitionStatesList:
        if item not in TSLogFilesList:
            print item

    print '\nExtra Files in the minima folder'
    for item in minimaLogFilesList:
        if item not in list(set(speciesList)):
            print item

    print '\nExtra Files in the TS folder'
    for item in TSLogFilesList:
        if item not in transitionStatesList:
            print item

def logFileData_minima(path,file,pathMinima):
    """
    Extract log file data for all minima and store it in a text file
    Keeps updating the text file as new species are added
    """
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    if not os.path.exists(path+file):
        header1 = '%-24s%5s%5s%5s%5s %15s %15s %15s %15s %15s %25s %25s' \
                   %('Species', '#C','#H','#N','#O', 'SCF(Hartrees)', 'ZPE(Hartrees)', 'H(Hartres)', 'G(Hartrees)', 'S(cal/mol-K)', 'Molar-Volume(cm**3/mol)', 'Molecular-Radius (a0)') + '\n'
        header2 = '%-24s%3s%3s%3s%3s %15s %15s %15s %15s %15s %25s %25s' \
                   %('-------', '-----','-----','-----','-----','-------------', '--------------', '---------- ', '------------','-------------','-------------------------', '-------------------------') + '\n' 
        
        with open(path+file, 'w') as ofile:
            ofile.writelines(header1+header2)                        
    
    fdata = open(os.path.join(path,path+file), "r").read()
    
    lines_to_write = []
    flag_for_entropy = 0
    counter_for_entropy = 0
    
    for file_name in os.listdir(pathMinima):
        if file_name.endswith('.log') and (file_name.split('.log')[0]+ ' ') not in fdata:
            f = open(os.path.join(pathMinima,file_name), "r")
            lines = f.readlines()
            
            number_of_carbon = 0
            number_of_hydrogen = 0
            number_of_nitrogen = 0
            number_of_oxygen = 0
            
            count_begin = 0
            count_end = 0
            
            for line in lines:
                # Read elemental composition
                if line.startswith(' Redundant internal coordinates found in file.'):
                    count_begin = 1
                if line.startswith(' Recover connectivity data from disk.'):
                    count_end = 1
                if count_begin==1 and count_end==0 and line.startswith(' C'):
                    number_of_carbon = number_of_carbon + 1
                if count_begin==1 and count_end==0 and line.startswith(' H'):
                    number_of_hydrogen = number_of_hydrogen + 1
                if count_begin==1 and count_end==0 and line.startswith(' N'):
                    number_of_nitrogen = number_of_nitrogen + 1
                if count_begin==1 and count_end==0 and line.startswith(' O'):
                    number_of_oxygen = number_of_oxygen + 1
                
                # Read SCF energy and ZPE
                if line.startswith(' SCF Done:'):
                    scf = float((line.strip().split()[4]))
            
                if line.startswith(' Zero-point correction='):
                    zpe = float((line.strip().split()[2]))
                    
                # Read enthalpy H, gibbs free energy G and entropy S
                if line.startswith(' Sum of electronic and thermal Enthalpies='):
                    H = line.split('=')[1]
                    H = float(H)
                if line.startswith(' Sum of electronic and thermal Free Energies='):
                    G = line.split('=')[1]
                    G = float(G)
                    flag_for_entropy = 1
                if counter_for_entropy==4:
                    entropy = line[50:70]
                    S = (float(entropy))
                if flag_for_entropy == 1:
                    counter_for_entropy = counter_for_entropy + 1
                    
                # Read molar volume and radius
                if line.startswith(' Molar volume ='):
                    vol = line.split('bohr**3/mol')[1].rstrip()
                if line.startswith(' Recommended a0'):
                    a0 = line.split('=')[1].rstrip()

            contents = '%-24s%-2s%3s%-2s%3s%-2s%3s%-2s%3s %15.6f %15.6f %15.6f %15.6f %15.3f %25s %25s' \
                        %(file_name.split('.log')[0], 'C',number_of_carbon,'H',number_of_hydrogen,'N',number_of_nitrogen,'O',number_of_oxygen, scf, zpe, H, G, S, vol, a0)
            lines_to_write.append(contents)
        
        flag_for_entropy = 0
        counter_for_entropy = 0
              
    #lines_to_write = sorted(list(set(lines_to_write)),key=str.lower)     
    text = "\n".join(lines_to_write) + "\n"
    with open(path+file, "a") as File:
        File.writelines(text)
        
def logFileData_TS(path,file,pathTS):
    """
    Extract log file data for all TS and store it in a text file
    Keeps updating the text file as new reactions are added
    """
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    if not os.path.exists(path+file):
        header1 = '%-50s%5s%5s%5s%5s %15s %15s %15s %15s %15s %20s %20s %20s' \
                   %('Species', '#C','#H','#N','#O', 'SCF(Hartrees)', 'ZPE(Hartrees)', 'H(Hartres)', 'G(Hartrees)', 'S(cal/mol-K)', 'Reduced mass(AMU)','Force const(mDyne/A)','Imag. freq.(cm-1)') + '\n'
        header2 = '%-50s%3s%3s%3s%3s %15s %15s %15s %15s %15s %20s %20s %20s' \
                   %('-------', '-----','-----','-----','-----','-------------', '--------------', '---------- ', '------------','-------------','-----------------','--------------------', '-----------------') + '\n' 
        
        with open(path+file, 'w') as ofile:
            ofile.writelines(header1+header2)                        
    
    fdata = open(os.path.join(path,path+file), "r").read()
    
    lines_to_write = []
    flag_for_entropy = 0
    counter_for_entropy = 0
    
    for file_name in os.listdir(pathTS):
        if file_name.endswith('.log') and (file_name.split('.log')[0]+ ' ') not in fdata:
            f = open(os.path.join(pathTS,file_name), "r")
            lines = f.readlines()
            
            number_of_carbon = 0
            number_of_hydrogen = 0
            number_of_nitrogen = 0
            number_of_oxygen = 0
            
            count_begin = 0
            count_end = 0
            
            flag = 0
            count = 0
            
            for line in lines:
                # Read elemental composition
                if line.startswith(' Redundant internal coordinates found in file.'):
                    count_begin = 1
                if line.startswith(' Recover connectivity data from disk.'):
                    count_end = 1
                if count_begin==1 and count_end==0 and line.startswith(' C'):
                    number_of_carbon = number_of_carbon + 1
                if count_begin==1 and count_end==0 and line.startswith(' H'):
                    number_of_hydrogen = number_of_hydrogen + 1
                if count_begin==1 and count_end==0 and line.startswith(' N'):
                    number_of_nitrogen = number_of_nitrogen + 1
                if count_begin==1 and count_end==0 and line.startswith(' O'):
                    number_of_oxygen = number_of_oxygen + 1
                
                # Read SCF energy and ZPE
                if line.startswith(' SCF Done:'):
                    scf = float((line.strip().split()[4]))
            
                if line.startswith(' Zero-point correction='):
                    zpe = float((line.strip().split()[2]))
                    
                # Read enthalpy H, gibbs free energy G and entropy S
                if line.startswith(' Sum of electronic and thermal Enthalpies='):
                    H = line.split('=')[1]
                    H = float(H)
                if line.startswith(' Sum of electronic and thermal Free Energies='):
                    G = line.split('=')[1]
                    G = float(G)
                    flag_for_entropy = 1
                if counter_for_entropy==4:
                    entropy = line[50:70]
                    S = (float(entropy))
                if flag_for_entropy == 1:
                    counter_for_entropy = counter_for_entropy + 1
                    
                # Read molar volume and radius
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
            
            contents = '%-50s%-2s%3s%-2s%3s%-2s%3s%-2s%3s %15.6f %15.6f %15.6f %15.6f %15.3f %20.6f %20.6f %20.6f' \
                        %(file_name.split('.log')[0], 'C',number_of_carbon,'H',number_of_hydrogen,'N',number_of_nitrogen,'O',number_of_oxygen, scf, zpe, H, G, S, mu, fc, img_freq)
            lines_to_write.append(contents)
        
        flag_for_entropy = 0
        counter_for_entropy = 0
              
    #lines_to_write = sorted(list(set(lines_to_write)),key=str.lower)     
    text = "\n".join(lines_to_write) + "\n"
    with open(path+file, "a") as File:
        File.writelines(text)
    
    
def writeCSV(path,file, fileMinima, fileTS):
    """
    Create output csv if it doesn't exist or update it if already exist
    """
    if not os.path.exists(path):
        os.mkdir(path)
    
    if not os.path.exists(path+file):
        with open(path+file, 'w') as ofile:
            writer = csv.writer(ofile)

            row0 = ['!','','','','','cm-1','','cal/mol-K','cal/mol-K','cal/mol','cal/mol','cal/mol','cal/mol','cal/mol','cal/mol']

            row1a = ['REACTIONS','for.sym.','conv.f.','back.sym.','conv.b.','Freq','Wig.Tunn.','DSf','DSb','DHf','DHb','DGf','DGb','DHr','DGr']
            row1b = ['H_A','H_B','H_TS','H_C','H_D','H_E','G_A','G_B','G_TS','G_C','G_D','G_E','S_A','S_B','S_TS','S_C','S_D','S_E']

            row1 = row1a + row1b
            writer.writerow(row0)
            writer.writerow(row1)
            writer.writerow(['NEW_FORMAT','','','','','','','','','','','','','',''])
            
    minima_data = open(path+fileMinima, "r").readlines()
    TS_data = open(path+fileTS, "r").readlines()
    
    with open(path+file, 'a') as ofile:
        writer = csv.writer(ofile)
    
        for reaction in reactionsList:
            reactants,products = reaction.split('=')
            TS = "-p-".join(reactants.split('+')) + '-TS-' + "-p-".join(products.split('+'))
            
            if len(reactants.split('+'))==1:
                A = reactants
                conv_f = 1
                
            if len(reactants.split('+'))==2:
                A,B = reactants.split('+')
                conv_f = 1000
                        
            if len(products.split('+'))==1:
                C = products
                conv_b = 1
                              
            if len(products.split('+'))==2:
                C,D = products.split('+')
                conv_b = 1000
                                
            if len(products.split('+'))==3:
                C,D,E = products.split('+')
                conv_b = 1000000
                
            if len(products.split('+'))==4:
                C,D,E,F = products.split('+')
                conv_b = 1000000000
                

            for line in minima_data:
                    if line.startswith(A+' '):
                        H_A = float(line.strip().split()[8])
                        G_A = float(line.strip().split()[9])
                        S_A = float(line.strip().split()[10])
                    
                    if conv_f == 1:
                        H_B = 0
                        G_B = 0
                        S_B = 0
                        
                    if conv_f == 1000:
                        if line.startswith(B+' '):
                            H_B = float(line.strip().split()[8])
                            G_B = float(line.strip().split()[9])
                            S_B = float(line.strip().split()[10])
                        
                    if line.startswith(C+' '):
                        H_C = float(line.strip().split()[8])
                        G_C = float(line.strip().split()[9])
                        S_C = float(line.strip().split()[10])
                        
                    if conv_b == 1:
                        H_D = 0
                        G_D = 0
                        S_D = 0
                        H_E = 0
                        G_E = 0
                        S_E = 0
                    
                    if conv_b == 1000:    
                        if line.startswith(D+' '):
                            H_D = float(line.strip().split()[8])
                            G_D = float(line.strip().split()[9])
                            S_D = float(line.strip().split()[10])
                            
                        H_E = 0
                        G_E = 0
                        S_E = 0
                        
                    if conv_b == 1000000:
                        if line.startswith(D+' '):
                            H_D = float(line.strip().split()[8])
                            G_D = float(line.strip().split()[9])
                            S_D = float(line.strip().split()[10])
                            
                        if line.startswith(E+' '):
                            H_E = float(line.strip().split()[8])
                            G_E = float(line.strip().split()[9])
                            S_E = float(line.strip().split()[10])
                                    
                    H_F = 0
                    G_F = 0
                    S_F = 0    
                    if conv_b == 1000000000:
                        if line.startswith(D+' '):
                            H_D = float(line.strip().split()[8])
                            G_D = float(line.strip().split()[9])
                            S_D = float(line.strip().split()[10])
                            
                        if line.startswith(E+' '):
                            H_E = float(line.strip().split()[8])
                            G_E = float(line.strip().split()[9])
                            S_E = float(line.strip().split()[10])
                            
                        if line.startswith(F+' '):
                            H_F = float(line.strip().split()[8])
                            G_F = float(line.strip().split()[9])
                            S_F = float(line.strip().split()[10])
                        
                        H_E = H_E + H_F
                        G_E = G_E + G_F
                        S_E = S_E + S_F
                    
            if conv_f == 1:
                if conv_b == 1:
                    if (G_A < G_C):
                        H_TS = H_C
                        G_TS = G_C
                        S_TS = S_C
                    else:
                        H_TS = H_A
                        G_TS = G_A
                        S_TS = S_A
                elif conv_b == 1000:
                    if (G_A < G_C + G_D):
                        H_TS = H_C + H_D
                        G_TS = G_C + G_D
                        S_TS = S_C + S_D
                    else:
                        H_TS = H_A
                        G_TS = G_A
                        S_TS = S_A
                else:
                     if (G_A < G_C + G_D + G_E):
                        H_TS = H_C + H_D + H_E
                        G_TS = G_C + G_D + G_E
                        S_TS = S_C + S_D + S_E
                     else:
                        H_TS = H_A
                        G_TS = G_A
                        S_TS = S_A
            
            if conv_f == 1000:
                if conv_b == 1:
                    if (G_A + G_B < G_C):
                        H_TS = H_C
                        G_TS = G_C
                        S_TS = S_C
                    else:
                        H_TS = H_A + H_B
                        G_TS = G_A + G_B
                        S_TS = S_A + S_B
                elif conv_b == 1000:
                    if (G_A + G_B < G_C + G_D):
                        H_TS = H_C + H_D
                        G_TS = G_C + G_D
                        S_TS = S_C + S_D
                    else:
                        H_TS = H_A + H_B
                        G_TS = G_A + G_B
                        S_TS = S_A + S_B
                else:
                     if (G_A + G_B < G_C + G_D + G_E):
                        H_TS = H_C + H_D + H_E
                        G_TS = G_C + G_D + G_E
                        S_TS = S_C + S_D + S_E
                     else:
                        H_TS = H_A + H_B
                        G_TS = G_A + G_B
                        S_TS = S_A + S_B        
                    
            Freq_TS = float(1)
            Wigner = 1 + ((0.0048366*Freq_TS)**2)/24        
            for line in TS_data:
                if line.startswith(TS):
                    H_TS = float(line.strip().split()[8])
                    G_TS = float(line.strip().split()[9])
                    S_TS = float(line.strip().split()[10])
                    Freq_TS = -1*float(line.strip().split()[13])
                    Wigner = 1 + ((0.0048366*Freq_TS)**2)/24   
                    
            dSf = S_TS - S_A - S_B
            dSb = S_TS - S_C - S_D - S_E
            
            dHf = (H_TS - H_A - H_B)*627.509*1000
            dHb = (H_TS - H_C - H_D - H_E)*627.509*1000
            
            dGf = (G_TS - G_A - G_B)*627.509*1000
            dGb = (G_TS - G_C - G_D - G_E)*627.509*1000
            
            dHr = (H_C + H_D + H_E - H_A - H_B)*627.509*1000
            dGr = (G_C + G_D + G_E - G_A - G_B)*627.509*1000
            
            for sym in symmetryList:
                if sym[0]== reaction:
                    for_sym = float(sym[1])
                    back_sym = float(sym[2])
                   
            row = [reaction,'%-8.2E' %(for_sym), '%-8.2E' %(conv_f), '%-8.2E' %(back_sym), '%-8.2E' %(conv_b), '%-8.2E' %(Freq_TS), '%-8.2E' %(Wigner), \
                    '%-10.4E' %(dSf),'%-10.4E' %(dSb), '%-10.4E' %(dHf), '%-10.4E' %(dHb), '%-10.4E' %(dGf),'%-10.4E' %(dGb), '%-10.4E' %(dHr), '%-10.4E' %(dGr), \
                    H_A, H_B, H_TS, H_C, H_D, H_E, G_A, G_B, G_TS, G_C, G_D, G_E, S_A, S_B, S_TS, S_C, S_D, S_E]
        
            writer.writerow(row)    
    

###############################################################################
# Change input output paths and set output file names here
###############################################################################

base_directory = '/gpfs/group/umt/default/HMX/nitramine_gas_phase_mechanism/improvements-to-caltech-mechanism/final-gas-phase-mechanism/'

pathCSV = base_directory + 'PSU-gas-phase-mechanism.csv'
pathMinima = base_directory + 'M062X-minima'
pathTS = base_directory + 'M062X-TS'
pathOutput = base_directory + 'Output-files/'

fileOutputCSV = 'PSU-gas-phase-mechanism-data.csv'

# Prefer not to change this. Otherwise subsequent programs, which uses these files, will need to be changed and give the correct file names                      
fileMinimaLogData = 'log-file-data-minima.txt'
fileTSLogData = 'log-file-data-TS.txt'
###############################################################################
               
readInputCSV(pathCSV)

showStats(pathMinima, pathTS)

logFileData_minima(pathOutput,fileMinimaLogData,pathMinima)

logFileData_TS(pathOutput,fileTSLogData,pathTS)

writeCSV(pathOutput,fileOutputCSV, fileMinimaLogData, fileTSLogData)
