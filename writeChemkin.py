import csv
import os

def write_elements_block(pathOutput, fileChemkin, elements_list):
    with open(os.path.join(pathOutput,fileChemkin), "w") as File:
        File.writelines('ELEMENTS\n')
        for item in elements_list:
            File.writelines(item)
            File.writelines('\n')
        File.writelines('END\n')
        File.writelines('!\n')

def write_species_block(pathCSV, fileCSV, pathOutput, fileChemkin):
    ifile  = open(os.path.join(pathCSV,fileCSV), "rb")
    reader = csv.reader(ifile)

    rownum = 0
    species_list = []

    for row in reader:
        
        if rownum < 3:
            pass
        else:
            reactants,products = row[0].split('=')
            species = reactants.split('+') + products.split('+')
            species_list.extend(species)

        rownum = rownum + 1
            
    with open(os.path.join(pathOutput,fileChemkin), "a") as File:
        File.writelines('SPECIES\n')
        for item in sorted(list(set(species_list)),key=str.lower):
            #File.writelines('%-20s %-50s' % (item,'! '+ species_renamed[item]))
            File.writelines('%-20s' % (item))
            File.writelines('\n')
        File.writelines('END\n')
        File.writelines('!\n')                
  
    ifile.close()

def write_thermo_block(pathCSV, fileCSV, pathOutput, fileChemkin):
    T_low = 200
    T_common = 1000
    T_high = 5000
    
    ifile  = open(os.path.join(pathCSV,fileCSV), "rb")
    reader = csv.reader(ifile)

    rownum = 0
    species_list = []

    for row in reader:
        
        if rownum < 3:
            pass
        else:
            reactants,products = row[0].split('=')
            species = reactants.split('+') + products.split('+')
            species_list.extend(species)

        rownum = rownum + 1
           
    with open(os.path.join(pathOutput, fileChemkin), "a") as File:
        File.writelines('THERMO\n')
        File.writelines('%10.2f %10.2f %10.2f' % (T_low,T_common,T_high))
        File.writelines('\n')
        for item in sorted(list(set(species_list)),key=str.lower):
            f = open(os.path.join(pathOutput+"thermo-data-minima/"+item+"/",item+'.thermo'), "r")
            lines = f.readlines()
            File.writelines(lines[0])
            File.writelines(lines[1])
            File.writelines(lines[2])
            File.writelines(lines[3])
        File.writelines('END\n')
        File.writelines('!\n')                
      
    ifile.close()
    
def write_reactions_block(pathOutput, fileChemkin):
    f1 = open(os.path.join(pathOutput,'kinetics.txt'), "r")
    #f2 = open(os.path.join(pathOutput,'kinetics-barrierless.txt'), "r")  
    with open(os.path.join(pathOutput,fileChemkin), "a") as File:
        File.writelines('REACTIONS\n')
        File.writelines(f1.readlines()[2:])
        #File.writelines(f2.readlines()[2:])
        File.writelines('END\n')
        File.writelines('\n')    
            
###############################################################################
# Change input output paths and set output file names here
###############################################################################
base_directory = '/gpfs/group/umt/default/HMX/nitramine_gas_phase_mechanism/improvements-to-caltech-mechanism/final-gas-phase-mechanism/'

pathCSV = base_directory
pathB3LYP = base_directory + 'M062X-minima'
pathTS = base_directory + 'M062X-TS'
pathOutput = base_directory + 'Output-files/'

fileCSV = 'PSU-gas-phase-mechanism.csv' 
fileChemkin = 'PSU-gas-phase-mechanism.txt'                  
###############################################################################

elements_list = ['C','H', 'N', 'O', 'Ar']
write_elements_block(pathOutput, fileChemkin, elements_list)

write_species_block(pathCSV, fileCSV, pathOutput, fileChemkin)

write_thermo_block(pathCSV, fileCSV, pathOutput, fileChemkin)

write_reactions_block(pathOutput, fileChemkin)


