# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 16:25:44 2023

@author: dundu
"""

import sys

# =============================================================================
# 
# =============================================================================

def Overwrite(i):

    # Load PARCELY environment input file
    ParcelyInput = open('PARCELY_V2/parcely_env_input_file.txt').read().splitlines()
    
    # Load variables table
    VarTable = open('InputTable.csv').read().splitlines()
    # Convert into list of lists
    for l, line in enumerate(VarTable):
        VarTable[l] = line.split(',')
                
    # First column is variable name, second column is variable line in input file
    # Column is given as input    
    NewFile = ParcelyInput[:4]
    
    # For each variable line
    for j in range(len(VarTable)):
        # Get previous value of variable   
        OldVal = ParcelyInput[4+j].split(':')[0]
        # Replace it with new value
        NewLine = ParcelyInput[4+j].replace(OldVal, 
                                        VarTable[j][2+i].replace(';', ',') +'\t\t')
        NewFile.append(NewLine)
    
    NewFile.append(ParcelyInput[-1])
        
    with open('PARCELY_V2/parcely_env_input_file.txt', 'w') as file:
        
        for line in NewFile:
            file.write(line + '\n')
    
        file.close()
    
    return

if __name__ == "__main__":
    i = int(sys.argv[1])
    Overwrite(i)