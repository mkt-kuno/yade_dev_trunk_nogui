''' A Helper class to setup the compilation of the Yade-OpenFOAM coupling module and OpenFOAM associated solvers. 
(c) 2024 Deepak Kunhappan, deepak.kn1990@gmail.com 
'''

import os
import argparse 

#parser = argparse.ArgumentParser() 
#parser.add_argument("build", default = "build", nargs = '?', help = "Compiles the coupling and solvers from scratch") 
#parser.add_argument("--clean", help = "Cleans existing build")
#

#args = parser.parse_args()

rootDir = os.getcwd()
print("Preparing Yade-OpenFOAM compilation")
supportedFoamVersions = ['v2312', 'v1906', '6'] #11 is not supported yet, maybe 2212 is also supported?

currentFoamVersion = os.environ['WM_PROJECT_VERSION']

# check if we are dealing with the foundation or the .com version of OpenFOAM. For .com versions, $WM_PROJECT_DIR
# returns a string starting with 'v' followd by the year and month of release. For the foundation versions, the
# release is denoted by an integer.
# TODO: handle dev versions of OF-foundation (?)
isFoundationVersion = False if currentFoamVersion[0] == 'v' else True
isVersion2312 = True if currentFoamVersion == 'v2312' else False

if (currentFoamVersion == ''): 
    print("Could not detect OpenFOAM version, have you sourced bashrc from $FOAM_SRC/etc/bashrc ?")
    exit() 
else:
    if not currentFoamVersion in supportedFoamVersions:
        print("Error : The OpenFOAM version", currentFoamVersion, " is not supported, but we will try compiling as if it was v2312 (or 6 if a foundation version is found)")
        isVersion2312 = !isFoundationVersion  # pretend we use 2312, in some cases it will work
    else: print("OpenFoam coupling will be compiled for OpenFOAM version", currentFoamVersion)



# modify the options make file for FoamYade if a foundation release is being used: 
if isFoundationVersion:
    os.chdir(rootDir + '/FoamYade/Make/')    
    with open('options', 'r') as makeOptions:        
        allOptions = makeOptions.readlines()
    allOptions.insert(3, '    -DFOUNDATION \ \n')        
    makeOptions = open('options', 'w')
    for line in allOptions: 
        makeOptions.write(line)
        #makeOptions.write(allOptions)
    # return to rootDir 
    os.chdir(rootDir)

# compile sources .. 
os.chdir(rootDir)

# compile commYade
print("Compiling communications module")
os.chdir(rootDir + '/FoamYade/commYade')
os.system('wclean')
os.system('wmake')

# mesh tree
print("Compiling mesh tree module")
os.chdir(rootDir + '/FoamYade/meshtree')
os.system('wclean')
os.system('wmake')

# FoamYade coupling
print("Compiling Yade-OpenFOAM coupling module")
os.chdir(rootDir + '/FoamYade')
os.system('wclean')
os.system('wmake')

# Solvers... 
print("Compiling icoFoamYade solver")
os.system('wclean')
os.chdir(rootDir + '/Solvers/icoFoamYade')
os.system('wclean')
os.system('wmake')

if isVersion2312:
    print("Compiling pimpleFoamYade solver")
    os.chdir(rootDir + '/Solvers/pimpleFoamYadev2312')
    os.system('wclean')
    os.system('wmake')
else :
    print("Compiling pimpleFoamYade solver")
    os.system('wclean')
    os.chdir(rootDir + '/Solvers/pimpleFoamYade')
    
