#!/bin/bash

#Input

appList=(advectorVoF advectorVoF) 
schemeList=(plicRDF plicRDFNoInterpol) 
meshList=(hex tri poly) 
CoList=(0.1 0.2 0.5)
#CoList=( 0.5)

#End of input

#VoFLib_ROOT_DIR=PATH TO VoFLib 

#Check if VoFLib_ROOT_DIR is set
if [ -z "$VoFLib_ROOT_DIR" ];
then
    echo " "
    echo "Warning: "
    echo "Please set and export VoFLib_ROOT_DIR to the "
    echo "root directory of the isoAdvector source code."
    echo "Aborting "
    echo " "
else

# Source utilities
. $VoFLib_ROOT_DIR/bin/dhiFoamTools

curDir=$PWD

#Looping through all mesh, solver, and CFL number combinations and generating cases
for nn in ${!meshList[*]}
do
    meshType=${meshList[$nn]}
    for mm in ${!appList[*]}
    do
        application=${appList[$mm]}
        scheme=${schemeList[$mm]}

        #Case location
        series=$PWD/$scheme/$meshType

        if [ "$meshType" = "hex" ];
        then
            #Domain  dimensions
            L=5
            H=3
		#Vertical velocity component
		Uz=0.5
		NzList=(15 30 60 120 240 480)
		NxList=(25 50 100 200 400 800)
		#NzList=(15 30 60)
		#NxList=(25 50 100)
	else
		NzList=(15 30 60 120 240)
		NxList=(75 150 300 600 1200)
		#NzList=(15 30 60)
		#NxList=(75 150 300)
		Uz=0.0
        fi

        mkdir --parents $series

        for n in ${!NzList[*]}
        do
            for m in ${!CoList[*]}
            do
                Co=${CoList[$m]}
                if [ "$meshType" = "hex" ];
                then
                    caseName=N${NzList[$n]}Co${Co}
                else
                    caseName=${NzList[$n]}Co${Co}
                fi
                caseDir=$series/$caseName
                echo $caseDir
                cp -r baseCase $caseDir

                foamParmSet maxAlphaCo "$Co" $caseDir/system/controlDict
                foamParmSet application "$application" $caseDir/system/controlDict

		mkdir $caseDir/logs
    
        if [ "$scheme" = "plicRDF" ];
        then
    		sed -i "s/INTERPOL/true/g" $caseDir/system/fvSolution
        else
            sed -i "s/INTERPOL/false/g" $caseDir/system/fvSolution
        fi
                if [ "$meshType" = "hex" ];
                then
                    nx=${NxList[$n]}
                    nz=${NzList[$n]}
                    foamParmSet L "$L" $caseDir/constant/polyMesh/blockMeshDict
                    foamParmSet H "$H" $caseDir/constant/polyMesh/blockMeshDict
                    foamParmSet nx "$nx" $caseDir/constant/polyMesh/blockMeshDict
                    foamParmSet nz "$nz" $caseDir/constant/polyMesh/blockMeshDict
                    foamParmSet 'internalField' "uniform (1 0 .5)" $caseDir/0.org/U
                    cp -r $caseDir/0.org $caseDir/0
		    blockMesh -case $caseDir > $caseDir/logs/blockMesh.log 2>&1
                else
                    #cp $triMeshDir/polyMesh_${NzList[$n]}/* $caseDir/constant/polyMesh/
                    foamParmSet 'internalField' "uniform (1 0 0)" $caseDir/0.org/U
		    nx=${NxList[$n]}
		    nz=${NzList[$n]}
		    #foamParmSet nx "$nx" $caseDir/triSquare.geo
         	    #foamParmSet nz "$nz" $caseDir/triSquare.geo
		    #./ofset2 replaceNx "$nx" $caseDir/triSquare.geo
			sed  -i "s/replaceNx/$nx/" $caseDir/triSquare.geo
			sed  -i "s/replaceNz/$nz/" $caseDir/triSquare.geo

		    cd $caseDir
		    gmsh -3 triSquare.geo > $caseDir/logs/gmsh.log 2>&1
		    gmshToFoam triSquare.msh > $caseDir/logs/gmshToFoam.log 2>&1
		    changeDictionary > $caseDir/logs/changeDictionary.log 2>&1
		    cd $curDir
                    cp -r $caseDir/0.org $caseDir/0
                    if [ "$meshType" = "poly" ];
                    then
                        #mkdir $caseDir/logs
                        #Convert from tet to poly mesh
                        polyDualMesh -case $caseDir -overwrite 160 > $caseDir/logs/polyDualMesh.log 2>&1
                        #Remove backmost  part of cells
                        topoSet -case $caseDir > $caseDir/logs/topoSet.log 2>&1
                        subsetMesh -case $caseDir -overwrite c0 -patch front > $caseDir/logs/subsetMesh.log 2>&1
                        transformPoints -case $caseDir -scale '(1 2 1)' > $caseDir/logs/transformPoints.log 2>&1 # thickness should be one
                        rm -rf *.obj
                    fi
                fi
            done
        done
    done
done

#echo " "
#echo "Cases generated. To run cases type"
#echo " "
#echo "foamRunCasesIn " $schemeList
#echo " "

fi
