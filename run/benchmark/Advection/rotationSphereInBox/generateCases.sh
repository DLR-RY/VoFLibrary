#!/bin/bash

#Input

#appList=(interFlow interFlow interFlow interFlow)
appList=(advectorVoF  )
#schemeList=(interFlow MULES HRIC CICSAM)
#schemeList=( perfectRDFPoints  RDFadvect isoInverseDistance )
schemeList=( plicRDF ) # isoInverseDistance RDFadvect
meshList=(poly hex tri) # hex tri poly
dtList=(0.01 0.005 0.0025 0.001  )
#dtList=( 0.0025 )
CoList=(  0.5) #0.1

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
wclean baseCase/generateURot
wmake  baseCase/generateURot

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
		NzList=(32 64 128 256) # 256) #128 256 )
    elif [ "$meshType" = "tri" ];
    then 
        NzList=(22 44 88 174) # 256) #128 256 )
	else
		
		NzList=(0.05 0.021 0.0105 0.0051) # 0.01) #7.8e-3  3.9e-3 )
	fi

        mkdir --parents $series

        for n in ${!NzList[*]}
        do
            for m in ${!CoList[*]}
            do
                Co=${CoList[$m]}
                caseName=N${NzList[$n]}Co${Co}

                caseDir=$series/$caseName
                echo $caseDir
                cp -r baseCase $caseDir

                foamParmSet maxAlphaCo "$Co" $caseDir/system/controlDict
                foamParmSet application "$application" $caseDir/system/controlDict
                foamParmSet deltaT "${dtList[$n]}" $caseDir/system/controlDict

		mkdir $caseDir/logs
		sed -i "s/RECONSCHEME/$scheme/g" $caseDir/system/fvSolution

        if [ "$meshType" = "hex" ];
        then
            nx=${NzList[$n]}
            ny=${NzList[$n]}
            nz=${NzList[$n]}
            foamParmSet nx "$nx" $caseDir/system/blockMeshDict
            foamParmSet ny "$ny" $caseDir/system/blockMeshDict
            foamParmSet nz "$nz" $caseDir/system/blockMeshDict
		    blockMesh -case $caseDir > $caseDir/logs/blockMesh.log 2>&1
            cp -r $caseDir/0.org $caseDir/0
		    initAlphaField -case $caseDir > $caseDir/logs/initAlphaField.log 2>&1
		    generateURot -case $caseDir > $caseDir/logs/generateU.log 2>&1
		elif [ "$meshType" = "poly" ]
		then
			nz=${NzList[$n]}
			sed  -i "s/replaceMaxSize/$nz/" $caseDir/system/meshDict
			cp -r $caseDir/0.org $caseDir/0
			cd $caseDir
				
				pMesh  > $caseDir/logs/pMesh.log 2>&1
				cd $curDir
			            changeDictionary > $caseDir/logs/changeDictionary.log 2>&1
				    initAlphaField -case $caseDir > $caseDir/logs/initAlphaField.log 2>&1
	               		    generateURot -case $caseDir > $caseDir/logs/generateU.log 2>&1

        else
                    #cp $tetMeshDir/polyMesh_${NzList[$n]}/* $caseDir/constant/polyMesh/
		    #foamParmSet 'internalField' "uniform (1 0 0)" $caseDir/0.org/U
                    nx=${NzList[$n]}
                    ny=${NzList[$n]}
                    nz=${NzList[$n]}
		    sed  -i "s/replaceNx/$nx/" $caseDir/triBlock.geo
		    sed  -i "s/replaceNz/$nz/" $caseDir/triBlock.geo
		    nz=${NzList[$n]}
		    #sed  -i "s/replaceMaxSize/$nz/" $caseDir/system/meshDict
            cp -r $caseDir/0.org $caseDir/0
		    cd $caseDir
		    gmsh -3 -optimize_netgen triBlock.geo > $caseDir/logs/gmsh.log 2>&1
		    gmshToFoam triBlock.msh > $caseDir/logs/gmshToFoam.log 2>&1
		    #tetMesh  > $caseDir/logs/tetMesh.log 2>&1
		    changeDictionary > $caseDir/logs/changeDictionary.log 2>&1
		    cd $curDir
		    initAlphaField -case $caseDir > $caseDir/logs/initAlphaField.log 2>&1
		    generateURot -case $caseDir > $caseDir/logs/generateU.log 2>&1
       fi

            done
        done
    done
done

fi
