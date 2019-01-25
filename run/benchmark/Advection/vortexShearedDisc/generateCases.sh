#!/bin/bash

#Input

#appList=(interFlow mulesFoam passiveAdvectionFoam passiveAdvectionFoam)
appList=(advectorVoF) #interFlow advectorVoF 
#schemeList=(isoAdvector MULES HRIC CICSAM)
schemeList=(plicRDF) # isoAlpha) # RDF isoInverseDistance  RDFCell  RDFCellAdvect RDFCell
meshList=( tri hex poly ) #
#meshList=(poly)
CoList=(0.5)
#CoList=(0.1 0.2 0.5)
#End of input
#VoFLib_ROOT_DIR=PATH TO VoFLib 
curDir=$PWD

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

wclean baseCase/generateUVortex2D/
wmake baseCase/generateUVortex2D/


#Looping through all mesh, solver, and CFL number combinations and generating cases
for nn in ${!meshList[*]}
do

    meshType=${meshList[$nn]}
    for mm in ${!appList[*]}
    do

	#echo "appList"
        application=${appList[$mm]}
        scheme=${schemeList[$mm]}

        #Case location
        series=$PWD/$scheme/$meshType

	if [ "$meshType" = "hex" ];
	then
		NzList=(64 128 256 512 1024) # 256) #128 256 )
    elif [ "$meshType" = "tri" ];
    then 
        NzList=(40 80 158 317 634) # 256) #128 256 )
	else
		
		NzList=(56 112 224 448 896) # 0.01) #7.8e-3  3.9e-3 )
	fi

        mkdir --parents $series

	#echo $NzList

        for n in ${!NzList[*]}
        do

	    #echo "NzList "
            for m in ${!CoList[*]}
            do

	        #echo "CoList "
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

		sed -i "s/RECONSCHEME/$scheme/g" $caseDir/system/fvSolution
                foamParmSet maxAlphaCo "$Co" $caseDir/system/controlDict
                foamParmSet application "$application" $caseDir/system/controlDict

		
		mkdir $caseDir/logs

                if [ "$meshType" = "hex" ];
                then
                    nx=${NzList[$n]}
                    nz=${NzList[$n]}
                    foamParmSet nx "$nx" $caseDir/system/blockMeshDict
                    foamParmSet nz "$nz" $caseDir/system/blockMeshDict
                    cp -r $caseDir/0.org $caseDir/0
		            blockMesh -case $caseDir > $caseDir/logs/blockMesh.log 2>&1
                else
                    #cp $triMeshDir/polyMesh_${NzList[$n]}/* $caseDir/constant/polyMesh/
		            nx=${NzList[$n]}
		            nz=${NzList[$n]}
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
                        
                        #Convert from tet to poly mesh
                        polyDualMesh -case $caseDir -overwrite 160 > $caseDir/logs/polyDualMesh.log 
                        #Remove backmost  part of cells
                        topoSet -case $caseDir > $caseDir/logs/topoSet.log
                        subsetMesh -case $caseDir -overwrite c0 -patch front > $caseDir/logs/subsetMesh.log
                        transformPoints -case $caseDir -translate '(0 -0.5 0)' > $caseDir/logs/translate.log 
                        transformPoints -case $caseDir -scale '(1 2 1)' > $caseDir/logs/scale.log 
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

fi
