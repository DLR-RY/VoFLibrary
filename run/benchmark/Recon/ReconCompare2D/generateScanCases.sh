#!/bin/bash

#appList=(reconstructError)
#reconSchemeList=(isoInverseDistance isoAdvection RDFCell RDFPoints gradAlphaSmoothed gradAlpha RDF perfectRDFPoints)
reconSchemeList=(plicRDF isoRDF isoAlpha)
#appList=(isoAdvector isoAdvectorRDF isoAdvectorPerfectRDF isoAdvectorRDFITER)
schemeList=(reconstructError reconstructError reconstructError)
#schemeList=(isoAdvector isoAdvectorRDF isoAdvectorPerfectRDF isoAdvectorRDFITER)
meshList=( tri ) # tri poly)
#CoList=(0)

#RDF


application=reconstructError
curDir=$PWD

#Location of tri meshes
triMeshDir=../triMeshes

for nn in ${!meshList[*]}
do
	meshType=${meshList[$nn]}
	for mm in ${!reconSchemeList[*]}
	do

		scheme=${schemeList[$mm]}
		reconScheme=${reconSchemeList[$mm]}

		#Case location
		series=$PWD/$reconScheme/$meshType

		if [ "$meshType" = "hex" ];
		then
			#Domain  dimensions
			L=2
			H=2
			#Vertical velocity component
			Uz=0.5
			NzList=(15 30 60 120 240 480)
			#NxList=(25 50 100 200 400 800)
		else
			NzList=(15 30 60 120 240)
			#NxList=(75 150 300 600 1200)
			Uz=0.0
		fi

		mkdir --parents $series

		for n in ${!NzList[*]} 
		do
			#for m in ${!CoList[*]}
			#do
				#Co=${CoList[$m]}
				caseName=N${NzList[$n]} #Co${Co}
				caseDir=$series/$caseName
				echo $caseDir
				cp -r baseCase $caseDir
				
				#./ofset maxAlphaCo "$Co" $caseDir/system/controlDict
				./ofset application "$application" $caseDir/system/controlDict
				sed -i "s/RECONSCHEME/$reconScheme/g" $caseDir/system/fvSolution
				#./ofset Uz "$Uz" $caseDir/0.org/U

				#Generating mesh
				mkdir $caseDir/logs
				cp -r $caseDir/0.org $caseDir/0
				touch $caseDir/case.foam
				if [ "$meshType" = "hex" ];
				then
					nx=${NzList[$n]}
					nz=${NzList[$n]}
					#nz=${NzList[$n]}
					./ofset L "$L" $caseDir/constant/polyMesh/blockMeshDict
					./ofset H "$H" $caseDir/constant/polyMesh/blockMeshDict
					./ofset nx "$nx" $caseDir/constant/polyMesh/blockMeshDict
					./ofset nz "$nz" $caseDir/constant/polyMesh/blockMeshDict
					blockMesh -case $caseDir > $caseDir/logs/blockMesh.log 2>&1
					setAlphaField -case $caseDir > $caseDir/logs/perfectVOF.log 2>&1
				else
					#cp $triMeshDir/N${NzList[$n]}/* $caseDir/constant/polyMesh/
					nx=${NzList[$n]}
					nz=${NzList[$n]}
					#nz=${NzList[$n]}
					./ofset2 replaceNx "$nx" $caseDir/triSquare.geo
					./ofset2 replaceNz "$nz" $caseDir/triSquare.geo
					cd $caseDir
					gmsh -3 triSquare.geo -optimize_netgen > $caseDir/logs/gmsh.log 2>&1
					gmshToFoam triSquare.msh > $caseDir/logs/gmshToFoam.log 2>&1
					changeDictionary > $caseDir/logs/changeDictionary.log 2>&1
					if [ "$meshType" = "poly" ];
					then
						#Convert from tet to poly mesh
						polyDualMesh -case $caseDir -overwrite 160 > $caseDir/logs/polyDualMesh.log 2>&1
						#Remove backmost  part of cells
						topoSet -case $caseDir > $caseDir/logs/topoSet.log 2>&1
						subsetMesh -case $caseDir -overwrite c0 -patch front > $caseDir/logs/subsetMesh.log 2>&1
						transformPoints -case $caseDir -translate '(0 -0.025 0 )' > /dev/null
						rm -rf *.obj
						setAlphaField > logs/perfectVOF.log 2>&1
					else
						setAlphaField > logs/perfectVOF.log 2>&1			
					fi

					cd $curDir
				fi
				
			#done
		done
	done
done


