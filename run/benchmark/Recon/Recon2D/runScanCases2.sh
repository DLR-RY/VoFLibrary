#!/bin/bash

reconSchemeList=(plicRDFIter2 plicRDFIter3 plicRDFIter5 isoRDFIter2 isoRDFIter3 isoRDFIter5 isoAlpha)
meshList=(hex tri poly)


curDir=$PWD

for nn in ${!reconSchemeList[*]}
do
	reconScheme=${reconSchemeList[$nn]}

	for mm in ${!meshList[*]}
	do
		
		mesh=${meshList[$mm]}
		meshType=${meshList[$mm]}
		casesDir=$curDir/$reconScheme/$mesh

		if [ "$meshType" = "hex" ];
		then
			#Domain  dimensions
			L=2
			H=2
			#Vertical velocity component
			Uz=0.5
			NzList=(  15 30 60 120 240 480)
			#NxList=(25 50 100 200 400 800)
		else
			NzList=( 10 15 30 60 120 240)
			#NxList=(75 150 300 600 1200)
			Uz=0.0
		fi



		for n in ${!NzList[*]} 
		do
			#for m in ${!CoList[*]}
			#do
				#Co=${CoList[$m]}
				caseName=N${NzList[$n]}
				caseDir=$casesDir/$caseName
				echo $caseDir
				cd $caseDir
				./Allrun
				cd $casesDir

			#done
		done
	

	done
done
