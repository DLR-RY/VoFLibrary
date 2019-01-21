#!/bin/bash

#schemeList=(isoAdvector HRIC CICSAM MULES)
#schemeList=(reconstructError reconstructError)
#reconSchemeList=(isoInverseDistance RDFPoints gradAlphaSmoothed gradAlpha perfectRDFPoints RDFPoints RDFCell)
reconSchemeList=(plicRDF isoRDF isoAlpha)
#schemeList=(isoAdvector isoAdvectorRDF isoAdvectorPerfectRDF isoAdvectorRDFITER)
meshList=( hex ) #hex poly
#meshList=(poly)
#CoList=( 0.1 0.2 0.5)
#CoList=(0.1 0.2 0.5)
#meshList=(hex)

curDir=$PWD

for nn in ${!reconSchemeList[*]}
do
	reconScheme=${reconSchemeList[$nn]}

	for mm in ${!meshList[*]}
	do
		
		mesh=${meshList[$mm]}
		meshType=${meshList[$mm]}
		casesDir=$curDir/$reconScheme/$mesh

		#NzList=(5 10 15 30 60 120 240)
		NzList=( 10 15 30 60 120 240 )
		#NxList=(75 150 300 600 1200)
		Uz=0.0

		#NzList=(480)

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
