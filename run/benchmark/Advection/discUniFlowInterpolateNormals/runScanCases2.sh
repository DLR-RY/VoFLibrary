#!/bin/bash

#schemeList=(isoAdvector HRIC CICSAM MULES)
#schemeList=(isoAdvector isoAdvectorRDF isoAdvectorPerfectRDF)
schemeList=(plicRDFNoInterpol) #(plicRDFN isoAlpha)
#meshList=( poly)
meshList=(hex tri poly) # hex tri poly
CoList=(0.1 0.2 0.5)

curDir=$PWD

#Looping through all mesh, solver, and CFL number combinations and generated cases

    
for mm in ${!schemeList[*]}
do
    for nn in ${!meshList[*]}
    do
	meshType=${meshList[$nn]}
        scheme=${schemeList[$mm]}

        #Case location
        casesDir=$curDir/$scheme/$meshType

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
                caseDir=$casesDir/$caseName
                echo $caseDir
		cd $caseDir
		./Allrun 
#		sleep 1m
                
            done
        done
    done
done
