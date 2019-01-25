#!/bin/bash

#schemeList=(isoAdvector HRIC CICSAM MULES)
#schemeList=(isoAdvector isoAdvectorRDF isoAdvectorPerfectRDF)
schemeList=( plicRDF) # isoInverseDistance RDFadvect
meshList=(hex tri poly) # hex tri poly
dtList=(0.01 0.005 0.0025 )
#dtList=( 0.0025 )
CoList=(  0.5)
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
		NxList=( 32 64 128 ) #256 ) #128 256 )
		NzList=( 64 128 256 ) #512 ) #128 256 )
    elif [ "$meshType" = "tri" ];
    then 
        NxList=(21 42 87 ) #174) # 256) #128 256 )
        NzList=(42 84 174 ) # 348) # 256) #128 256 )
	else
		
		NzList=(0.047 0.021 0.0103 ) # 0.0051) # 0.01) #7.8e-3  3.9e-3 )
	fi


        mkdir --parents $series

        for n in ${!NzList[*]}
        do
            for m in ${!CoList[*]}
            do
                Co=${CoList[$m]}

                caseName=N${NzList[$n]}Co${Co}

                caseDir=$casesDir/$caseName
                echo $caseDir
		cd $caseDir
		./Allrun 
#		sleep 1m
                
            done
        done
    done
done
