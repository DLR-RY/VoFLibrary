#!/bin/bash

#schemeList=(isoAdvector HRIC CICSAM MULES)
#schemeList=(isoAdvector isoAdvectorRDF isoAdvectorPerfectRDF)
schemeList=(plicRDF ) #isoAlpha
#meshList=( poly)
meshList=( hex poly tri ) #hex poly tri
CoList=(0.5)

curDir=$PWD

#Looping through all mesh, solver, and CFL number combinations and generating cases

    
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
		NzList=(64 128 256 512 1024) # 256) #128 256 )
    elif [ "$meshType" = "tri" ];
    then 
        NzList=(40 80 158 317 634) # 256) #128 256 )
	else
		
		NzList=(56 112 224 448 896) # 0.01) #7.8e-3  3.9e-3 )
	fi

        #mkdir --parents $series

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

                
            done
        done
    done
done
