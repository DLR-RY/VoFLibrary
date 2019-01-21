#!/bin/bash

#appList=(reconstructError)
#reconSchemeList=(isoInverseDistance isoAdvection RDFCell RDFPoints gradAlphaSmoothed gradAlpha RDF perfectRDFPoints)
reconSchemeList=(plicRDFN isoRDFN)
#appList=(isoAdvector isoAdvectorRDF isoAdvectorPerfectRDF isoAdvectorRDFITER)

#schemeList=(isoAdvector isoAdvectorRDF isoAdvectorPerfectRDF isoAdvectorRDFITER)
IterCount=(2 3 5) # tri poly)
#CoList=(0)

#RDF



for mm in ${!reconSchemeList[*]}
do
	reconScheme=${reconSchemeList[$mm]}
	for nn in ${!IterCount[*]}
	do
	
		iterNum=${IterCount[$nn]}
		newRecon=$reconScheme$iterNum

		cp -r $reconScheme $newRecon
		find $newRecon -name "fvSolution" -exec sed -i -e "s/ITERCOUNT/$iterNum/g" {} \;
		

		
	done
	rm -r $reconScheme
done




