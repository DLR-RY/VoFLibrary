# coding: utf-8
import pandas as pd
import numpy as np
hexplicRDF = pd.read_csv("results/hex-plicRDF",delim_whitespace=True,header=None)
hexplicRDFlatex = pd.DataFrame()
hexplicRDFlatex['Ev'] = hexplicRDF[3].map('{:,.3e}'.format)
hexplicRDFlatex['Ebound'] = hexplicRDF[4].map('{:,.3e}'.format)
hexplicRDFlatex['Emass'] = hexplicRDF[2].map('{:,.2e}'.format)
hexplicRDFlatex['O'] = np.log(hexplicRDF[2].shift(1)/hexplicRDF[2])/np.log(2)
hexplicRDFlatex['O'] = hexplicRDFlatex['O'].map('{:,.2f}'.format)
#hexplicRDFlatex.iloc[1:-1,3] = hexplicRDFlatex.iloc[:-1,2]/hexplicRDFlatex.iloc[1:,2]
hexplicRDFlatex['Te'] = hexplicRDF[5]+hexplicRDF[6]
hexplicRDFlatex['Tr'] = hexplicRDF[5]
hexplicRDFlatex.to_latex('results/hexData.tex')

triplicRDF = pd.read_csv("results/tri-plicRDF",delim_whitespace=True,header=None)
triplicRDFlatex = pd.DataFrame()
triplicRDFlatex['Ev'] = triplicRDF[3].map('{:,.3e}'.format)
triplicRDFlatex['Ebound'] = triplicRDF[4].map('{:,.3e}'.format)
triplicRDFlatex['Emass'] = triplicRDF[2].map('{:,.2e}'.format)
triplicRDFlatex['O'] = np.log(triplicRDF[2].shift(1)/triplicRDF[2])/np.log(2)
triplicRDFlatex['O'] = triplicRDFlatex['O'].map('{:,.2f}'.format)
#triplicRDFlatex.iloc[1:-1,3] = triplicRDFlatex.iloc[:-1,2]/triplicRDFlatex.iloc[1:,2]
triplicRDFlatex['Te'] = triplicRDF[5]+triplicRDF[6]
triplicRDFlatex['Tr'] = triplicRDF[5]
triplicRDFlatex.to_latex('results/triData.tex')

polyplicRDF = pd.read_csv("results/poly-plicRDF",delim_whitespace=True,header=None)
polyplicRDFlatex = pd.DataFrame()
polyplicRDFlatex['Ev'] = polyplicRDF[3].map('{:,.3e}'.format)
polyplicRDFlatex['Ebound'] = polyplicRDF[4].map('{:,.3e}'.format)
polyplicRDFlatex['Emass'] = polyplicRDF[2].map('{:,.2e}'.format)
polyplicRDFlatex['O'] = np.log(polyplicRDF[2].shift(1)/polyplicRDF[2])/np.log(2)
polyplicRDFlatex['O'] = polyplicRDFlatex['O'].map('{:,.2f}'.format)
#polyplicRDFlatex.iloc[1:-1,3] = polyplicRDFlatex.iloc[:-1,2]/polyplicRDFlatex.iloc[1:,2]
polyplicRDFlatex['Te'] = polyplicRDF[5]+polyplicRDF[6]
polyplicRDFlatex['Tr'] = polyplicRDF[5]
polyplicRDFlatex.to_latex('results/polyData.tex')
