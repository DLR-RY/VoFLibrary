import pandas as pd

Hex = pd.DataFrame(columns=("isoRDFIDW","isoRDFLS","Time isoRDFIDW","Time isoRDFLS"))
Tri = pd.DataFrame(columns=("isoRDFIDW","isoRDFLS","Time isoRDFIDW","Time isoRDFLS"))
Poly = pd.DataFrame(columns=("isoRDFIDW","isoRDFLS","Time isoRDFIDW","Time isoRDFLS"))

df1 = pd.read_csv("hex-isoRDF2",delim_whitespace=True,header=None)
df2 = pd.read_csv("hex-isoRDF",delim_whitespace=True,header=None)


dftri1 = pd.read_csv("tri-isoRDF2",delim_whitespace=True,header=None)
dftri2 = pd.read_csv("tri-isoRDF",delim_whitespace=True,header=None)


dfpoly1 = pd.read_csv("poly-isoRDF2",delim_whitespace=True,header=None)
dfpoly2 = pd.read_csv("poly-isoRDF",delim_whitespace=True,header=None)


Hex['res'] =  0.25/df1[8]
Hex['isoRDFIDW'] =  df1[2].map('{:,.3e}'.format)
Hex['isoRDFLS'] =  df2[2].map('{:,.3e}'.format)
Hex['Time isoRDFIDW'] =  df1[6]
Hex['Time isoRDFLS'] =  df2[6]

Tri['res'] =  0.25/dftri1[8]
Tri['isoRDFIDW'] =  dftri1[2].map('{:,.3e}'.format)
Tri['isoRDFLS'] =  dftri2[2].map('{:,.3e}'.format)
Tri['Time isoRDFIDW'] =  dftri1[6]
Tri['Time isoRDFLS'] =  dftri2[6]

Poly['res'] =  0.25/dfpoly1[8]
Poly['isoRDFIDW'] =  dfpoly1[2].map('{:,.3e}'.format)
Poly['isoRDFLS'] =  dfpoly2[2].map('{:,.3e}'.format)
Poly['Time isoRDFIDW'] =  dfpoly1[6]
Poly['Time isoRDFLS'] =  dfpoly2[6]


Hex.set_index("res",inplace=True)
Tri.set_index("res",inplace=True)
Poly.set_index("res",inplace=True)




result = pd.concat([Hex,Tri,Poly], keys=['hex', 'tri', 'poly'])



result.to_latex("latexTabisoRDF.tex",escape=False)
