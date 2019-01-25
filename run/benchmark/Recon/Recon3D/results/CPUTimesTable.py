import pandas as pd

cpuTimeHex = pd.DataFrame(columns=("isoAlpha","isoRDF","plicRDF","Young"))
cpuTimeTri = pd.DataFrame(columns=("isoAlpha","isoRDF","plicRDF","Young"))
cpuTimePoly = pd.DataFrame(columns=("isoAlpha","isoRDF","plicRDF","Young"))

df1 = pd.read_csv("hex-isoAlpha",delim_whitespace=True,header=None)
df2 = pd.read_csv("hex-isoRDF",delim_whitespace=True,header=None)
df3 = pd.read_csv("hex-plicRDF",delim_whitespace=True,header=None)
df4 = pd.read_csv("hex-gradAlpha",delim_whitespace=True,header=None)

dftri1 = pd.read_csv("tri-isoAlpha",delim_whitespace=True,header=None)
dftri2 = pd.read_csv("tri-isoRDF",delim_whitespace=True,header=None)
dftri3 = pd.read_csv("tri-plicRDF",delim_whitespace=True,header=None)
dftri4 = pd.read_csv("tri-gradAlpha",delim_whitespace=True,header=None)

dfpoly1 = pd.read_csv("poly-isoAlpha",delim_whitespace=True,header=None)
dfpoly2 = pd.read_csv("poly-isoRDF",delim_whitespace=True,header=None)
dfpoly3 = pd.read_csv("poly-plicRDF",delim_whitespace=True,header=None)
dfpoly4 = pd.read_csv("poly-gradAlpha",delim_whitespace=True,header=None)

cpuTimeHex['res'] =  0.25/df1[8]
cpuTimeHex['isoAlpha'] =  df1[6].map('{:,.3f}'.format)
cpuTimeHex['isoRDF'] =  df2[6].map('{:,.3f}'.format)
cpuTimeHex['plicRDF'] =  df3[6].map('{:,.3f}'.format)
cpuTimeHex['Young'] =  df4[6].map('{:,.3f}'.format)

cpuTimeTri['res'] =  0.25/dftri1[8]
cpuTimeTri['isoAlpha'] =  dftri1[6].map('{:,.3f}'.format)
cpuTimeTri['isoRDF'] =  dftri2[6].map('{:,.3f}'.format)
cpuTimeTri['plicRDF'] =  dftri3[6].map('{:,.3f}'.format)
cpuTimeTri['Young'] =  dftri4[6].map('{:,.3f}'.format)

cpuTimePoly['res'] =  0.25/dfpoly1[8]
cpuTimePoly['isoAlpha'] =  dfpoly1[6].map('{:,.3f}'.format)
cpuTimePoly['isoRDF'] =  dfpoly2[6].map('{:,.3f}'.format)
cpuTimePoly['plicRDF'] =  dfpoly3[6].map('{:,.3f}'.format)
cpuTimePoly['Young'] =  dfpoly4[6].map('{:,.3f}'.format)


cpuTimeHex.set_index("res",inplace=True)
cpuTimeTri.set_index("res",inplace=True)
cpuTimePoly.set_index("res",inplace=True)




result = pd.concat([cpuTimeHex,cpuTimeTri,cpuTimePoly], keys=['hex', 'tri', 'poly'])

#result['relisoRDF'] = result['isoRDF']/result['isoAlpha']
#result['relplicRDF'] = result['plicRDF']/result['isoAlpha']
#result['relisoRDFYoung'] = result['isoRDF']/result['Young']
#result['relplicRDFYoung'] = result['plicRDF']/result['Young']
#result['relisoAlphaYoung'] = result['isoAlpha']/result['Young']

result.to_latex("latexTab2.tex",escape=False)
