import re, analysisSeq
import matplotlib.pyplot as plt

fa_in = open("P-domain_pab1-all.fa","r")
#fa_in = open("Conserved-domain_pab1-all.fa","r")

# create the dictionary where the key is the name of the species and the value is the sequence
fa_dic = {}
name = ''
seq = ['W','F','Y','I','L','V','M','C','A','H','T','G','P','N','Q','S','D','E','K','R']

for line in fa_in.readlines():
  line = line.rstrip()
  if line[0] == ">":
    # get the name of the species
    line_lstrip = line.lstrip('>')
    line_split = re.split(r'/',line_lstrip)
    name_temp = line_split[0]
    name = name_temp
    fa_dic[name_temp] = ''
  else:
    fa_dic[name] = fa_dic[name] + line

fa_hydro = {}
# calculate the propotion of amino acids in a sequence according to the hydrophobicity
for i in fa_dic:
  fa_hydro[i] = []
  seq_gapfree = fa_dic[i].replace('-','')
  length = len(seq_gapfree)
  fa_hydro[i].append(seq_gapfree.count('W')/length)
  fa_hydro[i].append(seq_gapfree.count('F')/length)  
  fa_hydro[i].append(seq_gapfree.count('Y')/length)
  fa_hydro[i].append(seq_gapfree.count('I')/length)
  fa_hydro[i].append(seq_gapfree.count('L')/length)
  fa_hydro[i].append(seq_gapfree.count('V')/length)
  fa_hydro[i].append(seq_gapfree.count('M')/length)
  fa_hydro[i].append(seq_gapfree.count('C')/length)
  fa_hydro[i].append(seq_gapfree.count('A')/length)
  fa_hydro[i].append(seq_gapfree.count('H')/length)
  fa_hydro[i].append(seq_gapfree.count('T')/length)
  fa_hydro[i].append(seq_gapfree.count('G')/length)
  fa_hydro[i].append(seq_gapfree.count('P')/length)
  fa_hydro[i].append(seq_gapfree.count('N')/length)
  fa_hydro[i].append(seq_gapfree.count('Q')/length)
  fa_hydro[i].append(seq_gapfree.count('S')/length)
  fa_hydro[i].append(seq_gapfree.count('D')/length)
  fa_hydro[i].append(seq_gapfree.count('E')/length)
  fa_hydro[i].append(seq_gapfree.count('K')/length)
  fa_hydro[i].append(seq_gapfree.count('R')/length)

fa_hydro_sum = {}
# calculate the sum propotion of amino acids of which hydrophobicity is more than 0 according to Hopp.Woods
for i in fa_dic:
  fa_hydro_sum[i] = []
  temp = 0
  for j in range(0,11):
    temp += fa_hydro[i][j]
  fa_hydro_sum[i].append(temp)
  fa_hydro_sum[i].append(1-temp)

fa_hydro_sumP = {}
# calculate the sum propotion of amino acids of which hydrophobicity is no less than 0 according to Hopp.Woods
for i in fa_dic:
  fa_hydro_sumP[i] = []
  temp = 0
  for j in range(0,13):
    temp += fa_hydro[i][j]
  fa_hydro_sumP[i].append(temp)
  fa_hydro_sumP[i].append(1-temp)

fa_hydro_sumClass = {}
# calculate by the classficaiton of Alkyl(Hydrophobic), Aromatic(Hydrophobic), Neutral(Hydrophilic), Acidic(Hydrophilic), Basic(Hydrophilic)
for i in fa_dic:
  fa_hydro_sumClass[i] = []
  Al = [3,4,5,6,8,11,12]
  Ar = [0,1]
  Ne = [2,7,10,13,14,15]
  Ac = [16,17]
  Ba = [9,18,19]
  temp = [0,0,0,0,0]
  for j in Al:
    temp[0] += fa_hydro[i][j]
  fa_hydro_sumClass[i].append(temp[0])
  for j in Ar:
    temp[1] += fa_hydro[i][j]
  fa_hydro_sumClass[i].append(temp[1])
  for j in Ne:
    temp[2] += fa_hydro[i][j]
  fa_hydro_sumClass[i].append(temp[2])
  for j in Ac:
    temp[3] += fa_hydro[i][j]
  fa_hydro_sumClass[i].append(temp[3])
  for j in Ba:
    temp[4] += fa_hydro[i][j]
  fa_hydro_sumClass[i].append(temp[4])


#analysisSeq.draw(fa_dic, 10)
#analysisSeq.drawMax(fa_hydro, 4)
#analysisSeq.drawMin(fa_hydro, 4)


'''overall = analysisSeq.calculate(fa_dic)
print(overall)
'''
prop = analysisSeq.calculatePro(fa_dic)
print(prop)


#analysisSeq.draw2CP(fa_hydro_sumP)

#analysisSeq.draw2C(fa_hydro_sum)

#analysisSeq.draw5C(fa_hydro_sumClass)

'''{'A': 0.1425265359709829, 'P': 0.16448964151405238, 'Y': 0.029751440842030307, 'N': 0.05349473731057656, 'F': 0.030641536304768686, 'H': 0.009546273837869112, 'I': 0.021117514853468033, 'L': 0.027348183092636685, 'C': 0.00022252386568459467, 'E': 0.005985891986915597, 'M': 0.0584125147422061, 'S': 0.03633814726629431, 'T': 0.037005718863348094, 'R': 0.05876855292730145, 'Q': 0.1361623534124035, 'W': 0.003471372304679677, 'G': 0.13324729077193528, 'K': 0.006430939718284786, 'D': 0.003004072186742028, 'V': 0.042034758227819935}'''

'''{'A': 0.07956851124702585, 'P': 0.022800589720758463, 'Y': 0.03471178127782562, 'N': 0.06254835272307939, 'F': 0.06159954457208752, 'H': 0.019720612492154087, 'I': 0.049819726451311545, 'L': 0.06736537872042273, 'C': 0.011809012217729575, 'E': 0.07761250675113492, 'M': 0.0261725079188988, 'S': 0.06810982819273943, 'T': 0.04095931802589516, 'R': 0.05263695680733356, 'Q': 0.023603427386982352, 'W': 0.004174755864364226, 'G': 0.07734975987855255, 'K': 0.08305720583298057, 'D': 0.06749675215671391, 'V': 0.06888347176200972}'''


  
  




    


  
  
