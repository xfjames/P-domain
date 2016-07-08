import re

fa_in = open("P-domain_pab1-all.fa","r")

# create the dictionary where the key is the name of the species and the value is the sequence
fa_dic = {}
name = ''

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
  fa_hydro[i].append(seq_gapfree.count('I')/length)
  fa_hydro[i].append(seq_gapfree.count('V')/length)  
  fa_hydro[i].append(seq_gapfree.count('L')/length)
  fa_hydro[i].append(seq_gapfree.count('F')/length)
  fa_hydro[i].append(seq_gapfree.count('C')/length)
  fa_hydro[i].append(seq_gapfree.count('M')/length)
  fa_hydro[i].append(seq_gapfree.count('A')/length)
  fa_hydro[i].append(seq_gapfree.count('G')/length)
  fa_hydro[i].append(seq_gapfree.count('T')/length)
  fa_hydro[i].append(seq_gapfree.count('S')/length)
  fa_hydro[i].append(seq_gapfree.count('W')/length)
  fa_hydro[i].append(seq_gapfree.count('Y')/length)
  fa_hydro[i].append(seq_gapfree.count('P')/length)
  fa_hydro[i].append(seq_gapfree.count('H')/length)
  fa_hydro[i].append(seq_gapfree.count('E')/length)
  fa_hydro[i].append(seq_gapfree.count('Q')/length)
  fa_hydro[i].append(seq_gapfree.count('D')/length)
  fa_hydro[i].append(seq_gapfree.count('N')/length)
  fa_hydro[i].append(seq_gapfree.count('K')/length)
  fa_hydro[i].append(seq_gapfree.count('R')/length)



  
  
