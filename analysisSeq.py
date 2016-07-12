import matplotlib.pyplot as plt
import numpy as np

# scales of hydrophobicity according to the Hopps.Wood scales
seq = ['W','F','Y','I','L','V','M','C','A','H','T','G','P','N','Q','S','D','E','K','R']
hydroSeq = {'W':-3.4,'F':-2.5,'Y':2.3,'I':-1.8,'L':-1.8,'V':-1.5,'M':-1.3,'C':-1,'A':-0.5,'H':-0.5,'T':-0.4,'G':0,'P':0,'N':0.2,'Q':0.2,'S':0.2,'D':3,'E':3,'K':3,'R':3,'-':0}

# calculate all the amount of amino acids
def calculate(fa_dic):
  d = {'W':0,'F':0,'Y':0,'I':0,'L':0,'V':0,'M':0,'C':0,'A':0,'H':0,'T':0,'G':0,'P':0,'N':0,'Q':0,'S':0,'D':0,'E':0,'K':0,'R':0,'-':0}
  for i in fa_dic:
    for j in range(0,len(fa_dic[i])):
      d[fa_dic[i][j]] += 1
  del d['-']
  return d

# calculate overall propotion of amino acids in multiple sequences
def calculatePro(fa_dic):
  d = calculate(fa_dic)
  prop = {}
  length = 0
  for i in d:
    length += d[i]
  for i in d:
    prop[i] = d[i]/length
  return prop
  
# draw the line chart of different classes of amino acids of mulitiple sequences where first class scale is no more than 0
def draw2CP(hydro):
  dict = {'hydrophobic':[],'hydrophilic':[],'name':[]}
  t = len(hydro)
  for i in hydro:
    dict['hydrophobic'].append(hydro[i][0])
    dict['hydrophilic'].append(hydro[i][1])
    dict['name'].append(i)
  plot1 = plt.plot(dict['hydrophobic'],'b',label='scale<=0')
  plot2 = plt.plot(dict['hydrophilic'],'g',label='scale>0')
  plt.title('Propotion of different classes of amino acids')
  plt.xlabel('Species')
  plt.ylabel('Propotion')
  plt.xlim(0,t)
  plt.ylim(0,1)
  plt.legend()
  plt.tight_layout()
  plt.show()
  return 

# draw the line chart of different classes of amino acids of mulitiple sequences where first class scale is less than 0
def draw2C(hydro):
  dict = {'hydrophobic':[],'hydrophilic':[],'name':[]}
  t = len(hydro)
  for i in hydro:
    dict['hydrophobic'].append(hydro[i][0])
    dict['hydrophilic'].append(hydro[i][1])
    dict['name'].append(i)
  plot1 = plt.plot(dict['hydrophobic'],'b',label='scale<0')
  plot2 = plt.plot(dict['hydrophilic'],'g',label='scale>=0')
  plt.title('Propotion of different classes of amino acids')
  plt.xlabel('Species')
  plt.ylabel('Propotion')
  plt.xlim(0,t)
  plt.ylim(0,1)
  plt.legend()
  plt.tight_layout()
  plt.show()
  return 

# draw the line chart of different classes of amino acids of mulitiple sequences according to the normal standard of classification
def draw5C(hydro):
  dict = {'Alkyl':[],'Aromatic':[],'Neutral':[],'Acidic':[],'Basic':[],'name':[]}
  t = len(hydro)
  for i in hydro:
    dict['Alkyl'].append(hydro[i][0])
    dict['Aromatic'].append(hydro[i][1])
    dict['Neutral'].append(hydro[i][2])
    dict['Acidic'].append(hydro[i][3])
    dict['Basic'].append(hydro[i][4])
    dict['name'].append(i)
  plot1 = plt.plot(dict['Alkyl'],'b',label='Alkyl')
  plot2 = plt.plot(dict['Aromatic'],'g',label='Aromatic')
  plot3 = plt.plot(dict['Neutral'],'k',label='Neutral')
  plot4 = plt.plot(dict['Acidic'],'c',label='Acidic')
  plot5 = plt.plot(dict['Basic'],'y',label='Basic')
  plt.title('Propotion of different classes of amino acids')
  plt.xlabel('Species')
  plt.ylabel('Propotion')
  plt.xlim(0,t)
  plt.ylim(0,0.9)
  plt.legend(fontsize=8.5)
  plt.tight_layout()
  plt.show()
  return 


# return the n maximum amino acids 
def maxN(fa_hydro, n):
  if(n > 20):
    print('Warining: n is beyond index, please ensure n is no more than 20!')
    return 
  else:
    fa_hydro_max = {}
    for i in fa_hydro:
      fa_hydro_max[i] = []
      temp = list(fa_hydro[i])
      standard = list(fa_hydro[i])
      for j in range(0,n):
        fa_hydro_max[i].append(standard.index(max(temp)))
        standard[standard.index(max(temp))] = -1
        del temp[temp.index(max(temp))]
    return fa_hydro_max

# draw the bar graph of maximum amino acids in multiple sequences
def drawMax(fa_hydro, n):
  result = maxN(fa_hydro, n)
  seq = []
  amount = []
  for i in result:
    for j in range(0,n):
      seq.append(result[i][j])
  for k in range(0,20):
    amount.append(seq.count(k))
  # draw the bar graph
  n_group = 20
  fig, ax = plt.subplots()
  index = np.arange(n_group)
  bar_width = 0.7
  opacity = 0.7
  plt.bar(index,amount,bar_width,alpha=opacity,color='g')
  plt.hist(amount,20,normed=1, facecolor='b',alpha=0.75)
  plt.xlabel('Amino acids')
  plt.ylabel('Frequency of amino acids')
  plt.xticks(index + bar_width/2,('W','F','Y','I','L','V','M','C','A','H','T','G','P','N','Q','S','D','E','K','R'))
  plt.xlim(0,20)
  plt.title('Most popular amino acids')
  plt.tight_layout()
  plt.show()
 

# return the n minimum amino acids 
def minN(fa_hydro, n):
  if(n > 20):
    print('Warning: n is beyond index, please ensure n is no more than 20!')
    return 
  else:
    fa_hydro_min = {}
    for i in fa_hydro:
      fa_hydro_min[i] = []
      temp = list(fa_hydro[i])
      standard = list(fa_hydro[i])
      for j in range(0,n):
        fa_hydro_min[i].append(standard.index(min(temp)))
        standard[standard.index(min(temp))] = 1.7
        del temp[temp.index(min(temp))]
    return fa_hydro_min

# draw the bar graph of minimum amino acids in multiple sequences
def drawMin(fa_hydro, n):
  result = minN(fa_hydro, n)
  seq = []
  amount = []
  for i in result:
    for j in range(0,n):
      seq.append(result[i][j])
  for k in range(0,20):
    amount.append(seq.count(k))
  n_group = 20
  fig, ax = plt.subplots()
  index = np.arange(n_group)
  bar_width = 0.7
  opacity = 0.7
  plt.bar(index,amount,bar_width,alpha=opacity,color='g')
  plt.hist(amount,20,normed=1, facecolor='b',alpha=0.75)
  plt.xlabel('Amino acids')
  plt.ylabel('Frequency of amino acids')
  plt.xticks(index + bar_width/2,('W','F','Y','I','L','V','M','C','A','H','T','G','P','N','Q','S','D','E','K','R'))
  plt.xlim(0,20)
  plt.title('Most unpopular amino acids')
  plt.tight_layout()
  plt.show()

# repetitively calculate the average scores of amino acids according to the Hopps.Wood scale with gaps
def averageScores(fa_dic, n):
  lengthSeq = []
  fa_gapfree = {}
  for i in fa_dic:
    seq_gapfree = fa_dic[i].replace('-','')
    fa_gapfree[i] = seq_gapfree
    lengthSeq.append(len(seq_gapfree))
  if (n > max(lengthSeq)):
    print('Warning! n is out of range of sequences!')
    return 
  else:
    fa_average = {}
    for each in fa_gapfree:
      length = len(fa_gapfree[each])
      if (n <= length):
        fa_average[each] = []
        for j in range(0,length-n+1):
          temp = 0
          for k in range(j,j+n):
            temp += hydroSeq[fa_gapfree[each][k]]
          fa_average[each].append(temp)
    return fa_average

# repetitively calculate the average scores of amino acids according to the Hopps.Wood scale without gaps
def averageScoresGap(fa_dic, n):
  lengthSeq = []
  fa_gap = {}
  for i in fa_dic:
    seq_gap = fa_dic[i]
    fa_gap[i] = seq_gap
    lengthSeq.append(len(seq_gap))
  if (n > max(lengthSeq)):
    print('Warning! n is out of range of sequences!')
    return 
  else:
    fa_average = {}
    for each in fa_gap:
      length = len(fa_gap[each])
      fa_average[each] = []
      for j in range(0,length-n+1):
        temp = 0
        for k in range(j,j+n):
          temp += hydroSeq[fa_gap[each][k]]
        fa_average[each].append(temp)
    return fa_average

# draw the superimposed image of multiple sequences with gaps
def draw(fa_dic,n):
  fa_average_gap = averageScoresGap(fa_dic,n)
  for i in fa_average_gap:
    plt.figure(num="SuperimposedImage"+ str(n),figsize=(8,6))
    plt.plot(fa_average_gap[i])
  plt.show()

# draw and output separate image of multiple sequences
def drawSep(fa_dic,n):
  fa_average = averageScores(fa_dic,n)
  for i in fa_average:
    plt.figure(num=i,figsize=(8,8))
    plt.plot(fa_average[i])
    string = i + ".png"
    plt.savefig(string,dpi=72)


  
  
            
      
        
      
    
        



                
              
      
    

