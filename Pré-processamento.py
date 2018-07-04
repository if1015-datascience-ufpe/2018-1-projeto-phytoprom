
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


#Criando data frame para matriz de fatores de transcrição
df = pd.read_csv(
    'PWM.csv',
    sep=',', header=None, names=['ID', 'Name', 'A', 'C', 'G', 'T'])
print(type(df), df.shape)
df


# In[3]:


print(df['A'][0])
print(df['C'][0])
print(df['G'][0])
print(df['T'][0])
type(df['A'][0])


# Primeiro, é necessário "limpar" os dados para representar da forma correta em python. Atualmente existe apenas dados em texto.

# In[4]:


#Função trata os dados do data frame transformando strings da matriz em valores inteiros
def f_str_list (string):
    string = string.split()
    listed = []
    i = 2
    while(string[i] != ']'):
            listed.append(int (string[i]))
            i +=1
    
    return listed


# In[5]:


np_a = df['A']
np_a.count()
for i in range (489):
    np_a[i] = f_str_list(np_a[i])
    
df['A'] = np_a
type(df['A'][0][0])


# In[6]:


np_c = df['C']
np_c.count()
for i in range (489):
    np_c[i] = f_str_list(np_c[i])
    
df['C'] = np_c
type(df['C'][0][0])


# In[7]:


np_g = df['G']
np_g.count()
for i in range (489):
    np_g[i] = f_str_list(np_g[i])
    
df['G'] = np_g
type(df['G'][0][0])


# In[8]:


np_t = df['T']
np_t.count()
for i in range (489):
    np_t[i] = f_str_list(np_t[i])
    
df['T'] = np_t
type(df['T'][0][0])


# In[9]:


#Atualizando data frame
df


# É mais interessante representar todas as colunas de nucleotídeos (A, C, G e T) como uma única coluna de numpayarray de 2 dimensões

# In[10]:


df['Matrix'] = df.apply(lambda x: np.array([x['A'], x['C'], x['G'], x['T']]).transpose(), axis=1)
df['Matrix'][0]


# Podemos ainda transformar a matriz com os valores absolutos em uma matriz de distribuição de probabilidades

# In[11]:


def probability(matrix):
    matrix = matrix.astype(float)
    
    for line in matrix:
        maximum = sum(line)
        for i in range(len(line)):
            line[i] = line[i]/maximum
    
    return matrix


# In[12]:


df['Matrix'] = df['Matrix'].apply(lambda x: probability(x))
df['Matrix'] = df['Matrix'].apply(lambda x: x.transpose())
df


# Em uma primeira abordagem, tentei extrair todos os fatores de transcrição que poderiam ser obtidos a partir da matriz, para então fazer a busca nos promotores. Esse método logo mostrou ser inadequado devido ao grande número de colunas de algumas matrizes. No pior caso, a matriz possui 30 colunas e variam para todas as 4 possibilidades, resultando em uma combinação de 4^30

# In[13]:


import queue

"""A função toma como argumento o conjunto de nucleotídeos (a,c,g,t) em forma de matriz(4xn),
extrai os fatores de transcrição e retorna em forma de lista"""
def combinatoria (a, c, g, t):
    tam = len(a)
    ft = [''] # Lista de fatores de transcrição
    cont = 1 # Tamanho da lista de fatores de transcrição obs:tamanho mínimo==1
    q = queue.Queue() #fila auxiliar
    
    # Verifica se possui valor para o nucleotídeo e adiciona na fila
    for i in range(tam):
        if a[i] != 0:
            q.put('A')
        if c[i] != 0:
            q.put('C')
        if g[i] != 0:
            q.put('G')
        if t[i]!=0:
            q.put('T')
        
        tam_fila = q.qsize() #Guarda o tamanho da fila
        indice=cont #Guarda o número de iterações que acontecerá na lista de fatores de transcrição
        cont = cont*tam_fila #Aumenta o tamanho da lista de fatores de transcrição (cont==tamanho da lista de FT)
        
        # A combinação sempre ocorre sobre todos os elementos da lista de fatores de transcrição
        if tam_fila>1: #Se o tamanho da fila for 1, não é necessário copiar os elementos para aumentar a lista de fatores de transcrição
            copy = ft #Faz uma cópia dos fatores de transcrição
            for j in range(tam_fila-1): #Copia os elementos para os novos índices da lista
                ft = ft + copy
        
        #Extrai todos os nucleotídeos da fila
        fila = []
        while not q.empty():
            for j in range(q.qsize()):
                fila.append(q.get())
                
        j=0 #variável auxiliar para o tamanho da lista de fatores de transcrição
        k=0 #variável auxiliar para interação sobre a fila
        z=0 #variável auxiliar para interação sobre a fila na lista de fatores de transcrição
        while (j < cont):
            ft[j] += fila[k]
            j = j+1
            z = z+1
            
            if z == indice:
                k = k+1
                z=0
        
    return ft


# In[14]:


lista = combinatoria(df['A'][0], df['C'][0], df['G'][0], df['T'][0])
print (lista)


# A solução encontrada para este problema foi representar toda a lista de fatores de transcrição em uma única string. Felizmente pesquisadores já haviam utilizado da mesma estratégia na representação de aminoácidos e hoje existe um padrão para esta representação definido pelo International Union of Pure and Applied Chemistry (IUPAC)

# In[15]:


"""
A função toma como argumento uma matriz de distribuição de probabilidades dos nucleotídeos (a,c,g,t) 4xn 
e retorna uma única string que representa todos os fatores de transcrição em forma de consenso pelo padrão IUPAC
"""
def degeneration (matrix):
    lst = []
    a = matrix[0,]
    c = matrix[1,]
    g = matrix[2,]
    t = matrix[3,]
    
    
    tam = len(a)
    i=0
    while i < tam:
        ia = 0
        ic = 0
        ig = 0
        it = 0
        if a[i] >= 0.3:
            ia = 1
        if c[i] >= 0.3:
            ic = 1
        if g[i] >= 0.3:
            ig = 1
        if t[i] >= 0.3:
            it = 1
        
        
        if ia + ic + ig + it == 4:
            lst.append ('N')
        elif ia + ic + ig + it == 1:
            if ia == 1:
                lst.append ('A')
            elif ic == 1:
                lst.append ('C')
            elif ig == 1:
                lst.append ('G')
            elif it == 1:
                lst.append ('T')
        
        elif ia + ic + ig + it == 3:
            if ia & ic & ig:
                lst.append ('V')
            elif ia & ic & it:
                lst.append ('H')
            elif ia & ig & it:
                lst.append ('D')
            elif ic & ig & it:
                lst.append ('B')
        
        elif ia + ic + ig + it == 2:
            if ia & ic:
                lst.append('M')
            elif ia & ig:
                lst.append('R')
            elif ia & it:
                lst.append('W')
            elif ic & ig:
                lst.append ('S')
            elif ic & it:
                lst.append ('Y')
            elif ig & it:
                lst.append('K')
        
        else:
            lst.append('N')
        
        
        i +=1
        
    return lst


# In[16]:


df['Motifs'] = df.apply(lambda x: degeneration(x['Matrix']), axis=1)


# In[17]:


df


# Com esses motivos, já é possível iniciar testes na busca de ECRs

# In[18]:


#Conjunto de teste para mineração de ECRs
promotor = pd.read_csv(
    'promotores teste.txt',
    sep=',', header=None, names=['Id','Promotor'])
print(type(promotor), promotor.shape)
promotor


# In[19]:


"""
A função toma como argumento uma região gênica, denominada promotor, o tamanho dessa região
e um motivo (ECR) que será procurado dentro dessa região. A função retorna uma lista com todos
os locais que o elemento é encontrado.
"""
def searchFT (promotor, tamanho, motivo):
    promotor = promotor.lower()
    degeneracao = {
        'N' : ['a', 'c', 'g', 't'],
        'V' : ['a', 'c', 'g'],
        'H' : ['a', 'c', 't'],
        'D' : ['a', 'g', 't'],
        'B' : ['c', 'g', 't'],
        'M' : ['a', 'c'],
        'R' : ['a', 'g'],
        'W' : ['a', 't'],
        'S' : ['c', 'g'],
        'Y' : ['c', 't'],
        'K' : ['g', 't'],
        'A' : ['a'],
        'C' : ['c'],
        'G' : ['g'],
        'T' : ['t'],
    }
    
    lst = []
    tam = len(motivo)
    i = 0
    
    while (i < tamanho-tam-1):
        aminoacido = motivo[0]
        if (promotor[i] in degeneracao[aminoacido]):
            match = True
            j = 1
            while ((j < tam) & (match == True)):
                aminoacido = motivo[j]
                if (promotor[i+j] not in degeneracao[aminoacido]):
                    j=0
                    match = False
                    break
                else:
                    j += 1
            
            if ((match==True) & (j == tam)):
                lst.append(i)
        i += 1
        
    return lst
        
        


# In[20]:


print (searchFT(promotor.Promotor[0], 5000, 'AAAG'))


# In[21]:


FT = searchFT(promotor.Promotor[0], 5000, df['Motifs'][0])


# In[22]:


print (len(FT), FT)


# Como estamos fazendo a busca em uma única fita de RNA, é necessário buscar o mesmo ECR no complemento reverso, pois é possível que este mesmo ECR esteja representado na fita negativa.

# In[23]:


"""
 A função calcula o complemento reverso de um motivo para poder fazer a pesquisa do ECR através de
 uma única fita de RNA.
"""
def reverseComplement(lst):
    tam = len(lst)
    i=0
    while (i < tam):
        if lst[i] == 'V':
            lst[i] = 'B'
        elif lst[i] == 'H':
            lst[i] = 'D'
        elif lst[i] == 'D':
            lst[i] = 'H'
        elif lst[i] == 'B':
            lst[i] = 'V'
        elif lst[i] == 'M':
            lst[i] = 'K'
        elif lst[i] == 'R':
            lst[i] = 'Y'
        elif lst[i] == 'Y':
            lst[i] = 'R'
        elif lst[i] == 'A':
            lst[i] = 'T'
        elif lst[i] == 'T':
            lst[i] = 'A'
        elif lst[i] == 'C':
            lst[i] = 'G'
        elif lst[i] == 'G':
            lst[i] = 'C'
        else:
            pass
        
        i += 1
    
    return lst


# In[24]:


#Cria uma coluna no data frame para o complemento reverso (a coluna só recebe o motivo reverso)
df['ReverseComplement'] = df['Motifs'].apply(lambda x: x[::-1])
#Calcula o complemento do motivo (para resultar no complemento reverso)
df.ReverseComplement.apply(lambda x: reverseComplement(x))
df


# In[25]:


#Minerando ECRs no conjunto de promotores teste
procura = df.Motifs.apply(lambda x: searchFT(promotor.Promotor[0], 5000, x))
procura2 = df.ReverseComplement.apply(lambda x: searchFT(promotor.Promotor[0], 5000, x))
print(procura, len(procura))
print (procura2, len(procura2))


# Salvando o data frame em planilha

# In[26]:


df.to_csv(path_or_buf='df.csv',
          sep=',', header=True, index=False, mode='w',
         columns=['ID', 'Name', 'A', 'C', 'G', 'T', 'Matrix', 'Motifs', 'ReverseComplement'])


# In[27]:


genome = pd.read_csv(
    'promotores de vigna.txt',
    sep=' ', header=None, names=['ID','Size', 'Direction', 'Promoter'])
print(type(genome), genome.shape)
genome


# In[28]:


def ecrMiner (promoter, ecr):
    searchLog = pd.DataFrame(data=None, columns=['Promoter', 'ID', 'Name', 'Motifs', 'ReverseComplement', '5l3l', '3l5l', 'Mean', 'SUM'])
    searchLog
    
    for prom in promoter.itertuples():
        print (prom.ID)
        
        for motivo in ecr.itertuples():
            mean = 0.0
            positivo = searchFT (prom.Promoter, 2000, motivo.Motifs)
            negativo = searchFT(prom.Promoter, 2000, motivo.ReverseComplement)
            hits = len(positivo+negativo)
            if hits > 0:
                mean = sum(positivo+negativo)/hits
            searchLog = searchLog.append ({'Promoter':prom.ID, 'ID':motivo.ID ,'Name':motivo.Name, 'Motifs':motivo.Motifs, 'ReverseComplement':motivo.ReverseComplement, '5l3l': positivo, '3l5l':negativo, 'Mean':mean, 'SUM':hits}, ignore_index=True)
    return searchLog


# In[ ]:


searchLog = df.filter(['ID', 'Name', 'Motifs', 'ReverseComplement'], axis=1)
searchLog['Gene'] = ""
searchLog['5l3l'] = 0
searchLog['3l5l'] = 0
searchLog['Mean'] = 0
searchLog['SUM'] = 0
searchLog

for prom in genome.itertuples():
    print (prom.ID)

    for motivo in df.itertuples():
        mean = 0.0
        positivo = searchFT (prom.Promoter, 2000, motivo.Motifs)
        negativo = searchFT(prom.Promoter, 2000, motivo.ReverseComplement)
        hits = len(positivo+negativo)
        if hits > 0:
            mean = sum(positivo+negativo)/hits
            searchLog = searchLog.append ({'Gene':prom.ID, 'ID':motivo.ID ,'Name':motivo.Name, 'Motifs':motivo.Motifs, 'ReverseComplement':motivo.ReverseComplement, '5l3l': positivo, '3l5l':negativo, 'Mean':mean, 'SUM':hits}, ignore_index=True)


# In[ ]:


#searchLog = ecrMiner(genome, df)
searchLog


# In[ ]:


searchLog.to_csv(path_or_buf='Genome Log.csv',
          sep=',', header=True, index=False, mode='w',
         columns=['Promoter', 'ID', 'Name', 'Motifs', 'ReverseComplement', '5l3l', '3l5l'])


# In[ ]:


promoter = pd.read_csv(
    'promotores de fenilpropanoides.txt',
    sep=' ', header=None, names=['ID','Size', 'Direction', 'Promoter'])
print(type(promoter), promoter.shape)
promoter


# Executando a mineração de ecrs

# In[ ]:


searchLog = ecrMiner(promoter, df)
searchLog


# Salvando o log como planilha

# In[ ]:


searchLog.to_csv(path_or_buf='Log.csv',
          sep=',', header=True, index=False, mode='w',
         columns=['Promoter', 'ID', 'Name', 'Motifs', 'ReverseComplement', '5l3l', '3l5l'])


# In[ ]:


mean = searchLog.groupby(['ID']).mean()
mean


# In[ ]:


df['Fenilpropanoides'] = list(media)
df


# In[ ]:


media = analise.groupby(['ID']).mean()
media


# In[ ]:


df['Mean'] = list(media.Mean)
df


# In[ ]:


hits = analise.groupby(['ID'])['SUM'].sum()
df['SUM'] = list(hits)
df


# In[ ]:


analise = analise + df


# In[ ]:


media = analise.groupby(['ID']).mean()
media


# In[ ]:


df['Mean_GN'] = list(media.Mean)
df


# In[ ]:


hits = analise.groupby(['ID'])['SUM'].sum()
df['SUM_GN'] = list(hits)
df

