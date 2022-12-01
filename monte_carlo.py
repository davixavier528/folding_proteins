import numpy as np
import matplotlib.pyplot as plt

# Constantes utilizadas
R     = 2
A     = 2.6
delta = 1

# Definir a quantidade de aminoácidos na cadeia polipetídica
qtdAA = int ( input ( "\nDigite a quantidade de aminoácidos presentes na cadeia polipetídica: " ) )

# Definir coordenadas do primeiro aminoácido
print ( "\nDigite as coordenadas ( x , y , z ) da posição inicial do primeiro aminoácido da cadeia" )
x = float ( input ( "x: " ) )
y = float ( input ( "y: " ) )
z = float ( input ( "z: " ) )

# Definir intervalo comum entre todos os aminoácidos
intervalo = float ( input ( "\nDigite o intervalo que determinará a progressão linear das coordenadas dos demais aminoácidos: " ) )

# Definir número máximo de iterações que o algoritmo de Monte Carlo poderá fazer
maxIter = int ( input ( "\nQual o número máximo de iterações que o algoritmo Monte Carlo poderá realizar: " ) )

# Trajetórias das conformações
pos_trajectory     = [ ] # Registro das posições
dist_trajectory    = [ ] # Registro das distâncias
pot_trajectory     = [ ] # Registro dos potenciais
med_pot_trajectory = [ ] # Registro das médias dos potenciais
h_trajectory       = [ ] # Registro das hamiltonianas
delta_h_trajectory = [ ] # Registro dos delta H

# Calcula a posição inicial
posInit = [ ]
count   = 1
while count <= qtdAA :
    posInit.append ( [ x , y , z ] )
    x     += intervalo
    y     += intervalo
    z     += intervalo
    count += 1
pos_trajectory.append ( posInit )

# Define as distâncias iniciais
distInit = [ ]
for i in range ( qtdAA - 1 ) :
    dist         = ( ( posInit[i][0] - posInit[i+1][0] ) ** 2 + ( posInit[i][1] - posInit[i+1][1] ) ** 2 + ( posInit[i][2] - posInit[i+1][2] ) ** 2 ) ** 0.5
    distInit.append ( dist )
dist_trajectory.append ( distInit )

# Define os potenciais iniciais
potInit = [ ]
for i in range ( qtdAA - 1 ) :
    potencial = ( delta * distInit[i] ) ** 2 + ( R / ( distInit[i] ) ** 12 ) + ( A / ( distInit[i] ) ** 6 )
    potInit.append ( potencial )
pot_trajectory.append ( potInit )
med_pot_trajectory.append ( sum ( potInit ) )

# Define a hamiltoniana inicial
h_trajectory.append ( sum ( potInit ) )

# Função para definir distância, potencial e hamiltoniana dos estados posteriores
def hamiltonianaPosterior ( coordPost ) :

    distPost = [ ]
    for i in range ( qtdAA - 1 ) :
        dist = ( ( coordPost[i][0] - coordPost[i+1][0] ) ** 2 + ( coordPost[i][1] - coordPost[i+1][1] ) ** 2 + ( coordPost[i][2] - coordPost[i+1][2] ) ** 2 ) ** 0.5
        distPost.append ( dist )
    dist_trajectory.append ( distPost )

    potPost = [ ]
    for i in range ( qtdAA - 1 ) :
        potencial = ( delta * distPost[i] ) ** 2 + ( R / ( distPost[i] ) ** 12 ) + ( A / ( distPost[i] ) ** 6 )
        potPost.append ( potencial )
    pot_trajectory.append ( potPost )
    med_pot_trajectory.append ( sum ( potPost ) )

    hPost = sum ( potPost )
    h_trajectory.append ( sum ( potPost ) )

    return hPost

# Função para mover um aminoácido aleatoriamente
def moverAA ( coords ) :
    random_aa            = np.random.randint ( 0 , 9 )
    random_coord         = list ( np.random.uniform ( min ( pos_trajectory[0][0] ) , max ( pos_trajectory[0][9] ) , 3 ) )
    coordPost            = coords.copy()
    coordPost[random_aa] = random_coord
    pos_trajectory.append ( coordPost )
    
    return coordPost

# Monte Carlo
while maxIter > 0 :
    h1      = h_trajectory [-1]
    h2      = hamiltonianaPosterior ( moverAA ( pos_trajectory[-1] ) )
    delta_h = h2 - h1
    delta_h_trajectory.append ( delta_h )
    if delta_h < 0:
        continue
    else:
        del ( pos_trajectory[-1]     )
        del ( dist_trajectory[-1]    )
        del ( pot_trajectory[-1]     )
        del ( med_pot_trajectory[-1] )
        del ( h_trajectory[-1]       )
    
    maxIter -= 1

# print ( "\nTrajetória das posições: "                   ,     pos_trajectory )
# print ( "\nTrajetória das distâncias: "                 ,    dist_trajectory )
# print ( "\nTrajetória das energias potenciais médias: " , med_pot_trajectory )
# print ( "\nTrajetória dos valores de hamiltoniana: "    ,       h_trajectory )
# print ( "\nTrajetória dos valores de delta H: "         , delta_h_trajectory )

# Gerar plot da evolução da energia potencial
x = list ( range ( len ( med_pot_trajectory ) ) )
y = med_pot_trajectory
plt.plot ( x , y)
plt.xlabel ('Decorrer da simulação')
plt.ylabel ('Energia potencial')
plt.xticks ( [ ] )
plt.show ( )

"""
D.S. Xavier (davi_xavier@icb.usp.br)
This scripts is licensed under the GNU General Public License v. 3.
Please see http://www.fsf.org/licensing/licenses/gpl.html for details.
"""