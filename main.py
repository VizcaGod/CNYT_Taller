import numpy as np
import math

'''El sistema debe calcular la probabilidad de encontrarlo en una posición en particular'''
def superposicion (v,p):
    norma = 0
    for i in range (len(v)):        #norma
        n = abs(v[i])**2
        norma += n
    c = p*p                         #complejo al cuadradp
    return c/norma

'''El sistema si se le da otro vector Ket debe buscar la probabilidad de transitar del primer vector al segundo.'''
def prob_ket(v1,v2):

    pro_int = np.inner(v1,v2)  #saca el producto interno de los dos vectores
    norma1 = 0
    norma2= 0
    for i in range(len(v1)):    #norma del primer vector
        n = abs(v1[i])
        norma1 += n
    for i in range(len(v2)):   #norma del segundo vector
        n = abs(v2[i])
        norma2 += n
    normas = norma1*norma2      #producto de las 2 normas
    return pro_int/normas       #resultado

'''Amplitud de transición. El sistema puede recibir dos vectores y calcular la probabilidad de transitar de el uno al otro después de hacer la observación'''

def AmplitudTransition(v1, v2):

    norma1 = 0
    for i in range(len(v1)):  # norma primer vector
        n = abs(v1[i])
        norma1 += n
    norma2 = 0
    for i in range(len(v2)):  # norma segundo vector
        n = abs(v2[i])
        norma2 += n
    normas = norma1 * norma2  # producto de las 2 normas
    pro_int = np.inner(v1, v2)  # producto interno de los dos vectores
    return pro_int / normas  # resultado


'''Ahora con una matriz que describa un observable y un vector ket, el sistema revisa que la matriz sea hermitiana, y si lo es, calcula la
   media y la varianza del observable en el estado dado'''
#Funciones necesarias

def mat_vec (mat, v):                               #Acción matriz por vector
    mat = np.array(mat)
    v = np.array(v)
    matriz_vector = mat.dot(v)
    return matriz_vector

def es_hermitiana(mat):                             #Retorna si es hermitiana o no
    matConjug = np.conjugate(mat)
    return np.array_equal(mat, np.transpose(matConjug))

def producto_inner (v1, v2):
    prod = 0
    for i in range(len(v1)):
        prod += v2[i] * v1[i].conjugate()
    return prod

def normV(v):                                       #normaliza vector
    for i in range(len(v)):
        v[i] = v[i] / (np.linalg.norm(v))
    return v

'''media y la varianza del observable en el estado dado'''

def media(obs, est):
    vp_norm = np.linalg.norm(est)
    if es_hermitiana(obs) == True:                                          #Revisa que la matriz sea hermitiana
        if vp_norm != 1:
            est = normV(est)                                                #normaliza el vector, si no lo está
        media = producto_inner((mat_vec(obs, est)), est)                    #saca la media
    return media


def varianza (obs, est):
    identity_por_media = (media(obs, est)) * (np.identity(len(est)))        #Multiplicación media por la identidad
    rest =  np.array(obs) - np.array(identity_por_media)                    #Al observable se le resta la identidad por la media
    resultado = media((np.array(rest).dot(np.array(rest))), est)            #media del resultado final de la varianza
    return resultado



'''El sistema calcula los valores propios del observable y la probabilidad de que el sistema transite a alguno de los vectores propios después de la observación'''

def valores_propios(mat):
    mat = np.array(mat)                                                     #Valores propios
    return (np.linalg.eigvals(mat))

'''Se considera la dinámica del sistema. Ahora con una serie de matrices Un el sistema calcula el estado final a partir de un estado inicial'''

def vect_propios(mat):
    mat = np.array(mat)
    valoresPropios, vectoresPropios = np.linalg.eig(mat)                    #Vectores propios
    return vectoresPropios

'''Modele en su librería los problemas'''

'''4.3.1'''

def ej1 (obs):
    spin_up = np.array([1,0])                                               #Vector spinup
    accionV_Obs = mat_vec(obs, spin_up)
    result = superposicion(accionV_Obs,0)
    return result

'''4.3.2'''

def ej2 (obs):
    spin_up = np.array([1, 0])                                              #Vector spinup
    valoresProp, vectoresProp = valores_propios(obs), vect_propios(obs)
    vectoresPropios = normV(vectoresProp)
    accion_obsV = mat_vec(obs, spin_up)
    producto1 = producto_inner(accion_obsV, vectoresPropios[0])
    n1 = np.linalg.norm(producto1)                                          #Calcula uno de los vectores norma del producto
    producto2 = producto_inner(accion_obsV, vectoresPropios[1])
    n2 = np.linalg.norm(producto2)                                          #Calcula uno de los vectores norma del producto
    res = n1 * valoresProp[0] + n2 * valoresProp[1]
    return res

'''4.4.1'''

def ej3(mat1, mat2):
    mat1_C, mat2_C = (np.array(mat1)).conjugate(), (np.array(mat2)).conjugate()             #Conjugada de las matrices
    uni_1  =  np.array(mat1_C) * np.array(mat1)                                             #Unitarias de las matrices
    uni_2 = np.array(mat2_C) * np.array(mat2)
    identity = np.identity(len(mat1))                                                       #Matriz Identidad
    if (uni_1 == identity).all() and (uni_2 == identity).all():                             #Verifica que las matrices unitarias sean iguales
        mat_Prod = np.array(mat1) * np.array(mat2)
        mat_ProdUnitaria = mat_Prod.conjugate()                                             #Producto de las matrices
        if mat_ProdUnitaria == identity:
            return True
        else:
            return False
    return False

'''4.4.2'''

def ej4(mat, est):
    mat = np.array(mat)
    matDe3 = mat**3
    Acción = mat_vec(matDe3,est)                                                #Acción matriz por vector de la matriz por su vector de estado
    result = superposicion(Acción, 3)                                           #Se saca la probabilidad
    return result

#Vectores y Matrices
v = [2+1j,-1+2j, 1j, 1+0j, 3-1j, 2+0j, -2j, -2+1j, 1-3j, -1j]
v3 = [(5+2j),(-3j)]
v4 = [(5+2j),(1-1j)]
p = (10)
vec1 = [(-1j), (1)]
vec2 = [(1), (-1j)]
obs = [[0,1/2],[1/2,0]]
matriz1 =[[0, 1],[1,0]]
matriz2 =[[math.sqrt(2)/2, math.sqrt(2)/2,],[math.sqrt(2)/2,-(math.sqrt(2)/2)]]
est = [1,0,0,0]
matriz =[[0, 1/math.sqrt(2),1/math.sqrt(2),0],[1j/math.sqrt(2),0,0,1/math.sqrt(2)],[1/math.sqrt(2),0,0,1j/math.sqrt(2)],[0, 1/math.sqrt(2),-1/math.sqrt(2),0]]
mtrx = [[1,0,1,0],[0,1,1,0],[1,0,1,1], [0,1,1,0]]

def main():
    print(superposicion(v, p))
    print(prob_ket(v3, v4))
    print(AmplitudTransition(vec1, vec2))
    print(media([[1, -1j], [1j, 2]], [math.sqrt(2) / 2 + 0j, (math.sqrt(2) / 2) * 1j]))
    print(varianza([[1, -1j], [1j, 2]], [math.sqrt(2) / 2 + 0j, (math.sqrt(2) / 2) * 1j]))
    print(valores_propios(matriz1))
    print(vect_propios(matriz2))
    print(ej1(obs))
    print(ej2(obs))
    print(ej3(matriz1, matriz2))
    print(ej4(matriz, est))

main()



