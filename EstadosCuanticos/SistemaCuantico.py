import numpy as np


def position_probablity(num_positions, amplitudes, position):
    if position >= 2 ** num_positions:
        raise ValueError("Posición invalida.")
    probability = np.abs(amplitudes[position]) ** 2
    return probability


def transition_probability(num_positions, vectorA, vectorB):
    if len(vectorA) != 2 ** num_positions or len(vectorB) != 2 ** num_positions:
        raise ValueError("Numero de amplitud no coincide con el numero de posiciones")
    producto_interno = np.vdot(vectorA, vectorB)
    probability = np.abs(producto_interno) ** 2
    return probability


def main():
    num_positions = 2
    v1 = np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0, 0])
    v2 = np.array([0, 1 / np.sqrt(2), 1 / np.sqrt(2), 0])

    punto1 = np.round((position_probablity(num_positions, v1, 1)), 2)
    punto2 = np.round((transition_probability(num_positions, v1, v2)), 2)
    print("vector estado inicial:", v1)
    print("vector estado final:", v2)
    print(punto1)
    print(punto2)



main()
