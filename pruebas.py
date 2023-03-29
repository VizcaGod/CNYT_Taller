import unittest
import math
import main as tcb


class Test_libtcb(unittest.TestCase):

    def test_superposicion(self):
        self.assertEqual(tcb.superposicion([2+1j,-1+2j, 1j, 1+0j, 3-1j, 2+0j, -2j, -2+1j, 1-3j, -1j], 10), 2.1739130434782608)

    def test_prob_ket(self):
        self.assertEqual(tcb.prob_ket([(5+2j),(-3j)],[(5+2j),(1-1j)]),(0.3157124386617289+0.29817285873607724j))

    def test_AmplitudTransition(self):
        self.assertEqual(tcb.AmplitudTransition([(-1j), (1)],[(1), (-1j)]),-0.5j)

    def test_media(self):
        self.assertEqual(tcb.media([[1, -1j], [1j, 2]], [math.sqrt(2)/2+0j, (math.sqrt(2)/2)*1j]),(2.5000000000000004+0j))

    def test_varianza(self):
        self.assertEqual(tcb.varianza([[1, -1j], [1j, 2]], [math.sqrt(2) / 2, (math.sqrt(2) / 2) * 1j]), (0.25+0j))

    '''def test_valores_propios(self):
            arr1 = np.array([[0, 1],[1,0]])
            print (arr1)
            self.assertEqual(tcb.valores_propios(np.all(arr1)),[ 1. -1.])'''

    def test_ej1(self):
        self.assertEqual(tcb.ej1([[0,1/2],[1/2,0]]),0.0)

    def test_ej2(self):
        self.assertEqual(tcb.ej2([[0,1/2],[1/2,0]]),-0.019337567297406433)

    def test_ej3(self):
        matriz1 = [[0, 1], [1, 0]]
        matriz2 = [[math.sqrt(2) / 2, math.sqrt(2) / 2, ], [math.sqrt(2) / 2, -(math.sqrt(2) / 2)]]
        self.assertEqual(tcb.ej3(matriz1, matriz2),False)

    def test_ej4(self):
        est = [1, 0, 0, 0]
        matriz = [[0, 1 / math.sqrt(2), 1 / math.sqrt(2), 0], [1j / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
                  [1 / math.sqrt(2), 0, 0, 1j / math.sqrt(2)], [0, 1 / math.sqrt(2), -1 / math.sqrt(2), 0]]
        self.assertEqual(tcb.ej4(matriz, est),36.000000000000014)

if __name__ == '__main__':
    unittest.main()

