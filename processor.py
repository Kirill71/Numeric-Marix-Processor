from enum import IntEnum


class DeterminantResolver:

    def __init__(self, matrix):
        if matrix.get_cols() != matrix.get_rows():
            raise ArithmeticError
        self._matrix = matrix

    def determine(self):
        return self.__determine(self._matrix)

    @staticmethod
    def generate_cofactor_matrix(result, senior_minor, cofactor_i=0, cofactor_j=0):

        if cofactor_j >= result.get_cols():
            cofactor_j = 0
            cofactor_i += 1

        if cofactor_i >= result.get_rows():
            return

        minor = DeterminantResolver.calculate_minor(senior_minor, cofactor_i, cofactor_j)
        sign = 1 if (cofactor_i + cofactor_j) % 2 == 0 else -1
        result[cofactor_i][cofactor_j] = sign * DeterminantResolver.__determine(minor)

        DeterminantResolver.generate_cofactor_matrix(result, senior_minor, cofactor_i, cofactor_j + 1)

    @staticmethod
    def __determine(senior_minor, cofactor_i=0, cofactor_j=0, is_toplevel_minor=False):
        if senior_minor.get_rows() <= 2 and senior_minor.get_cols() <= 2:
            return DeterminantResolver.determaine_2_on_2(senior_minor)

        # it means that we process all junior minors from 0 minor to n
        if is_toplevel_minor:
            cofactor_j = 0

        junior_minor = DeterminantResolver.calculate_minor(senior_minor, cofactor_i, cofactor_j)
        sign = 1 if (cofactor_i + cofactor_j) % 2 == 0 else -1

        if cofactor_j >= senior_minor.get_cols():
            return 0.0

        return sign * senior_minor[cofactor_i][cofactor_j] * \
               DeterminantResolver.__determine(junior_minor, cofactor_i, cofactor_j, True) + \
               DeterminantResolver.__determine(senior_minor, cofactor_i,
                                                cofactor_j + 1)

    @staticmethod
    def calculate_minor(matrix, cofactor_i, cofactor_j):
        minor = Matrix(matrix.get_rows() - 1, matrix.get_cols() - 1)
        new_i = 0
        new_j = 0
        if cofactor_j >= matrix.get_cols():
            cofactor_j = matrix.get_cols() - 1
        if cofactor_i >= matrix.get_rows():
            cofactor_i = matrix.get_rows() - 1

        for i in range(matrix.get_rows()):
            for j in range(matrix.get_cols()):
                if i != cofactor_i and j != cofactor_j:
                    minor[new_i][new_j] = matrix[i][j]
                    if new_j == minor.get_cols() - 1:
                        # it means: new line will be started
                        new_i += 1
                        new_j = 0
                    else:
                        new_j += 1
        return minor

    @staticmethod
    def determaine_2_on_2(minor):
        if minor.get_cols() == 1 and minor.get_rows() == 1:
            return minor[0][0]

        return minor[0][0] * minor[1][1] - minor[0][1] * minor[1][0]


class TransposeAlg(IntEnum):
    MAIN = 1
    SIDE = 2
    VERT = 3
    HORZ = 4


class Matrix:

    def __init__(self, n, m):
        self._n = int(n)
        self._m = int(m)
        if self._n <= 0 or self._m <= 0:
            raise IndexError
        self._matrix = [[0.0 for _ in range(self._m)] for _ in range(self._n)]

    def get_rows(self):
        return self._n

    def get_cols(self):
        return self._m

    def __getitem__(self, key):
        return self._matrix[key]

    def __setitem__(self, key, value):
        self._matrix[key] = value

    def __iadd__(self, other):
        return self + other

    def __add__(self, other):
        if self._n == other.get_rows() \
                and self._m == other.get_cols():
            result = Matrix(self._n, self._m)
            for i in range(self._n):
                for j in range(self._m):
                    result[i][j] = self[i][j] + other[i][j]
            return result
        else:
            raise ArithmeticError

    def __mul__(self, other):
        if type(other) == float or type(other) == int:
            result = Matrix(self._n, self._m)
            for i in range(self._n):
                for j in range(self._m):
                    result[i][j] = self[i][j] * other
            return result

        if type(other) == Matrix:
            result = Matrix(self._n, other.get_cols())
            for i in range(self._n):
                for j in range(other.get_cols()):
                    for k in range(self._m):
                        result[i][j] += self[i][k] * other[k][j]

            return result

        return None

    def __imul__(self, other):
        types = {float, int, Matrix}
        if type(other) in types:
            return self * other
        return None

    def __rmul__(self, other):
        return self * other

    def transpose(self, alg: TransposeAlg):
        n = self._n
        m = self._m
        if alg == TransposeAlg.MAIN:
            for i in range(n):
                for j in range(m):
                    if i > j:
                        temp = self[i][j]
                        self[i][j] = self[j][i]
                        self[j][i] = temp
        elif alg == TransposeAlg.SIDE:
            for i in range(n):
                for j in range(m):
                    if i + j < m:
                        temp = self[i][j]
                        self[i][j] = self[n - j - 1][m - i - 1]
                        self[n - j - 1][m - i - 1] = temp
        elif alg == TransposeAlg.VERT:
            for i in range(n):
                for j in range(m // 2):
                    temp = self[i][j]
                    self[i][j] = self[i][m - j - 1]
                    self[i][m - j - 1] = temp
        else:
            for i in range(n // 2):
                for j in range(m):
                    temp = self[i][j]
                    self[i][j] = self[n - i - 1][j]
                    self[n - i - 1][j] = temp

    def determine(self):
        resolver = DeterminantResolver(self)
        return resolver.determine()

    def generate_cofactor_matrix(self):
        matrix = Matrix(self._n, self._m)
        resolver = DeterminantResolver(self)
        resolver.generate_cofactor_matrix(matrix, self)
        return matrix


class Processor:

    def __init__(self):
        self._matrix1 = None
        self._matrix2 = None

    @staticmethod
    def prepare(n, m):
        matrix = Matrix(n, m)
        for i in range(int(n)):
            matrix[i] = [float(number) for number in input().split()]
        return matrix

    @staticmethod
    def help():
        print('1. Add matrices')
        print('2. Multiply matrix by a constant')
        print('3. Multiply matrices')
        print('4. Transpose matrix')
        print('5. Calculate a determinant')
        print('6. Inverse matrix')
        print('0. Exit')

    @staticmethod
    def transpose_help():
        print('1. Main diagonal')
        print('2. Side diagonal')
        print('3. Vertical line')
        print('4. Horizontal line')

    @staticmethod
    def print(matrix: Matrix):
        for row in matrix:
            print(*row)
        print()

    @staticmethod
    def get_transpose_alg(alg):
        if alg == 1:
            return TransposeAlg.MAIN
        elif alg == 2:
            return TransposeAlg.SIDE
        elif alg == 3:
            return TransposeAlg.VERT
        else:
            return TransposeAlg.HORZ

    def run(self):
        while True:
            Processor.help()
            number = int(input('Your choice:'))
            if number == 0:
                break
            try:
                if number == 1:
                    n, m = input('Enter size of first matrix:').split()
                    print('Enter first matrix:')
                    self._matrix1 = Processor.prepare(int(n), int(m))

                    n1, m1 = input('Enter size of second matrix:').split()
                    print('Enter second matrix:')
                    self._matrix2 = Processor.prepare(int(n1), int(m1))
                    print('The result is:')
                    Processor.print(self._matrix1 + self._matrix2)

                elif number == 2:
                    n, m = input('Enter size of first matrix:').split()
                    print('Enter first matrix:')
                    self._matrix1 = Processor.prepare(int(n), int(m))

                    constant = int(input('Enter constant:'))
                    print('The result is:')
                    Processor.print(self._matrix1 * constant)

                elif number == 3:
                    n, m = input('Enter size of first matrix:').split()
                    print('Enter first matrix:')
                    self._matrix1 = Processor.prepare(int(n), int(m))

                    n1, m1 = input('Enter size of second matrix:').split()
                    print('Enter second matrix:')
                    self._matrix2 = Processor.prepare(int(n1), int(m1))
                    print('The result is:')
                    Processor.print(self._matrix1 * self._matrix2)
                elif number == 4:
                    Processor.transpose_help()
                    alg = int(input('Your choice:'))
                    n, m = input('Enter size of first matrix:').split()
                    print('Enter first matrix:')
                    self._matrix1 = Processor.prepare(int(n), int(m))
                    self._matrix1.transpose(Processor.get_transpose_alg(alg))
                    print('The result is:')
                    Processor.print(self._matrix1)
                elif number == 5:
                    n, m = input('Enter size of first matrix:').split()
                    print('Enter first matrix:')
                    self._matrix1 = Processor.prepare(int(n), int(m))
                    print(self._matrix1.determine())
                elif number == 6:
                    n, m = input('Enter size of first matrix:').split()
                    print('Enter first matrix:')
                    self._matrix1 = Processor.prepare(int(n), int(m))
                    det = self._matrix1.determine()
                    if det == 0.0:
                        print("This matrix doesn't have an inverse.")
                    else:
                        cofactor_matrix = self._matrix1.generate_cofactor_matrix()
                        cofactor_matrix.transpose(TransposeAlg.MAIN)
                        Processor.print((1.0 / det) * cofactor_matrix)
                else:
                    pass

            except ArithmeticError:
                print('ERROR')


if __name__ == '__main__':
    processor = Processor()
    processor.run()
