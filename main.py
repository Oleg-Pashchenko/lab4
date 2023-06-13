from math import sqrt, log, exp

from prettytable import PrettyTable
def solve_minor(matrix, i, j):
    """ Найти минор элемента матрицы """
    n = len(matrix)
    minor = []
    for row in range(n):
        if row != i:
            minor_row = []
            for col in range(n):
                if col != j:
                    minor_row.append(matrix[row][col])
            minor.append(minor_row)

    return minor


def solve_det(matrix):
    """ Найти определитель матрицы """
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    det = 0
    sgn = 1
    for j in range(n):
        det += sgn * matrix[0][j] * solve_det(solve_minor(matrix, 0, j))
        sgn *= -1
    return det


def calc_s(dots, f):
    """ Найти меру отклонения """
    n = len(dots)
    x = [dot[0] for dot in dots]
    y = [dot[1] for dot in dots]

    return sum([(f(x[i]) - y[i]) ** 2 for i in range(n)])


def calc_stdev(dots, f):
    """ Найти среднеквадратичное отклонение """
    n = len(dots)

    return sqrt(calc_s(dots, f) / n)


def lin_func(dots):
    """ Линейная аппроксимация """
    data = {}

    n = len(dots)
    x = [dot[0] for dot in dots]
    y = [dot[1] for dot in dots]

    sx = sum(x)
    sx2 = sum([xi ** 2 for xi in x])
    sy = sum(y)
    sxy = sum([x[i] * y[i] for i in range(n)])

    d = solve_det([[sx2, sx],
                   [sx, n]])
    d1 = solve_det([[sxy, sx],
                    [sy, n]])
    d2 = solve_det([[sx2, sxy],
                    [sx, sy]])

    try:
        a = d1 / d
        b = d2 / d
    except ZeroDivisionError:
        return 'Linear error input'
    data['a'] = a
    data['b'] = b

    f = lambda z: a * z + b
    data['f'] = f

    data['str_f'] = "fi = a*x + b"

    data['s'] = calc_s(dots, f)
    data['title'] = 'Линейная аппроксимация'
    data['stdev'] = calc_stdev(dots, f)

    return data


def sqrt_func(dots):
    """ Квадратичная аппроксимация """
    data = {}

    n = len(dots)
    x = [dot[0] for dot in dots]
    y = [dot[1] for dot in dots]

    sx = sum(x)
    sx2 = sum([xi ** 2 for xi in x])
    sx3 = sum([xi ** 3 for xi in x])
    sx4 = sum([xi ** 4 for xi in x])
    sy = sum(y)
    sxy = sum([x[i] * y[i] for i in range(n)])
    sx2y = sum([(x[i] ** 2) * y[i] for i in range(n)])

    d = solve_det([[n, sx, sx2],
                   [sx, sx2, sx3],
                   [sx2, sx3, sx4]])
    d1 = solve_det([[sy, sx, sx2],
                    [sxy, sx2, sx3],
                    [sx2y, sx3, sx4]])
    d2 = solve_det([[n, sy, sx2],
                    [sx, sxy, sx3],
                    [sx2, sx2y, sx4]])
    d3 = solve_det([[n, sx, sy],
                    [sx, sx2, sxy],
                    [sx2, sx3, sx2y]])

    try:
        c = d1 / d
        b = d2 / d
        a = d3 / d
    except ZeroDivisionError:
        return 'Sqrt error input'
    data['c'] = c
    data['b'] = b
    data['a'] = a

    f = lambda z: a * (z ** 2) + b * z + c
    data['f'] = f

    data['str_f'] = "fi = a*x^2 + b*x + c"

    data['s'] = calc_s(dots, f)
    data['title'] = 'Квадратичная аппроксимация'
    data['stdev'] = calc_stdev(dots, f)

    return data


def cubic_func(dots):
    """ Кубическая аппроксимация """
    data = {}

    n = len(dots)
    x = [dot[0] for dot in dots]
    y = [dot[1] for dot in dots]

    sx = sum(x)
    sx2 = sum([xi ** 2 for xi in x])
    sx3 = sum([xi ** 3 for xi in x])
    sx4 = sum([xi ** 4 for xi in x])
    sx5 = sum([xi ** 5 for xi in x])
    sx6 = sum([xi ** 6 for xi in x])
    sy = sum(y)
    sxy = sum([x[i] * y[i] for i in range(n)])
    sx2y = sum([(x[i] ** 2) * y[i] for i in range(n)])
    sx3y = sum([(x[i] ** 3) * y[i] for i in range(n)])

    # print(sx, sx2, sx3, sx4, sx5, sx6)
    # print(sy, sxy, sx2y, sx3y)

    d = solve_det([[n, sx, sx2, sx3],
                   [sx, sx2, sx3, sx4],
                   [sx2, sx3, sx4, sx5],
                   [sx3, sx4, sx5, sx6]])
    d1 = solve_det([[sy, sx, sx2, sx3],
                    [sxy, sx2, sx3, sx4],
                    [sx2y, sx3, sx4, sx5],
                    [sx3y, sx4, sx5, sx6]])
    d2 = solve_det([[n, sy, sx2, sx3],
                    [sx, sxy, sx3, sx4],
                    [sx2, sx2y, sx4, sx5],
                    [sx3, sx3y, sx5, sx6]])
    d3 = solve_det([[n, sx, sy, sx3],
                    [sx, sx2, sxy, sx4],
                    [sx2, sx3, sx2y, sx5],
                    [sx3, sx4, sx3y, sx6]])
    d4 = solve_det([[n, sx, sx2, sy],
                    [sx, sx2, sx3, sxy],
                    [sx2, sx3, sx4, sx2y],
                    [sx3, sx4, sx5, sx3y]])

    try:
        a = d1 / d
        b = d2 / d
        c = d3 / d
        d = d4 / d
    except ZeroDivisionError:
        return 'Cube error input'

    data['a'] = a
    data['b'] = b
    data['c'] = c
    data['d'] = d
    data['title'] = 'Кубическая аппроксимация'
    f = lambda x: a + b * x + c * x ** 2 + d * x ** 3
    data['f'] = f

    data['str_f'] = "f(x) = a*x^3 + b*x^2 + c*x + d"

    data['s'] = calc_s(dots, f)

    data['stdev'] = calc_stdev(dots, f)

    return data


def exp_func(dots):
    """ Экспоненциальная аппроксимация """
    data = {}

    n = len(dots)
    x = [dot[0] for dot in dots]
    y = []
    for dot in dots:
        if dot[1] <= 0:
            return 'Exp error input'
        y.append(dot[1])

    lin_y = [log(y[i]) for i in range(n)]
    lin_result = lin_func([(x[i], lin_y[i]) for i in range(n)])

    a = exp(lin_result['b'])
    b = lin_result['a']
    data['a'] = a
    data['b'] = b

    f = lambda z: a * exp(b * z)
    data['f'] = f

    data['str_f'] = "fi = a*e^(b*x)"

    data['s'] = calc_s(dots, f)
    data['title'] = 'Экспоненциальная аппроксимация'
    data['stdev'] = calc_stdev(dots, f)

    return data


def log_func(dots):
    """ Логарифмическая аппроксимация """
    data = {}

    n = len(dots)
    x = []
    for dot in dots:
        if dot[0] <= 0:
            return "Значения должны быть > 0"
        x.append(dot[0])
    y = [dot[1] for dot in dots]

    lin_x = [log(x[i]) for i in range(n)]
    lin_result = lin_func([(lin_x[i], y[i]) for i in range(n)])

    a = lin_result['a']
    b = lin_result['b']
    data['a'] = a
    data['b'] = b

    f = lambda z: a * log(z) + b
    data['f'] = f

    data['str_f'] = "fi = a*ln(x) + b"

    data['s'] = calc_s(dots, f)
    data['title'] = 'Логарифмическая аппроксимация'
    data['stdev'] = calc_stdev(dots, f)

    return data


def pow_func(dots):
    """ Степенная аппроксимация """
    data = {}

    n = len(dots)
    x = []
    for dot in dots:
        if dot[0] <= 0:
            return 'Step error input'
        x.append(dot[0])
    y = []
    for dot in dots:
        if dot[1] <= 0:
            return 'Step error input'
        y.append(dot[1])

    lin_x = [log(x[i]) for i in range(n)]
    lin_y = [log(y[i]) for i in range(n)]
    lin_result = lin_func([(lin_x[i], lin_y[i]) for i in range(n)])

    a = exp(lin_result['b'])   # объяснить почему
    b = lin_result['a']  # объяснить
    data['a'] = a
    data['b'] = b

    f = lambda z: a * (z ** b)
    data['f'] = f

    data['str_f'] = "fi = a*x^b"
    data['title'] = 'Степенная аппроксимация'
    data['s'] = calc_s(dots, f)

    data['stdev'] = calc_stdev(dots, f)

    return data


inp = [[1.1, -2.5], [2.2, 4.9],
       [3.3, 8.2], [4.4, 17.4]
       ]
#print('Ввндите точки: (для завершения введите END)')
#while True:
#    i = input()
#    if i == 'END':
#        break
#    inp.append(i.split(','))
#



funcs = [lin_func, sqrt_func, cubic_func, exp_func, log_func, pow_func]
responses = []
unique_keys = []
for f in funcs:
    response = f(inp)
    if type(response) == str:
        print(response)
        continue
    if response:
        responses.append(response)
        for key in response.keys():
            if key not in unique_keys:
                if key != 'f':
                    unique_keys.append(key)
unique_keys.sort()
unique_keys.remove('title')
unique_keys.insert(0, 'title')
arr = []
top_v = []
for resp in responses:
    r = []
    for key in unique_keys:
        if key in resp.keys():
            r.append(resp[key])
        else:
            r.append('')
    arr.append(r)
    top_v.append([resp['stdev'], resp['title']])

top_v.sort()

mytable = PrettyTable()
mytable.field_names = unique_keys
mytable.add_rows(arr)
print(mytable)
best_v = 9999999999999
for i in top_v:
    if i[0] <= best_v:
        best_v = i[0]
        print('Лучший метод:', i[1])
    else:
        break
