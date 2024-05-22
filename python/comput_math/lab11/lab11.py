from matplotlib import pyplot as plt
from math import e as exp


class method():
    def __init__(self, a=None, b=None, n=None, x0=None, y0=None, eps=None):
        self.a = a
        self.b = b
        self.x0 = x0
        self.y0 = y0
        self._n = n
        self.h = (self.b - self.a) / self._n
        self.eps = eps
        self._x = []
        self._y = []
        self._bn = self._n
        self._bh = self.h
        self._ex = []
        self._sy1 = []
        self._sdy = []

    @staticmethod
    def dydx(x, y):
        return (0 + 1) * exp ** (-x ** 2) + 2 * (2 - x) * y

    @staticmethod
    def f(x):
        return exp ** (-x ** 2 + 4 * x) * (.1 - (4 * x + 1) * exp ** (-4 * x) / 16 - exp ** (-4 * x) / 4)

    def _cor_y(self, i, h):
        return self.f(self._x[i])

    def _runge_y(self, i, h):
        k1 = h * self.dydx(self._x[i], self._y[i])
        k2 = h * self.dydx(self._x[i] + 0.5 * h, self._y[i] + 0.5 * k1)
        k3 = h * self.dydx(self._x[i] + 0.5 * h, self._y[i] + 0.5 * k2)
        k4 = h * self.dydx(self._x[i] + h, self._y[i] + k3)
        return self._y[i] + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    def _create_grid(self, f, n, h):
        self._x = []
        self._x.append(self.a)
        self._y = []
        self._y.append(self.y0)
        for i in range(n):
            self._x.append(self._x[i] + self.h)
        for i in range(len(self._x) - 1):
            self._y.append(f(i, h))
        return self._x, self._y

    def cor(self, n, h):
        _cor = dict()
        _cor['x'], _cor['y'] = self._create_grid(self._cor_y, self._n, self.h)
        _cor['y'].remove(self.y0)
        _cor['y'].insert(len(_cor['x']) - 1, self.f(self._x[len(_cor['x']) - 1]))
        return _cor

    def rungeK(self, n, h):
        _rung = dict()
        _rung['x'], _rung['y'] = self._create_grid(self._runge_y, n, h)
        return _rung

    def sup(self):
        self._bh = self.h
        self._bn = self._n
        g = self.rungeK(self._n, self.h)
        self.sy1 = g['y']
        maxx = 0
        while True:
            self._n *= 2
            self.h /= 2
            g = self.rungeK(self._n, self.h)
            self.sdy = g['y']
            for i in range(1, int(self._n / 2) + 1):
                m = abs(self.sy1[i] - self.sdy[i * 2])
                if m > maxx:
                    maxx = m
            if maxx < self.eps:
                self._ex = g['x']
                break
            else:
                maxx = 0
                self.sy1 = self.sdy
                self.sdy = []
        return maxx

    def out(self):

        maxx = this.sup()
        r = self.rungeK(self._n, self.h)
        c = self.cor(self._n, self.h)
        mm = []
        for i in range(this._n + 1):
            mm.append(abs(c['y'][i] - r['y'][i]))
        maxxx = max(mm)
        print('\n\n|max: ', maxx, '\n\nmax |y - dy|: ', maxxx, '| n begin: ', this._bn, '-|- n end: ', this._n,
              '| h end: ',
              this.h, '-|- h begin', this._bh, '|\n')
        print('Table:')
        print('----------------------------------------------------------------------------------')
        print("|     X     |       Y       |       dY      |     |Y - dY|     |      Runge      |")
        print('----------------------------------------------------------------------------------')
        for i in range(0, this._n + 1):
            if i % 2 == 0:

                j = int(i / 2)
                print(
                    '| {: ^9.4f} | {: ^13.9f} | {: ^13.9f} | {: ^16.14f} | {: ^15.12f} |'.format(this._ex[i], c['y'][i],
                                                                                                 r['y'][i],
                                                                                                 abs(c['y'][i] - r['y'][
                                                                                                     i]),
                                                                                                 abs(r['y'][i] -
                                                                                                     this.sy1[j])))

                print('----------------------------------------------------------------------------------')
            else:
                print('| {: ^9.4f} | {: ^13.9f} | {: ^13.9f} | {: ^16.14f} | {: ^15s} |'.format(r['x'][i], c['y'][i],
                                                                                                r['y'][i],
                                                                                                abs(c['y'][i] - r['y'][
                                                                                                    i]),
                                                                                                'None'))
                print('----------------------------------------------------------------------------------')
        plt.plot(r["x"], r["y"], label="RK_4")
        plt.plot(c["x"], c["y"], label="Result")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    a, b = 0, 10
    x_0, y_0 = 1, 10
    n = 10
    eps = .01
    this = method(a, b, n, x_0, y_0, eps)
    this.out()
