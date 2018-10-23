import random
import numpy as np
from bisect import bisect_left


class RandomGen(object):
    """ Weighted Random Sampling
    There are three algorithms available:
    -. 'RouletteWheel', essentially a binary search of the cumulative distribution.
    -. 'BiasedCoin', throw coins with varied odd
    -. 'AliasMethod', specially designed for randon weighted sampling.

    Usage:
    values = [...]
    probabilities = [...]
    generator =  RandomGen(values, probabilities, <algo>)
    # where <algo> can be 'RouletteWheel', 'BiasedCoin', or 'AliasMethod'.
    choice = generator.next()
    """
    def __init__(self, vals, probs, algo='RouletteWheel'):
        # validity check
        if not vals or abs(sum(probs) - 1) > 1e-5:
            raise ValueError('Invalid input.')

        self.vals = vals
        self.probs = probs
        self.n = len(vals)

        self.algo = algo

        # helper function to choose from the three algorithms
        def choose_algorithm(x):
            return {
                'RouletteWheel':
                    [self.init_roulete_wheel, self.gen_roulete_wheel],
                'BiasedCoin':
                    [self.init_biased_coin, self.gen_biased_coin],
                'AliasMethod':
                    [self.init_alias_method, self.gen_alias_method]
            }[x]

        self.init_method, self.gen_method = choose_algorithm(algo)
        self.init_method()

    def next(self):
        return self.gen_method()

    # Algorithm for a Roulette Wheel
    def init_roulete_wheel(self):
        self.A = np.cumsum(self.probs)

    def gen_roulete_wheel(self):
        return self.vals[bisect_left(self.A, random.random())]

    # Algorithm for Biased Coin
    def init_biased_coin(self):
        pass

    def gen_biased_coin(self):
        mass = 1
        for i in range(self.n):
            if random.random() < self.probs[i] / mass:
                return self.vals[i]
            else:
                mass -= self.probs[i]

    # Algorithm for Alias Method
    def init_alias_method(self):
        self.alias = [None] * self.n
        self.prob = [None] * self.n

        small = []
        large = []

        self.probs = [self.n * e for e in self.probs]
        [small.append(i) if e < 1 else large.append(i)
            for i, e in enumerate(self.probs)]

        while small and large:
            l = small.pop()
            g = large.pop()

            self.prob[l] = self.probs[l]
            self.alias[l] = g

            self.probs[g] = (self.probs[g] + self.probs[l] - 1)
            if self.probs[g] < 1:
                small.append(g)
            else:
                large.append(g)

        def changeToOne(l, e):
            l[e] = 1
        [changeToOne(self.prob, large.pop())
            for e in reversed(large)]
        [changeToOne(self.prob, small.pop())
            for e in reversed(small)]

    def gen_alias_method(self):
        i = int(np.floor(random.random() * self.n))
        return self.vals[i] \
            if random.random() < self.prob[i] \
            else self.vals[self.alias[i]]

if __name__ == '__main__':
    print("Weighted Random Selection.");
    pass
