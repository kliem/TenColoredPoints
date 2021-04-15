from sage.all import *
from itertools import combinations, product, chain
from sage.ext.memory_allocator cimport MemoryAllocator
from libc.string            cimport memcmp, memcpy, memset
from libc.stdio                     cimport FILE, fopen, fclose, fwrite, fread


class PseudoOrderSet(SageObject):
    def __init__(self, dic, n):
        for i, j, k in combinations(range(n), 3):
            assert (i, j, k) in dic
            val = dic[i, j, k]
            dic[i, k, j] = -val
            dic[j, i, k] = -val
            dic[k, j, i] = -val
            dic[j, k, i] = val
            dic[k, i, j] = val
        self.dic = dic
        self.n = n

    def _repr_(self):
        return {i: self.dic[i] for i in combinations(range(self.n), 3)}.__repr__()

    @cached_method
    def lambda_matrix(self):
        M = Matrix([[0]*self.n]*self.n)
        for i, j, k in combinations(range(self.n), 3):
            if self.dic[i, j, k] == 1:
                M[i, j] += 1
                M[j, k] += 1
                M[k, i] += 1
            else:
                M[j, i] += 1
                M[k, j] += 1
                M[i, k] += 1
        return M

    def is_minimal(self):
        M = self.lambda_matrix()
        Mt = self.lambda_matrix().transpose()
        for dic, dic2 in self.relabels():
            M1 = Matrix([[M[dic2[j], dic2[i]] for i in range(self.n)] for j in range(self.n)])
            if M1 < M:
                return False
        for dic, dic2 in self.relabels_op():
            M1 = Matrix([[Mt[dic2[j], dic2[i]] for i in range(self.n)] for j in range(self.n)])
            if M1 < M:
                return False
        return True

    def __hash__(self):
        return (tuple(self.dic[a, b, c] for a, b, c in combinations(range(self.n), 3))).__hash__()

    def relabels(self):
        M = self.lambda_matrix()
        correct = [0] + list(range(self.n - 1))
        for i in range(1, self.n):
            a = M.row(i)
            if sorted(a) == correct:
                yield {j: a[j] + 1 if j != i else 0
                       for j in range(self.n)}, {a[j] + 1 if j != i else 0: j for j in range(self.n)}

    def relabels_op(self):
        M = self.lambda_matrix().transpose()
        correct = [0] + list(range(self.n - 1))
        for i in range(self.n):
            a = M.row(i)
            if sorted(a) == correct:
                yield {j: a[j] + 1 if j != i else 0
                       for j in range(self.n)}, {a[j] + 1 if j != i else 0: j for j in range(self.n)}

    def last_is_valid(self):
        for i in range(self.n - 1):
            for one in range(self.n - 1):
                if one == i:
                    continue
                if all(self.dic[x, i, one] == self.dic[x, i, self.n - 1] for x in range(self.n - 1) if x not in (one, i)) or \
                        all(self.dic[x, i, one] != self.dic[x, i, self.n - 1] for x in range(self.n - 1) if x not in (one, i)):
                    break
            else:
                return False
        return True

        for i in range(self.n - 1):
            other = tuple(x for x in range(self.n - 1) if x is not i)
            likewise = []
            opposite = []
            for one in other:
                if all(self.dic[x, i, one] == self.dic[x, i, self.n - 1] for x in other if x is not one):
                    likewise.append(one)
                elif all(self.dic[x, i, one] != self.dic[x, i, self.n - 1] for x in other if x is not one):
                    opposite.append(one)
            if len(likewise) + len(opposite) != 2:
                return False
            continue
            all_others = tuple(x for x in range(self.n - 1) if x not in likewise + opposite + [i])
            for a, b in combinations(all_others, 2):
                vals = [self.dic[a, i, x] for x in likewise] \
                        + [self.dic[i, b, x] for x in likewise] \
                        + [-self.dic[a, i, x] for x in opposite] \
                        + [-self.dic[b, i, x] for x in opposite]
                if all(x == self.dic[a, b, i] for x in vals):
                    if not self.dic[a, b, i] == self.dic[a, b, self.n - 1]:
                        # This section is outside the permitted range.
                        return False

        return True

    def all_regions(self, goal=None):
        if goal is None:
            goal = self.n + 1
        lines = tuple(x + (self.n,) for x in combinations(range(self.n), 2))
        choice = [0]*len(lines)
        pos = 0
        dic = copy(self.dic)
        while True:
            if pos == -1:
                break
            if 0 in lines[pos] and choice[pos] == 0:
                choice[pos] = 1
            if choice[pos] == 0:
                choice[pos] = 1
                dic[lines[pos]] = 1
            elif choice[pos] == 1:
                choice[pos] = -1
                dic[lines[pos]] = -1
            else:
                pos -= 1
                continue
            if pos == len(lines) - 1:
                x = PseudoOrderSet(dic, self.n+1)
                if not x.last_is_valid():
                    continue
                if goal == self.n + 1:
                    if x.is_minimal():
                        yield x
                else:
                    for y in x.all_regions(goal=goal):
                        yield y
            else:
                pos += 1
                choice[pos] = 0


cdef struct impli:
    int typ
    int max_index
    int other
    int length_same
    int length_opposite
    int** same
    int** opposite


def pseudo_order_type_writer(int n, path):
    cdef int i, j, k
    combs = tuple(combinations(range(n), 3))
    combs_inv = {comb: i for i, comb in enumerate(combs)}
    implications_helper = [[] for _ in range(len(combs))]
    for positive, neg in all_implications(n):
        maxi = ()
        prev_max = ()
        for x, y in chain(positive, neg):
            if y > maxi:
                if maxi > prev_max:
                    prev_max = maxi
                maxi = y
            if maxi > y > prev_max:
                prev_max = y
            if maxi > x > prev_max:
                prev_max = x
        max_index = combs_inv[maxi]
        prevmax_index = combs_inv[prev_max]

        # Read like this: if all pairs in same are the same and all pairs in opposite are opposite,
        # then max_index must be likewise/opposite to other.
        typ = 0  # 1 is likewise, -1 is opposite
        other = None
        same = []
        opposite = []

        for x, y in positive:
            if y == maxi:
                other = combs_inv[x]
                typ = -1
            else:
                same.append((combs_inv[x], combs_inv[y]))

        for x, y in neg:
            if y == maxi:
                other = combs_inv[x]
                typ = 1
            else:
                opposite.append((combs_inv[x], combs_inv[y]))

        implications_helper[prevmax_index].append((typ, max_index, other, same, opposite))

    cdef int minimalpos = combs_inv[1, 2, 3]
    cdef int end = len(combs)
    cdef int pos = minimalpos
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef char** choices = <char**> mem.allocarray(end, sizeof(char*))
    for i in range(end):
        choices[i] = <char*> mem.allocarray(end, sizeof(char))
    for i in range(minimalpos):
        choices[pos][i] = -1
        choices[pos-1][i] = -1
    for i in range(minimalpos, end):
        choices[pos][i] = 0
        choices[pos-1][i] = 0

    cdef impli** implications = <impli**> mem.allocarray(end, sizeof(impli*))
    cdef int *n_implications = <int*> mem.allocarray(end, sizeof(int))
    for i in range(end):
        implications[i] = <impli*> mem.allocarray(len(implications_helper[i]), sizeof(impli))
        n_implications[i] = len(implications_helper[i])
        for j in range(n_implications[i]):
            typ, max_index, other, same, opposite = implications_helper[i][j]
            implications[i][j].typ = typ
            implications[i][j].max_index = max_index
            implications[i][j].other = other
            implications[i][j].length_same = len(same)
            implications[i][j].length_opposite = len(opposite)
            implications[i][j].same = <int**> mem.allocarray(len(same), sizeof(int*))
            implications[i][j].opposite = <int**> mem.allocarray(len(opposite), sizeof(int*))
            for k in range(len(same)):
                implications[i][j].same[k] = <int*> mem.allocarray(2, sizeof(int))
                implications[i][j].same[k][0] = same[k][0]
                implications[i][j].same[k][1] = same[k][1]
            for k in range(len(opposite)):
                implications[i][j].opposite[k] = <int*> mem.allocarray(2, sizeof(int))
                implications[i][j].opposite[k][0] = opposite[k][0]
                implications[i][j].opposite[k][1] = opposite[k][1]

    cdef bint *skips = <bint*> mem.calloc(end, sizeof(bint))
    cdef int designated
    cdef impli implication
    for i in range(end):
        skips[i] = False

    cdef int **lambda_matrix = <int**> mem.allocarray(n, sizeof(int*))
    for i in range(n):
        lambda_matrix[i] = <int*> mem.allocarray(n, sizeof(int))

    # For each combinations we save which entries of the lambda matrix should be set, if the entry is positive.
    cdef int *lambda_matrix_entries = <int*> mem.allocarray(3*end, sizeof(int))
    for i, (a, b, c) in enumerate(combs):
        lambda_matrix_entries[i*3] = a
        lambda_matrix_entries[i*3+1] = b
        lambda_matrix_entries[i*3+2] = c

    cdef int foo

    cdef int* new_order = <int*> mem.allocarray(n, sizeof(int))

    cdef FILE* fp
    fp = fopen(path.encode('utf-8'), "w")
    if (fp==NULL):
        raise IOError("cannot open file {}".format(path))

    while pos >= minimalpos:
        if not skips[pos]:
            if choices[pos][pos] == 0:
                choices[pos][pos] = 1
            elif choices[pos][pos] == 1:
                memcpy(choices[pos], choices[pos-1], end*sizeof(char))
                choices[pos][pos] = -1
            else:
                pos -= 1
                while skips[pos] and pos >= minimalpos:
                    pos -= 1
                continue
        for i in range(n_implications[pos]):
            implication = implications[pos][i]

            if (all(choices[pos][implication.same[j][0]] == choices[pos][implication.same[j][1]]
                    for j in range(implication.length_same)) and
                all(choices[pos][implication.opposite[j][0]] != choices[pos][implication.opposite[j][1]]
                    for j in range(implication.length_opposite))):
                designated = choices[pos][implication.other]*implication.typ
                if choices[pos][implication.max_index] == 0:
                    choices[pos][implication.max_index] = designated
                elif choices[pos][implication.max_index] != designated:
                    # Not a valid choice.
                    while skips[pos] and pos >= minimalpos:
                        pos -= 1
                    break
        else:
            # Valid choice so far.
            if pos == end - 1:
                for i in range(n):
                    memset(lambda_matrix[i], 0, n*sizeof(int))

                for i in range(end):
                    if choices[pos][i] == 1:
                        lambda_matrix[lambda_matrix_entries[i*3]][lambda_matrix_entries[i*3+1]] += 1
                        lambda_matrix[lambda_matrix_entries[i*3+1]][lambda_matrix_entries[i*3+2]] += 1
                        lambda_matrix[lambda_matrix_entries[i*3+2]][lambda_matrix_entries[i*3]] += 1
                    else:
                        lambda_matrix[lambda_matrix_entries[i*3+2]][lambda_matrix_entries[i*3+1]] += 1
                        lambda_matrix[lambda_matrix_entries[i*3]][lambda_matrix_entries[i*3+2]] += 1
                        lambda_matrix[lambda_matrix_entries[i*3+1]][lambda_matrix_entries[i*3]] += 1

                # See, if the lambda matrix is minimal.
                for i in range(1, n):
                    new_order[n-1] = -1
                    for j in range(n):
                        if j == i:
                            new_order[0] = i
                        else:
                            new_order[lambda_matrix[i][j] + 1] = j

                    if new_order[n-1] == -1:
                        # the corresponding point does not lie on the boundary of the convex hull
                        continue

                    for k in range(n):
                        for j in range(n):
                            foo = lambda_matrix[new_order[j]][new_order[k]] - lambda_matrix[j][k]
                            if foo:
                                break
                        else:
                            continue
                        break
                    if foo < 0:
                        # Turns out our matrix is not the smallest one.
                        break
                else:
                    # Now the transposed.
                    for i in range(n):
                        new_order[n-1] = -1
                        for j in range(n):
                            if j == i:
                                new_order[0] = i
                            else:
                                new_order[lambda_matrix[j][i] + 1] = j

                        if new_order[n-1] == -1:
                            # the corresponding point does not lie on the boundary of the convex hull
                            continue

                        for k in range(n):
                            for j in range(n):
                                foo = lambda_matrix[new_order[k]][new_order[j]] - lambda_matrix[j][k]
                                if foo:
                                    break
                            else:
                                continue
                            break
                        if foo < 0:
                            # Turns out our matrix is not the smallest one.
                            break

                    else:
                        # yield {comb: choices[pos][i] for i,comb in enumerate(combs)}
                        # yield tuple(choices[pos][i] for i in range(end))
                        fwrite(choices[pos], end, sizeof(char), fp)
                '''
                x = PseudoOrderSet({comb: choices[pos][i] for i,comb in enumerate(combs)}, n)
                if x.is_minimal():
                    yield x
                '''
                while skips[pos] and pos >= minimalpos:
                    pos -= 1
            else:
                pos += 1
                memcpy(choices[pos], choices[pos-1], end*sizeof(char))
                skips[pos] = bool(choices[pos][pos] != 0)

    fclose(fp)


def pseudo_order_type_iterator(int n, path):
    cdef FILE* fp
    fp = fopen(path.encode('utf-8'), "r")
    if (fp == NULL):
        raise IOError("cannot open file {}".format(path))
    combs = tuple(combinations(range(n), 3))
    cdef size_t end = len(combs)
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef char* choices = <char*> mem.allocarray(len(combs), sizeof(char))
    while fread(choices, end, sizeof(char), fp):
        yield tuple(choices[i] for i in range(len(combs)))

    fclose(fp)


def pseudo_order_type_set(int n, path):
    cdef FILE* fp
    fp = fopen(path.encode('utf-8'), "r")
    if (fp == NULL):
        raise IOError("cannot open file {}".format(path))
    combs = tuple(combinations(range(n), 3))
    cdef size_t end = len(combs)
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef char* choices = <char*> mem.allocarray(len(combs), sizeof(char))
    x = set()
    while fread(choices, end, sizeof(char), fp):

        x.add(tuple(choices[i] for i in range(len(combs))))

    fclose(fp)
    return x


def all_implications(n):
    return tuple(set(all_implications_iter(n)))


def all_implications_iter(n):
    def sign(a, b, c):
        if a == b or b == c or a == c:
            return 0
        if a < b < c or b < c < a or c < a < b:
            return 1
        return -1

    def sort(a, b, c):
        return tuple(sorted((a, b, c)))

    for a, b, c in product(range(n), repeat=3):
        if a == b or b == c or a == c:
            continue
        for d1, e1, f1 in combinations(range(n), 3):
            for d, e, f in Permutations((d1, e1, f1)):
                if a in (d, e, f):
                    continue
                if b in (d, e, f) and c in (d, e, f):
                    continue
                pos = []
                neg = []
                typ = None
                if sign(a, b, c)*sign(d, e, f) == 1:
                    neg.append(tuple(sorted((sort(a, b, c), sort(d, e, f)))))
                else:
                    pos.append(tuple(sorted((sort(a, b, c), sort(d, e, f)))))

                for a1, b1, c1, d1, e1, f1 in ((d, b, c, a, e, f), (e, b, c, d, a, f), (f, b, c, d, e, a)):
                    sig = sign(a1, b1, c1)*sign(d1, e1, f1)
                    if sig == 0:
                        continue
                    if sig == 1:
                        pos.append(tuple(sorted((sort(a1, b1, c1), sort(d1, e1, f1)))))
                    else:
                        neg.append(tuple(sorted((sort(a1, b1, c1), sort(d1, e1, f1)))))
                yield tuple(sorted(pos)), tuple(sorted(neg))
