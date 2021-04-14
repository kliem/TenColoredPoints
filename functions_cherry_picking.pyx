# cython: cdivision=True
# distutils: depends = KPartiteKClique/kpkc.cpp KPartiteKClique/kpkc.h
# distutils: sources = KPartiteKClique/kpkc.cpp
# distutils: include_dirs = KPartiteKClique
# distutils: extra_compile_args=-O3 -march=native -std=c++11
# distutils: language = c++


from sage.misc.cachefunc import cached_method, cached_function
from sage.structure.sage_object cimport SageObject
from itertools import combinations, permutations, product, chain
from sage.ext.memory_allocator cimport MemoryAllocator
from sage.rings.all import ZZ

cimport cython

from sage.data_structures.bitset cimport Bitset
from sage.data_structures.bitset_base cimport bitset_next
from libc.stdio                     cimport FILE, fopen, fclose, fwrite, fread

from libcpp cimport bool

cdef extern from "KPartiteKClique/kpkc.h":
    cdef cppclass KPartiteKClique:
        KPartiteKClique()
        void init(bool **, int n_vertices, int* first_per_part, int k) except +
        bool next() except +
        const int* k_clique()

def color_iterator():
    """
    Iterate through all possible colors of points for the colorful Tverberg with 10 points.
    """
    # 3,3,3,1
    for i in range(10):
        d = (i,)
        other = tuple(j for j in range(10) if not j == i)
        other_new = other[1:]
        k = other[0]
        for l,m in combinations(other_new, 2):
            a = (k,l,m)
            other2 = tuple(j for j in range(10) if not (j == i or j in a))
            n = other2[0]
            other2_new = other2[1:]
            for o,p in combinations(other2_new, 2):
                b = (n,o,p)
                other3 = tuple(j for j in range(10) if not (j ==i or j in a or j in b))
                yield (a,b,other3,d)

    # 3,3,2,2
    for four in combinations(range(10),4):
        six = tuple(i for i in range(10) if not i in four)
        a = six[0]
        five = six[1:]
        g = four[0]
        three = four[1:]
        for h in three:
            C = (g,h)
            D = tuple(i for i in three if not i == h)
            for b,c in combinations(five,2):
                A = (a,b,c)
                B = tuple(i for i in five if not i in A)
                yield (A,B,C,D)

    # 2,2,2,2,2
    a = 0
    for b in range(1,10):
        A = (a,b)
        A1 = tuple(i for i in range(1,10) if not i == b)
        c = A1[0]
        for d in A1[1:]:
            B = (c,d)
            B1 = tuple(i for i in A1 if not i in B)
            e = B1[0]
            for f in B1[1:]:
                C = (e,f)
                C1 = tuple(i for i in B1 if not i in C)
                g = C1[0]
                for h in C1[1:]:
                    D = (g,h)
                    E = tuple(i for i in C1 if not i in D)
                    yield (A,B,C,D,E)

@cached_function
def all_colors():
    """
    A tuple consisiting of a dictionary per color.

    The dictionary assigns each point to a color.
    """
    def to_dic(colors):
        dic = {i:-1 for i in range(10)}
        for i in range(len(colors)):
            for x in colors[i]:
                dic[x] = i
        return dic
    return tuple(to_dic(c) for c in color_iterator())

@cached_function
def n_colors():
    return len(all_colors())

@cached_function
def colors_for_partition(partition):
    r"""
    Given a partition.

    Return a bitset with each color set, where this partition is colorful.
    """
    cap = ((n_colors()-1)//256 + 1)*256

    return Bitset(chain(_colors_for_partition(partition), range(n_colors(),cap)) , capacity=cap)

def _colors_for_partition(partition):
    for i,color in enumerate(all_colors()):
        partition2 = [[color[x] for x in y] for y in partition]
        # Each part in the partition must contain each color at most once.
        if all(len(set(x)) == len(x) for x in partition2):
            yield i

def OrderTypesIterator():
    r"""
    Iterate through all order types in general position for 10 points.

    See http://www.ist.tugraz.at/staff/aichholzer/research/rp/triangulations/ordertypes/
    """
    filename = "/srv/public/kliem/OrderTypes/otypes10.b16"
    with open(filename, "rb") as f:
        for i in range(14309547):
            data = f.read(40)
            yield tuple((int.from_bytes(data[4*i:4*i+2], "little"), int.from_bytes(data[4*i+2:4*i+4], "little")) for i in range(10))

def OrderType(n):
    r"""
    Return the n-th order type.

    See ``OrderTypesIterator``.
    """
    filename = "/srv/public/kliem/OrderTypes/otypes10.b16"
    assert n in range(14309547)
    with open(filename, "rb") as f:
        f.seek(40*n)
        data = f.read(40)
        return tuple((int.from_bytes(data[4*i:4*i+2], "little"), int.from_bytes(data[4*i+2:4*i+4], "little")) for i in range(10))

def pseudo_order_type_iterator(int n=10, path="/storage/mi/kliem/OrderSets/10"):
    r"""
    Iterate over all pseudo order sets generated from ``get_all_order_sets.py``.

    Each pseudo order set is given as a 120 ``char``, for the signitarues of the chirotopes.
    """
    cdef FILE *fp
    fp = fopen(path.encode('utf-8'), "r")
    if (fp==NULL):
        raise IOError("cannot open file {}".format(path))
    combs = tuple(combinations(range(n),3))
    cdef size_t end = len(combs)
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef char* choices = <char*> mem.allocarray(len(combs), sizeof(char))
    while fread(choices, end, sizeof(char), fp):
        yield tuple(choices[i] for i in range(len(combs)))

    fclose(fp)

def check_all_colors(start=0, end=2**20):
    r"""
    Check for each order set the possible counter example color indices, if any.
    """
    return _check_all_colors(start, end, OrderTypesIterator())

def check_all_colors_pseudo(start=0, end=2**20):
    r"""
    Check for each pseoduo order set the possible counter example color indices, if any.
    """
    return _check_all_colors(start, end, pseudo_order_type_iterator(10, "/srv/public/kliem/OrderSets/10"))

def _check_all_colors(start, end, iterator):
    """
    See above.
    """
    for i,O1 in enumerate(iterator):
        if i < start:
            continue
        if i >= end:
            break
        O = OrientedMatroid(O1)
        if i % 1000 == 0:
            print("currently doing: ", i)
        try:
            _ = poss_color_finder(O)
        except AssertionError:
            # There appears to be a counter example.
            raise AssertionError("found counterexample at index {}".format(i))

def orientation(p1, p2, p3):
    # to find the orientation of
    # an ordered triplet (p1,p2,p3)
    # function returns the following values:
    # 0 : Colinear points
    # -1 : Clockwise points
    # 1 : Counterclockwise
    # Taken from https://www.geeksforgeeks.org/orientation-3-ordered-points/
    val = ((p2[1] - p1[1]) * (p3[0] - p2[0])) - \
           ((p2[0] - p1[0]) * (p3[1] - p2[1]))

    if (val > 0):

        # Clockwise orientation
        return -1
    elif (val < 0):

        # Counterclockwise orientation
        return 1
    else:

        # Colinear orientation
        return 0

@cached_function
def combs():
    return tuple(combinations(range(10), 3))

@cython.final
cdef class OrientedMatroid(SageObject):
    r"""
    An Oriented Matroid with 10 points and methods that help to analyze
    the Tverberg situation.
    """
    cdef dict __dict__
    cdef int ***dic
    cdef MemoryAllocator __mem__

    def __init__(self, data):
        r"""
        Input can be one of the following:

        - complete dictionary of all orientations
        - an index of the OrderType
        - 120 integers in `\{1,-1\}` for each orientation of the chirotopes in order of
          ``itertools.combinations(range(10), 3)``
        - 10 points in the plane
        """
        self.__mem__ = MemoryAllocator()

        # Initialize self.dic.
        self.dic = <int***> self.__mem__.allocarray(10, sizeof(int**))
        cdef int i,j,k
        for i in range(10):
            self.dic[i] = <int**> self.__mem__.allocarray(10, sizeof(int*))
            for j in range(10):
                self.dic[i][j] = <int*> self.__mem__.allocarray(10, sizeof(int))


        if isinstance(data, dict):
            for i,j,k in data:
                self.dic[i][j][k] = data[i,j,k]
            return
        if data in ZZ:
            # Input is an index of the OrderType
            points = OrderType(data)
        elif len(data) == 120:
            for a,(i,j,k) in enumerate(combs()):
                o = data[a]
                self.dic[i][j][k] = o
                self.dic[j][k][i] = o
                self.dic[k][i][j] = o
                self.dic[i][k][j] = -o
                self.dic[j][i][k] = -o
                self.dic[k][j][i] = -o
            return
        else:
            points = data

        for i,j,k in combinations(range(10), 3):
            o = orientation(points[i],points[j],points[k])
            assert o != 0
            self.dic[i][j][k] = o
            self.dic[j][k][i] = o
            self.dic[k][i][j] = o
            self.dic[i][k][j] = -o
            self.dic[j][i][k] = -o
            self.dic[k][j][i] = -o

    def _repr_(self):
        r"""
        Maybe not pretty, but at least something.
        """
        return {(i,j,k): self.dic[i][j][k] for i,j,k in combinations(range(10), 3)}

    @cached_method
    def n_intersection_points(self):
        return len(self.intersection_points())

    @cached_method
    def intersection_points(self):
        r"""
        A tuple containing all the intersection points.

        The intersection points with the least number of possibilities first.
        """
        def my_len(i):
            (a,b), (c,d) = i
            return len(self.possibilities_for_intersection(a,b,c,d))
        return tuple(sorted(tuple(self._intersection_points()), key=my_len))

    def _intersection_points(self):
        r"""
        Iterate over all pairs (a,b), (c,d)
        such that the lines ab and cd intersect
        and there could be a Tverberg partition with this intersection point.
        """
        cdef int a,b,c,d
        cdef tuple y
        cdef int counter, counter2
        for a,b in combinations(range(10), 2):
            y = tuple(i for i in range(a+1,10) if i != b)
            for c,d in combinations(y, 2):

                # Check if the points intersect.
                if self.dic[a][b][c] != self.dic[a][b][d] and self.dic[c][d][a] != self.dic[c][d][b]:

                    # Check if there are between 2 and 4 points above and below each line (of the remaining 6).
                    # If is only one points below a line, then we cannot partition the 6 points into 3,3
                    # where the convex hull of each three contains ab \cap cd.
                    counter = 0
                    counter2 = 0
                    for i in range(10):
                        if i not in (a,b,c,d):
                            if self.dic[a][b][i] == 1:
                                counter += 1
                            if self.dic[c][d][i] == 1:
                                counter2 += 1
                    if 2 <= counter <= 4 and 2 <= counter2 <= 4:
                        yield (a,b), (c,d)

    @cached_method
    def possibilities_for_intersection(self, int a, int b, int c, int d):
        return tuple(self._possibilities_for_intersection(a,b,c,d))

    def _possibilities_for_intersection(self, int a, int b, int c, int d):
        r"""
        Iterate over all possible dictionaries for ab \cap cd.

        Possible means, we can extend this oriented matroid by the intersection ab \cap cd,
        located according to the dictionary.
        """
        # We put each of the remaining 6 points in the correct region induced by ab and cd.
        # Note that the opposite regions are 0,3 and 1,2.
        cdef int i,j,x,y,w,v,k,l, n_oposed, one, two, one_val, two_val, i1, i2, j1, j2

        sections, section_dic, _ = self.sections(a,b,c,d)

        # For each points i,j distinct and distinct from a,b,c,d
        # the orientation of (i,j,(ab \cap cd)) is clear unless i and j are in opposite regions.
        # So we iterate over all i,j in opposite regions.
        dic = {}
        z = ((a,b), (c,d))
        for x,y in ((0,3),(1,2)):
            for i in sections[x]:
                for j in sections[y]:

                    if self.dic[i][j][a] == self.dic[i][j][b]:
                        # The line i,j doesn't intersect the line a,b.
                        dic[(i,j,z)] = self.dic[i][j][a]
                    elif self.dic[i][j][d] == self.dic[i][j][c]:
                        # The line i,j doesn't intersect the line c,d.
                        dic[(i,j,z)] = self.dic[i][j][c]
                    else:
                        # The line i,j intersects a,b and c,d.
                        w,v = (0,3) if (x,y) == (1,2) else (1,2)

                        # We can choose the orientation of i,j,z iff the points
                        # in the sections w,v lie on the correct sides.

                        if self.dic[a][b][c] == 1:
                            # In this case section 1 lies between a and c etc.

                            # Case 1: x == 0, y == 3:
                            # In this case our situation looks like this:
                            #
                            #       c
                            #    w  |  j
                            # a-----+----b
                            #    i  |  v
                            #       d

                            # This means that i-j-k oriented different than i-j-a for any k in w
                            # implies that i-j-z is oriented different like i-j-k.

                            # Likewise i-j-l oriented the same as i-j-a for any l in v
                            # means that i-j-z is oriented just like i-j-a.

                            # Case 2: x == 1, y == 2:
                            # In this case our situation looks like this:
                            #
                            #       c
                            #    i  |  v
                            # a-----+----b
                            #    w  |  j
                            #       d

                            # This means that i-j-k oriented different than i-j-a for any k in w
                            # implies that i-j-z is oriented different like i-j-k.

                            # Likewise i-j-l oriented the same as i-j-a for any l in v
                            # means that i-j-z is oriented just like i-j-a.

                            for k in sections[w]:
                                if self.dic[i][j][k] != self.dic[i][j][a]:
                                    # the point k is "in between" the intersection point
                                    # and i,j
                                    dic[(i,j,z)] = self.dic[i][j][k]
                                    break
                            else:
                                for l in sections[v]:
                                    if self.dic[i][j][l] == self.dic[i][j][a]:
                                        # the point l is "in between" the intersection point
                                        # and i,j
                                        dic[(i,j,z)] = self.dic[i][j][l]
                                        break
                        else:
                            # In this case section 2 lies between a and c etc.

                            # Case 1: x == 0, y == 3:
                            # In this case our situation looks like this:
                            #
                            #       d
                            #    j  |  w
                            # a-----+----b
                            #    v  |  i
                            #       c

                            # This means that i-j-k oriented different than i-j-b for any k in w
                            # implies that i-j-z is oriented different like i-j-k.

                            # Likewise i-j-l oriented the same as i-j-b for any l in v
                            # means that i-j-z is oriented just like i-j-b.

                            # Case 2: x == 1, y == 2:
                            # In this case our situation looks like this:
                            #
                            #       d
                            #    v  |  i
                            # a-----+----b
                            #    j  |  w
                            #       c

                            # This means that i-j-k oriented different than i-j-b for any k in w
                            # implies that i-j-z is oriented different like i-j-k.

                            # Likewise i-j-l oriented the same as i-j-b for any l in v
                            # means that i-j-z is oriented just like i-j-ab
                            for k in sections[w]:
                                if self.dic[i][j][k] != self.dic[i][j][b]:
                                    # the point k is "in between" the intersection point
                                    # and i,j
                                    dic[(i,j,z)] = self.dic[i][j][k]
                                    break
                            else:
                                for l in sections[v]:
                                    if self.dic[i][j][l] == self.dic[i][j][b]:
                                        # the point l is "in between" the intersection point
                                        # and i,j
                                        dic[(i,j,z)] = self.dic[i][j][l]
                                        break

                    if (i,j,z) not in dic:
                        # In any remaining case we have full choice.
                        dic[(i,j,z)] = 0

        n_opposed = len(dic.keys())
        opposed = list(dic.keys())

        # For each possibility we write down the implications for each orientation.
        implications = {i: {1: [], -1: []} for i in range(n_opposed)}
        for (one, two) in combinations(range(n_opposed), 2):
            one_val = dic[opposed[one]]
            i1, j1,_ = opposed[one]
            two_val = dic[opposed[two]]
            i2, j2,_ = opposed[two]
            if dic[(i1,j1,z)] or dic[(i2,j2,z)]:
                # We skip the trivial cases from above.
                # So both of them must have a choice.
                continue

            if i1 != i2 and j1 != j2:
                # All 4 points i1,i2,j1,j2 are distinct.
                if section_dic[i1] != section_dic[i2]:
                    # This case has already been handled above.
                    continue

                x = self.dic[i1][j1][i2]
                y = self.dic[i1][j1][j2]
                if x == y:
                    # i2,j2 are both on the same side of i1,j1.

                    # i1 and i2 in same region; (j1 and j2 in same region)
                    # This means that if the intersection point is on the other side of (i1,j1), than the same must hold for (i2,j2).
                    implications[one][-x].append(two)

                x = self.dic[i2][j2][i1]
                y = self.dic[i2][j2][j1]
                if x == y:
                    # i1,j1 are both on the same side of i2,j2.

                    # i1 and i2 in same region; (j1 and j2 in same region)
                    # That means that if the intersection point is on the other side of (i2,j2), than the same must hold for (i1,j1).
                    implications[two][-x].append(one)

            elif i1 != i2:
                # This means that j1 == j2.
                x = self.dic[i1][j1][i2]
                # If the intersection point and i2 are on different sides of i1-j1 this implies that
                # the orientation i1-j1-intersection is the same as i2-j2-intersection (j1 == j2).
                implications[one][-x].append(two)
            else:
                # This means that i1 == i2.
                y = self.dic[i1][j1][j2]
                # If the intersection point and j2 are on different sides of i1-j1 this implies that
                # the orientation i1-j1-intersection is the same as i2-j2-intersection (i1 == i2).
                implications[one][-y].append(two)

        for chosen in product([1,-1], repeat=n_opposed):
            if any(chosen[i] == -dic[(opposed[i])] for i in range(n_opposed)):
                # There are some intersections, for which we cannot choose.
                continue
            for i in range(n_opposed -1):
                val = chosen[i]

                if not all(chosen[k] == val for k in  implications[i][val]):
                    break
            else:
                yield {opposed[k]: chosen[k] for k in range(n_opposed)}

    @cached_method
    def _n_poss(self, i):
        r"""
        Number of possibilities for intersection number i.
        """
        (a,b), (c,d) = self.intersection_points()[i]
        return len(self.possibilities_for_intersection(a,b,c,d))

    @cached_method
    def sections(self, int a, int b, int c, int d):
        r"""
        The 2 lines ab and cd give 4 quadrants/sections.

        Sort the other points in those sections.

        Return the sections, an inverse dictionary and a tuple of the other points.
        """
        sections = [[], [], [], []]
        sections_op = {}
        others = tuple(i for i in range(10) if not i in (a,b,c,d))
        for i in others:
            count = 0
            # i belongs to an odd region iff a,b,i is oriented counter-clock-wise.
            count += int(self.dic[a][b][i] == 1)
            # i belongs to region >= 2 iff c,d,i is oriented counter-clock-wise.
            count += 2*int(self.dic[c][d][i] == 1)
            sections[count].append(i)
            sections_op[i] = count
        return sections, sections_op, others

    @cached_method
    def obtain_dic(self, int i, int j):
        r"""
        Obtain the j-th possibility for the i-th intersection point.
        """
        (a,b), (c,d) = self.intersection_points()[i]
        dic = self.possibilities_for_intersection(a,b,c,d)[j]
        return (dic, a,b,c,d)

    cdef int ****__consistency_cache__

    cpdef inline bint check_for_consistency_cached(self, int i, int j, int k, int l):
        r"""
        See below.
        """
        cdef int a,b,c,d,a1,b1,c1,d1,n

        # Initialize cache.
        if self.__consistency_cache__ is NULL:
            n = self.n_intersection_points()
            self.__consistency_cache__ = <int****> self.__mem__.allocarray(n, sizeof(int ***))
            for a in range(n):
                self.__consistency_cache__[a] = <int***> self.__mem__.allocarray(self._n_poss(a), sizeof(int **))
                for b in range(self._n_poss(a)):
                    self.__consistency_cache__[a][b] = <int**> self.__mem__.allocarray(n, sizeof(int *))
                    for c in range(n):
                        self.__consistency_cache__[a][b][c] = <int*> self.__mem__.calloc(self._n_poss(c), sizeof(int))

        cdef dict dic1, dic2

        if self.__consistency_cache__[i][j][k][l]:
            return self.__consistency_cache__[i][j][k][l] -1
        dic1, a,b,c,d = self.obtain_dic(i,j)
        dic2, a1,b1,c1,d1 = self.obtain_dic(k,l)
        foo = self.check_for_consistency(dic1, dic2, i,j,k,l, a,b,c,d,a1,b1,c1,d1)
        self.__consistency_cache__[i][j][k][l] = int(foo) + 1
        return foo

    cpdef inline bint check_for_consistency(self, dict dic1, dict dic2, int i, int j, int k, int l, int a, int b, int c, int d, int a1, int b1, int c1, int d1):
        r"""
        Check if the two choices for intersection points are consistent.

        If they are intersection points on the same section, they both define
        what is above and below this intersection point.
        Possibly those definitions are not consistent.

        Otherwise, we say that the are consistent.
        """
        if a == -1:
            intersections2 = []
            for dic in (dic1, dic2):
                for w in dic:
                    intersections2.append(w[2])
                    break
            intersections = tuple(intersections2)
            ((a,b), (c,d)), ((a1,b1), (c1,d1)) = intersections

        cdef int* above_below1
        cdef int* above_below2
        cdef int one, two
        cdef int index

        if (a == a1 and b == b1) or (a == c1 and b == d1):
            # (a,b) is the common section.
            if i >= 0:
                above_below1 = self.get_above_below_cached(i, j, a, b)
                one = above_below1[55]
                above_below2 = self.get_above_below_cached(k, l, a, b)
                two = above_below2[55]
            else:
                above_below1x, one = self.get_above_below(dic1, a, b)
                above_below2x, two = self.get_above_below(dic2, a, b)
                # TODO initialize as int*
                assert above_below1 is not NULL

        elif (c == a1 and d == b1) or (c == c1 and d == d1):
            # (c,d) is the common section.
            if i >= 0:
                above_below1 = self.get_above_below_cached(i, j, c, d)
                one = above_below1[55]
                above_below2 = self.get_above_below_cached(k, l, c, d)
                two = above_below2[55]
            else:
                above_below1x, one = self.get_above_below(dic1, c, d)
                above_below2x, two = self.get_above_below(dic2, c, d)
                # TODO initialize as int*
                assert above_below1 is not NULL
        else:
            # No common section
            return True


        if above_below2[one] == 1:
            # one lies above two
            if above_below1[two] == 1:
                return False
            #assert two in below1
            # The situation is not consistent, if there is anything above ``one`` that is below ``two``
            for index in range(55):
                if above_below1[index] == 1 and above_below2[index] == -1:
                    return False
            return True
        else:
            # one lies below two
            if above_below1[two] == -1:
                return False
            # The situation is not consistent, if there is anything below ``one`` that is above ``two``
            for index in range(55):
                if above_below1[index] == -1 and above_below2[index] == 1:
                    return False
            return True

        return True

    cdef int**** __above_below_cache__

    cdef inline int* get_above_below_cached(self, int i, int j, int a, int b):
        cdef int x,y,z, n

        # Initialize cache
        if self.__above_below_cache__ is NULL:
            n = self.n_intersection_points()
            self.__above_below_cache__ = <int****> self.__mem__.allocarray(n, sizeof(int ***))
            for x in range(n):
                self.__above_below_cache__[x] = <int***> self.__mem__.allocarray(self._n_poss(x), sizeof(int **))
                for y in range(self._n_poss(x)):
                    self.__above_below_cache__[x][y] = <int**> self.__mem__.calloc(10, sizeof(int *))

        cdef dict dic
        cdef tuple foo
        if self.__above_below_cache__[i][j][a] is NULL:
            # Note that ``b`` is determined by a, as the dictionary i,j already inscribes an intersection of two lines.
            self.__above_below_cache__[i][j][a] = <int*> self.__mem__.allocarray(56, sizeof(int))

            dic = self.obtain_dic(i,j)[0]
            foo = self.get_above_below(dic, a,b)
            for x in range(55):
                self.__above_below_cache__[i][j][a][x] = foo[0][x]
            self.__above_below_cache__[i][j][a][55] = foo[1]
        return self.__above_below_cache__[i][j][a]

    cpdef inline tuple get_above_below(self, dict dic, int a, int b):
        r"""
        Determine which sections are above the intersection point and which are below.

        Return a triple ``above``, ``below``, ``(c,d)``.

        ``dic is a dictionary`` containing the placement of a intersection point of ``(a,b)``.

        We define ``a`` to be on top.

        All tuples are sections are given, as (x,y) such that (x,y,a) is oriented counter-clockwise.
        """
        # Let c1,d1 be the other section of the intersection point.
        for k in dic:
            (i,j), (k,l) = k[2]
            break
        intersection = ((i,j), (k,l))
        c1,d1 = (i,j) if a == k else (k,l)

        # Reorient c1,d1 such that c,d,a is counter-clockwise.
        if self.dic[c1][d1][a] == 1:
            c,d = c1,d1
        else:
            c,d = d1,c1

        above_below = [0]*55

        for index,(i1,j1) in enumerate(combinations(range(10),2)):
            if i1 in (a,b) or j1 in (a,b):
                continue
            if (i1,j1) == (c1,d1):
                # This is the section defining the intersection.
                index_cd = index
                continue
            if self.dic[a][b][i1] == self.dic[a][b][j1] or self.dic[i1][j1][a] == self.dic[i1][j1][b]:
                # The lines don't even intersect.
                continue

            # Orient i,j,a counter-clock-wise.
            if self.dic[i1][j1][a] == -1:
                j = i1
                i = j1
            else:
                i = i1
                j = j1

            # If the location of i-j is up to a choice,
            # we look it up in the dictionary.
            # Note that ``a`` is defined to be on top and that i,j,a is counter-clockwise.
            if (i,j,intersection) in dic:
                val = dic[(i,j,intersection)]
                if val == 1:
                    above_below[index] = -1
                else:
                    above_below[index] = 1
            elif (j,i,intersection) in dic:
                val = -dic[(j,i,intersection)]
                if val == 1:
                    above_below[index] = -1
                else:
                    above_below[index] = 1

            elif i in (c,d):
                if self.dic[c][d][j] == 1:
                    # j lies above c,d, (i is c or d) therefore
                    # i-j intersects a,b, above c,d.
                    above_below[index] = 1
                else:
                    above_below[index] = -1
            elif j in (c,d):
                if self.dic[c][d][i] == 1:
                    # As above.
                    above_below[index] = 1
                else:
                    above_below[index] = -1
            else:
                # In all remaining cases i,j must both lie above or both below c,d,
                # otherwise the location of the section i,j with respect to the intersection point (a,b,c,d)
                # would be up to a choice and would therefore be needed to be in the dictionary.
                assert self.dic[c][d][i] == self.dic[c][d][j]
                if self.dic[c][d][i] == 1:
                    above_below[index] = 1
                else:
                    above_below[index] = -1
        return tuple(above_below), index_cd

    def shortest_final_realization(self, *args, start=0, stop=1000):
        r"""
        Currently not used.
        """
        rainbows = [self.rainbow_partitions_cached(i,args[i]) for i in range(len(args))]
        c = rainbows[0]
        for x in rainbows:
            c = c.intersection(x)
        choices = list(c)
        assert choices
        minimum = 1000
        for seed in range(start, min(stop,len(choices))):
            choice = choices[seed]

            def is_interesting(int_points, color_dic):
                (a,b, c,d) = int_points
                left = [color_dic[x] for x in range(10) if not x in (a,b,c,d)]
                return color_dic[a] != color_dic[b] and color_dic[c] != color_dic[d] and all(2 >= sum(1 for x in left if x == i) for i in range(3))
            new_args = []
            for i,arg in enumerate(args):
                int_points = self.obtain_dic(i,arg)[1:5]
                if is_interesting(int_points, all_colors()[choice]):
                    new_args.append(self.obtain_dic(i,arg)[0])
            minimum = min(len(new_args), minimum)

        return minimum

    def final_realization(self, *args, verbose=True, seed=0):
        r"""
        Find a realization with SCIP including the exact location of intersection points.

        Currently not needed.

        Seed is to select one of the possible colors. Irrelevant intersection points are then ignored to reduce complexity.
        """
        rainbows = [self.rainbow_partitions_cached(i,args[i]) for i in range(len(args))]
        c = rainbows[0]
        for x in rainbows:
            c = c.intersection(x)
        choices = list(c)
        assert choices
        if seed >= len(choices):
            raise ValueError("seed to large")
        choice = choices[seed]

        def is_interesting(int_points, color_dic):
            (a,b, c,d) = int_points
            left = [color_dic[x] for x in range(10) if not x in (a,b,c,d)]
            return color_dic[a] != color_dic[b] and color_dic[c] != color_dic[d] and all(2 >= sum(1 for x in left if x == i) for i in range(3))
        new_args = []
        for i,arg in enumerate(args):
            int_points = self.obtain_dic(i,arg)[1:5]
            if is_interesting(int_points, all_colors()[choice]):
                new_args.append(self.obtain_dic(i,arg)[0])

        print(len(new_args), choice, all_colors()[choice])
        return self.find_realization(*new_args, verbose=verbose, color_choice=choice)

    def find_realization(self, *args, verbose=False, color_choice=None):
        r"""
        Find a realization with SCIP, currently not needed.
        """
        from pyscipopt import Model, SCIP_PARAMSETTING
        if len(args) > 0 and not isinstance(args[0], dict):
            new_args = []
            for i,arg in enumerate(args):
                new_args.append(self.obtain_dic(i,arg)[0])
            args = new_args

        '''
        for dic1,dic2 in combinations(args,2):
            if not self.check_for_consistency(dic1,dic2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1):
                return None
        '''
        model = Model("OM")
        if not verbose:
            model.hideOutput()
        model.setPresolve(SCIP_PARAMSETTING.OFF)
        xval = [model.addVar(name="x{}".format(i), ub=100, vtype="C") for i in range(10)]
        yval = [model.addVar(name="x{}".format(i), ub=100, vtype="C") for i in range(10)]

        intersections = []
        for dic in args:
            for k in dic:
                intersections.append(k[2])
                break

        xint = [model.addVar(name="x1{}".format(i), ub=100, vtype="C") for i in intersections]
        yint = [model.addVar(name="y1{}".format(i), ub=100, vtype="C") for i in intersections]

        if color_choice is not None:
            colors = all_colors()[color_choice]

        def is_rainbow_choice(*args):
            if color_choice is None:
                return True
            return len(args) == len(set(colors[i] for i in args))

        bounds = {}
        for i,j,k in combinations(range(10), 3):
            o = self.dic[i][j][k]
            if o == 1:
                bounds[(i,j,k)] = model.addCons((yval[j] - yval[i])*(xval[k] - xval[j]) - (xval[j] - xval[i])*(yval[k] - yval[j]) <= -0.1)
            else:
                bounds[(i,j,k)] = model.addCons((yval[j] - yval[i])*(xval[k] - xval[j]) - (xval[j] - xval[i])*(yval[k] - yval[j]) >= 0.1)

        for i,dic in enumerate(args):
            (a,b), (c,d) = intersections[i]
            intersection = intersections[i]
            j = a
            k = b
            bounds[(intersection,1,1)] = model.addCons((yval[j] - yint[i])*(xval[k] - xval[j]) - (xval[j] - xint[i])*(yval[k] - yval[j]) == 0)
            j = c
            k = d
            bounds[(intersection,1,2)] = model.addCons((yval[j] - yint[i])*(xval[k] - xval[j]) - (xval[j] - xint[i])*(yval[k] - yval[j]) == 0)
            for s,t,intersection in dic:
                if not is_rainbow_choice(s,t):
                    # We don't care about how those not colorful lines are chosen.
                    continue
                val = dic[(s,t,intersection)]
                if val == 1:
                    j = s
                    k = t
                else:
                    j = t
                    k = s
                bounds[(intersection,s,t,1)] = model.addCons((yval[j] - yint[i])*(xval[k] - xval[j]) - (xval[j] - xint[i])*(yval[k] - yval[j]) <= -0.1)

        #model.setObjective(xval[0], "minimize")

        model.optimize()
        self.model = model
        if not model.getStatus() == "optimal":
            raise ValueError("appears to be infeasible")
        else:
            return tuple((model.getVal(xval[i]), model.getVal(yval[i])) for i in range(10))

    @cached_method
    def partitions_one(self):
        return tuple(self._partitions_one())

    @cached_method
    def partitions_one_bitset(self):
        r"""
        Return a bitset that has each color index set to 1, for which there IS Tverberg partition of type 3,3,3,1
        """
        cap = ((n_colors()-1)//256 + 1)*256
        cdef Bitset foo = Bitset(range(n_colors(),cap), capacity = cap)
        for part in self._partitions_one():
            foo = foo.union(colors_for_partition(part))
        return foo

    def _partitions_one(self):
        r"""
        Obtain all tverberg partitions that are determined by oriented matroid on 10 points, i.e. partitions with one point in three triangles.
        """
        for x in range(10):
            others = tuple(i for i in range(10) if i != x)
            a = others[0]
            for b,c in combinations(others[1:], 2):
                if self.dic[a][b][x] == self.dic[b][c][x] == self.dic[c][a][x]:
                    # x in (convex hull of) a,b,c
                    others1 = tuple(i for i in range(10) if not i in (x,a,b,c))
                    d = others1[0]
                    for e,f in combinations(others1[1:], 2):
                        # x in d,e,f
                        if self.dic[d][e][x] == self.dic[e][f][x] == self.dic[f][d][x]:
                            g,h,i = tuple(i for i in range(10) if not i in (x,a,b,c,d,e,f))
                            if self.dic[g][h][x] == self.dic[h][i][x] == self.dic[i][g][x]:
                                # x in g,h,i
                                yield ((x,), (a,b,c), (d,e,f), (g,h,i))

    def partitions_two(self, dict dic, int a, int b, int c, int d):
        """
        Yield all Tverberg partitions given by a certain dictionary.
        """
        sections, sections_op, others = self.sections(a,b,c,d)
        intersection = ((a,b), (c,d))
        cdef int x,y,w,e,f,g

        def is_triangle(x,y,z):
            """
            Check if x,y,z is a triangle containing the intersection of ab and cd.
            """
            for i,j in ((0,3),(1,2)):
                if (any(w in sections[i] for w in (x,y,z)) and any(w in sections[j] for w in (x,y,z))):
                    k,l = (0,3) if i == 1 else (1,2)
                    # Now the opposite sections i,j contain one point of x,y,z.

                    if any(w in sections[k] or w in sections[l] for w in (x,y,z)):
                        # the vertices are in three sections, so it is all up to one line
                        v1,v2 = tuple(w for w in (x,y,z) if not (w in sections[k] or w in sections[l]))
                        v3 = tuple(w for w in (x,y,z) if not w in (v1,v2))[0]
                        if (v1,v2,intersection) in dic:
                            return dic[(v1,v2,intersection)] == self.dic[v1][v2][v3]
                        return dic[(v2,v1,intersection)] == self.dic[v2][v1][v3]
                    else:
                        # the vertices are only in two sections
                        one = tuple(w for w in (x,y,z) if w in sections[i])
                        two = tuple(w for w in (x,y,z) if w in sections[j])
                        v2,v3 = one if len(one) == 2 else two
                        v1 = one[0] if len(one) == 1 else two[0]
                        if (v1,v3,intersection) in dic:
                            val = dic[(v1,v3,intersection)]
                        else:
                            val = -dic[(v3,v1,intersection)]
                        if (v1,v2,intersection) in dic:
                            val2 = dic[(v1,v2,intersection)]
                        else:
                            val2 = -dic[(v2,v1,intersection)]
                        return val != val2
            return False

        x = others[0]
        for y,z in combinations(others[1:],2):
            if is_triangle(x,y,z):
                e,f,g = tuple(i for i in range(10) if not i in (a,b,c,d,x,y,z))
                if is_triangle(e,f,g):
                    yield ((a,b), (c,d), (x,y,z), (e,f,g))

    @cached_method
    def rainbow_partitions_cached(self, i,j):
        r"""
        Return a bitset for intersection i,j that contains a 0 for each color,
        that has a Tverberg partion for this choice (1 otherwise).
        """
        cdef int a,b,c,d
        cdef dict dic
        dic, a,b,c,d = self.obtain_dic(i,j)
        return self.rainbow_partitions(dic, a,b,c,d)

    def rainbow_partitions(self, dict dic, int a, int b, int c, int d):
        cdef Bitset start = self.partitions_one_bitset()
        for part in self.partitions_two(dic,a,b,c,d):
            start = start.union(colors_for_partition(part))
        return ~start

    @cached_method
    def poss_graph(self):
        """
        currently not used
        """
        from sage.graphs.graph import Graph
        G = Graph()
        ls = self.intersection_points()
        n = len(ls)
        for i,j in combinations(range(n), 2):
            x = ls[i]
            y = ls[j]
            for k,a in enumerate(self.possibilities_for_intersection(*x[0], *x[1])):
                for l,b in enumerate(self.possibilities_for_intersection(*y[0], *y[1])):
                    if self.check_for_consistency_cached(i,k,j,l):
                        G.add_edge((i,k), (j,l))
        V = G.vertices()
        for v in V:
            rainbows = self.rainbow_partitions_cached(*v)
            for l in rainbows:
                G.add_edge(v, (n, l))
        return G

# -------- Main algorithm

def poss_color_finder(OrientedMatroid O):
    r"""
    Determine all colors for which this OrientedMatroid has a solution:

    There exists some setup with the intersections so that the colors don't have a Tverberg Partition.
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef Bitset somea  # type awareness of cython
    cdef int i, j

    cdef int k = O.n_intersection_points() + 1
    cdef int* first_per_part = <int*> mem.allocarray(k, sizeof(int))
    cdef int counter = 0
    first_per_part[0] = 0
    for i in range(k-1):
        counter += O._n_poss(i)
        first_per_part[i+1] = counter

    cdef int num_colors = n_colors()

    cdef int n = counter + num_colors

    cdef bool ** incidences = <bool **> mem.allocarray(n, sizeof(bool *))
    for i in range(n):
        incidences[i] = <bool *> mem.calloc(n, sizeof(bool))

    # Obstructions
    # Obtain for each possible constellation the colors for which there is no Tverberg Partition.

    cdef int offset_colors = first_per_part[k-1]
    cdef int offset_i, ind_i
    cdef long ind_color

    for i in range(k-1):
        offset_i = first_per_part[i]
        for ind_i in range(O._n_poss(i)):
            somea = O.rainbow_partitions_cached(i, ind_i)
            ind_color = bitset_next(somea._bitset, 0)
            while ind_color != -1:
                incidences[offset_i + ind_i][offset_colors + ind_color] = True
                incidences[offset_colors + ind_color][offset_i + ind_i] = True
                ind_color = bitset_next(somea._bitset, ind_color + 1)

    cdef int l

    # All connections
    # For each constellation we store, which other constellations are still possible, if this is picked.
    cdef int ind_j, offset_j
    for i in range(k-1):
        offset_i = first_per_part[i]
        for ind_i in range(O._n_poss(i)):
            for j in range(i):
                offset_j = first_per_part[j]
                for ind_j in range(O._n_poss(j)):
                    # There is an arc between those vertices if and only if i,ind_i and j,ind_j are consistent.
                    if O.check_for_consistency_cached(i, ind_i, j, ind_j):
                        incidences[offset_i + ind_i][offset_j+ ind_j] = True
                        incidences[offset_j + ind_j][offset_i+ ind_i] = True

    cdef KPartiteKClique K = KPartiteKClique()
    K.init(incidences, n, first_per_part, k)
    foo = K.next()
    if foo:
        # There is a counter example.
        raise AssertionError

    return []
