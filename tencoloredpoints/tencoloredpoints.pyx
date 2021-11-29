# cython: cdivision=True
# distutils: depends = KPartiteKClique/kpkc.cpp KPartiteKClique/kpkc.h
# distutils: include_dirs = KPartiteKClique
# distutils: extra_compile_args=-O3 -march=native -std=c++11
# distutils: language = c++


from functools import lru_cache
from itertools import combinations, permutations, product, chain
from memory_allocator cimport MemoryAllocator

cimport cython

from libc.stdio                       cimport FILE, fopen, fclose, fwrite, fread

from libcpp cimport bool

cdef extern from "KPartiteKClique/kpkc.cpp" namespace "kpkc":
    cdef cppclass KPartiteKClique "kpkc::KPartiteKClique<kpkc::kpkc>":
        KPartiteKClique(bool **, int n_vertices, int* first_per_part, int k) except +
        const int* k_clique()

        # Do NOT wrap this in cysignals sig_on/sig_off, it has its own interrupt handling.
        # Its own interrupt handling exits at safe states, so that this can be called again
        # to continue.
        bool next() except +

    cdef cppclass Bitset:
        Bitset(int)
        int first(int start)
        void set(int)

cdef extern from "bitset2.cpp" namespace "tencoloredpoints":
    cdef cppclass Bitset2(Bitset):
        Bitset2(int)
        void clear()
        void operator=(const Bitset2&)
        void flip_inplace()
        void union_assign(const Bitset2& l, const Bitset2& r)


cdef class Bitset_wrapper:
    cdef Bitset2* ptr
    cdef int n_bits

    def __cinit__(self, int n_bits, iterable=[]):
        """
        Initialize a bitset with ``n_bits`` bits.

        The number of bits are assumed to be a multiple of 256.

        Set exactly the bits in ``iterable``.
        """
        if n_bits % 256:
            raise ValueError("n_bits must be a multiple of 256")
        self.ptr = new Bitset2(n_bits)
        self.n_bits= n_bits
        cdef int i
        self.ptr.clear()
        for i in iterable:
            if not 0 <= i < n_bits:
                raise IndexError("can't set that bit")
            self.ptr.set(i)

    def __dealloc__(self):
        del self.ptr

    cdef void clear(self):
        """
        Set all bits to zero.
        """
        self.ptr.clear()

    cdef Bitset_wrapper copy(self):
        cdef Bitset_wrapper foo = Bitset_wrapper(self.n_bits)
        foo.ptr[0] = self.ptr[0]
        return foo

    cdef void flip_inplace(self):
        r"""
        Note that one needs to consider trailing bits for this to make sense.

        In :meth:`__cinit__` we assure that we do not have trailing bits.
        """
        self.ptr.flip_inplace()

    cdef int union_assign(self, Bitset_wrapper l, Bitset_wrapper r) except -1:
        """
        Set ``self`` to be the union of ``l`` and ``r``.
        """
        if not self.n_bits == l.n_bits == r.n_bits:
            raise ValueError("the lengths of the bitsets do not agree")
        self.ptr.union_assign(l.ptr[0], r.ptr[0])

    cdef inline int set(self, int i) except -1:
        if not 0 <= i < self.n_bits:
            raise IndexError("can't set that bit")
        self.ptr.set(i)

    cdef inline int first(self, int i):
        return self.ptr.first(i)


def color_iterator():
    r"""
    Iterate over all partitions of 10 points `C_1, \dots, C_m` such that
    `3 \geq |C_1| \geq |C_2| \geq \dots \geq |C_m|` and such that `|C_{m-1}| + |C_m| > 3`.
    """
    # Type (3, 3, 3, 1)
    for i in range(10):
        D = (i,)
        a0, *other = (j for j in range(10) if j != i)
        for a12 in combinations(other, 2):
            A = a0, *a12
            b0, *other2 = (j for j in other if j not in A)
            for b12 in combinations(other2, 2):
                B = b0, *b12
                C = tuple(j for j in other2 if j not in B)
                yield (A, B, C, D)

    # Type (3, 3, 2, 2)
    for four in combinations(range(10), 4):
        c0, *three = four
        a0, *five = tuple(i for i in range(10) if i not in four)
        for c1 in three:
            C = (c0, c1)
            D = tuple(i for i in three if not i == c1)
            for a12 in combinations(five, 2):
                A = a0, *a12
                B = tuple(i for i in five if i not in A)
                yield (A, B, C, D)

    # Type (2, 2, 2, 2, 2)
    a0 = 0
    for a1 in range(1, 10):
        A = (a0, a1)
        b0, *other = (i for i in range(1, 10) if i != a1)
        for b1 in other:
            B = (b0, b1)
            c0, *other2 = (i for i in other if i != b1)
            for c1 in other2:
                C = (c0, c1)
                d0, *other3 = (i for i in other2 if i != c1)
                for d1 in other3:
                    D = (d0, d1)
                    E = tuple(i for i in other3 if i != d1)
                    yield (A, B, C, D, E)


@lru_cache(None)
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


@lru_cache(None)
def n_colors():
    return len(all_colors())


@lru_cache(None)
def colors_for_partition(partition):
    r"""
    Given a partition.

    Return a bitset with each color set, where this partition is rainbow.

    In addition set trailing bits.
    """
    cap = ((n_colors() - 1) // 256 + 1) * 256
    return Bitset_wrapper(cap, chain(_colors_for_partition(partition), range(n_colors(), cap)))


def _colors_for_partition(partition):
    """
    Yield the indices of the color partitions, where the given partition is rainbow.
    """
    for i, color in enumerate(all_colors()):
        partition2 = [[color[x] for x in y] for y in partition]
        # Each part in the partition must contain each color at most once.
        if all(len(set(x)) == len(x) for x in partition2):
            yield i


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


@lru_cache(None)
def combs():
    return tuple(combinations(range(10), 3))


@cython.final
cdef class Chirotope:
    r"""
    A Chirotope on 10 points and methods that help to analyze
    the Tverberg situation.
    """
    cdef dict __dict__
    cdef int ***chi
    cdef MemoryAllocator __mem__
    cdef object __valid_intersection_points
    cdef dict __extensions
    cdef dict _obtain_dic
    cdef dict _regions
    cdef dict _n_extensions_cache
    cdef Bitset_wrapper _partitions_one_bitset

    def __init__(self, data):
        r"""
        Input can be one of the following:

        - complete dictionary of all orientations
        - an index of the OrderType from
          http://www.ist.tugraz.at/staff/aichholzer/research/rp/triangulations/ordertypes/
        - 120 integers in `\{1,-1\}` for each orientation of the chirotopes in order of
          ``itertools.combinations(range(10), 3)``
        - 10 points in the plane
        """
        self.__mem__ = MemoryAllocator()
        self.__valid_intersection_points = None
        self.__extensions = {}
        self._obtain_dic = {}
        self._regions = {}
        self._n_extensions_cache = {}
        self._partitions_one_bitset = None

        # Initialize self.chi.
        self.chi = <int***> self.__mem__.allocarray(10, sizeof(int**))
        cdef int i,j,k
        for i in range(10):
            self.chi[i] = <int**> self.__mem__.allocarray(10, sizeof(int*))
            for j in range(10):
                self.chi[i][j] = <int*> self.__mem__.allocarray(10, sizeof(int))


        if isinstance(data, dict):
            for i,j,k in data:
                self.chi[i][j][k] = data[i,j,k]
            return
        is_integer = isinstance(data, int)
        try:
            from sage.rings.all import ZZ
            if data in ZZ:
                is_integer = True
        except ModuleNotFoundError:
            pass
        if is_integer:
            # Input is an index of the OrderType
            points = OrderType(data)
        elif len(data) == 120:
            for a,(i,j,k) in enumerate(combs()):
                o = data[a]
                self.chi[i][j][k] = o
                self.chi[j][k][i] = o
                self.chi[k][i][j] = o
                self.chi[i][k][j] = -o
                self.chi[j][i][k] = -o
                self.chi[k][j][i] = -o
            return
        else:
            points = data

        for i,j,k in combinations(range(10), 3):
            o = orientation(points[i],points[j],points[k])
            assert o != 0
            self.chi[i][j][k] = o
            self.chi[j][k][i] = o
            self.chi[k][i][j] = o
            self.chi[i][k][j] = -o
            self.chi[j][i][k] = -o
            self.chi[k][j][i] = -o

    def _repr_(self):
        r"""
        Maybe not pretty, but at least something.
        """
        return {(i,j,k): self.chi[i][j][k] for i,j,k in combinations(range(10), 3)}

    def n_valid_intersection_points(self):
        return len(self.valid_intersection_points())

    def valid_intersection_points(self):
        r"""
        A tuple containing all the intersection points.

        The intersection points with the least number of possibilities first.
        """
        if self.__valid_intersection_points:
            return self.__valid_intersection_points
        def my_len(i):
            (a,b), (c,d) = i
            return len(self.extensions(a,b,c,d))
        self.__valid_intersection_points = tuple(sorted(tuple(self._valid_intersection_points()), key=my_len))
        return self.__valid_intersection_points

    def _valid_intersection_points(self):
        r"""
        Iterate over all pairs (a,b), (c,d)
        such that the lines ab and cd intersect
        and there could be a Tverberg partition with this intersection point.

        See Lemma 2.4.

        TODO: Relabel to Lem:IntersectionPoints.
        """
        cdef int a,b,c,d
        cdef tuple y
        cdef int counter, counter2
        for a,b in combinations(range(10), 2):
            y = tuple(i for i in range(a+1,10) if i != b)
            for c,d in combinations(y, 2):

                # Check if the points intersect.
                if self.chi[a][b][c] != self.chi[a][b][d] and self.chi[c][d][a] != self.chi[c][d][b]:

                    # Check if there are between 2 and 4 points above and below each line (of the remaining 6).
                    # A triangle containing ab \cap cd must have at least one point above and below the lines ab and cd.
                    # Thus at least two points are needed, to construct two such triangle.
                    counter = 0
                    counter2 = 0
                    for i in range(10):
                        if i not in (a,b,c,d):
                            if self.chi[a][b][i] == 1:
                                counter += 1
                            if self.chi[c][d][i] == 1:
                                counter2 += 1
                    if 2 <= counter <= 4 and 2 <= counter2 <= 4:
                        yield (a,b), (c,d)

    def extensions(self, int a, int b, int c, int d):
        dic = self.__extensions
        if not (a, b, c, d) in dic:
            dic[a, b, c, d] = tuple(self._extensions(a, b, c, d))
        return dic[a, b, c, d]

    def _extensions(self, int a, int b, int c, int d):
        r"""
        Iterate over all possible extensions of the chirotope with `y := ab \cap cd`.

        However, we will only determine those orientations ``(i, j, y)``, where ``i``
        and ``j`` are in opposite regions. According to Lemma 3.5 and Proposition 3.6
        this determines all Tverberg partitons of type (3,3,2,2) in y. This might
        identify two or more extensions.

        Also, we do not check all axioms. So some extensions might not even be valid.

        TODO: Relabel to Lem:OppositeRegions and Prop:OppositeRegions.
        """
        # We put each of the remaining 6 points in the correct region induced by ab and cd.
        # i and j lie in opposite regions, if
        # ``self.chi[a][b][i] != self.chi[a][b][j]`` and
        # ``self.chi[c][d][i] != self.chi[c][d][j]``.
        # Note that the opposite regions are indexed by 0,3 and 1,2.
        cdef int i,j,x,y,w,v,k,l, n_oposed, one, two, one_val, two_val, i1, i2, j1, j2

        regions, region_dic, _ = self.regions(a,b,c,d)

        # For each points i,j distinct and distinct from a,b,c,d
        # the orientation of (i,j,(ab \cap cd)) is clear unless i and j are in opposite regions.
        # So we iterate over all i,j in opposite regions.
        dic = {}
        z = ((a,b), (c,d))
        for x,y in ((0,3),(1,2)):
            for i in regions[x]:
                for j in regions[y]:

                    if self.chi[i][j][a] == self.chi[i][j][b]:
                        # The line i,j doesn't intersect the line a,b.
                        dic[(i,j,z)] = self.chi[i][j][a]
                    elif self.chi[i][j][d] == self.chi[i][j][c]:
                        # The line i,j doesn't intersect the line c,d.
                        dic[(i,j,z)] = self.chi[i][j][c]
                    else:
                        # The line i,j intersects a,b and c,d.
                        w,v = (0,3) if (x,y) == (1,2) else (1,2)

                        # We can choose the orientation of i,j,z iff the points
                        # in the regions w,v lie on the correct sides.

                        if self.chi[a][b][c] == 1:
                            # In this case region 1 lies between a and c etc.

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

                            for k in regions[w]:
                                if self.chi[i][j][k] != self.chi[i][j][a]:
                                    # the point k is "in between" the intersection point
                                    # and i,j
                                    dic[(i,j,z)] = self.chi[i][j][k]
                                    break
                            else:
                                for l in regions[v]:
                                    if self.chi[i][j][l] == self.chi[i][j][a]:
                                        # the point l is "in between" the intersection point
                                        # and i,j
                                        dic[(i,j,z)] = self.chi[i][j][l]
                                        break
                        else:
                            # In this case region 2 lies between a and c etc.

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
                            for k in regions[w]:
                                if self.chi[i][j][k] != self.chi[i][j][b]:
                                    # the point k is "in between" the intersection point
                                    # and i,j
                                    dic[(i,j,z)] = self.chi[i][j][k]
                                    break
                            else:
                                for l in regions[v]:
                                    if self.chi[i][j][l] == self.chi[i][j][b]:
                                        # the point l is "in between" the intersection point
                                        # and i,j
                                        dic[(i,j,z)] = self.chi[i][j][l]
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
                if region_dic[i1] != region_dic[i2]:
                    # This case has already been handled above.
                    continue

                x = self.chi[i1][j1][i2]
                y = self.chi[i1][j1][j2]
                if x == y:
                    # i2,j2 are both on the same side of i1,j1.

                    # i1 and i2 in same region; (j1 and j2 in same region)
                    # This means that if the intersection point is on the other side of (i1,j1), than the same must hold for (i2,j2).
                    implications[one][-x].append(two)

                x = self.chi[i2][j2][i1]
                y = self.chi[i2][j2][j1]
                if x == y:
                    # i1,j1 are both on the same side of i2,j2.

                    # i1 and i2 in same region; (j1 and j2 in same region)
                    # That means that if the intersection point is on the other side of (i2,j2), than the same must hold for (i1,j1).
                    implications[two][-x].append(one)

            elif i1 != i2:
                # This means that j1 == j2.
                x = self.chi[i1][j1][i2]
                # If the intersection point and i2 are on different sides of i1-j1 this implies that
                # the orientation i1-j1-intersection is the same as i2-j2-intersection (j1 == j2).
                implications[one][-x].append(two)
            else:
                # This means that i1 == i2.
                y = self.chi[i1][j1][j2]
                # If the intersection point and j2 are on different sides of i1-j1 this implies that
                # the orientation i1-j1-intersection is the same as i2-j2-intersection (i1 == i2).
                implications[one][-y].append(two)

        for chosen in product([1,-1], repeat=n_opposed):
            if any(chosen[i] == -dic[(opposed[i])] for i in range(n_opposed)):
                # There are some intersections, for which we cannot choose.
                continue
            for i in range(n_opposed -1):
                val = chosen[i]

                if not all(chosen[k] == val for k in implications[i][val]):
                    break
            else:
                yield {opposed[k]: chosen[k] for k in range(n_opposed)}

    def _n_extensions(self, i):
        r"""
        Number of possibilities for intersection number i.
        """
        if i not in self._n_extensions_cache:
            (a,b), (c,d) = self.valid_intersection_points()[i]
            self._n_extensions_cache[i] = len(self.extensions(a,b,c,d))
        return self._n_extensions_cache[i]

    def regions(self, int a, int b, int c, int d):
        r"""
        The 2 lines ab and cd induce 4 regions parametrized by
        ``(self.chi[a][b][i], self.chi[c][d][i])`` for each ``i`` different from
        ``a, b, c, d``.

        Sort the other points in those regions.

        Return the regions, an inverse dictionary and a tuple of the other points
        (all points different from ``a, b, c, d``).
        """
        if (a, b, c, d) not in self._regions:
            regions = [[], [], [], []]
            regions_op = {}
            others = tuple(i for i in range(10) if not i in (a,b,c,d))
            for i in others:
                count = 0
                # i belongs to an odd region iff a,b,i is oriented counter-clock-wise.
                count += int(self.chi[a][b][i] == 1)
                # i belongs to region >= 2 iff c,d,i is oriented counter-clock-wise.
                count += 2*int(self.chi[c][d][i] == 1)
                regions[count].append(i)
                regions_op[i] = count
            self._regions[a, b, c, d] = (regions, regions_op, others)
        return self._regions[a, b, c, d]

    def obtain_dic(self, int i, int j):
        r"""
        Obtain the j-th possibility for the i-th intersection point.
        """
        if (i, j) not in self._obtain_dic:
            (a,b), (c,d) = self.valid_intersection_points()[i]
            dic = self.extensions(a,b,c,d)[j]
            self._obtain_dic[i, j] = (dic, a,b,c,d)
        return self._obtain_dic[i, j]

    cpdef inline bint check_for_consistency(self, int i, int j, int k, int l):
        r"""
        Check if the two choices for intersection points are consistent.

        If they are intersection points on the same section, they both define
        what is above and below this intersection point.
        Possibly those definitions are not consistent.

        Otherwise, we say that the are consistent.
        """
        cdef int a, b, c, d, a1, b1, c1, d1, n

        cdef dict dic1, dic2

        dic1, a,b,c,d = self.obtain_dic(i,j)
        dic2, a1,b1,c1,d1 = self.obtain_dic(k,l)
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
            n = self.n_valid_intersection_points()
            self.__above_below_cache__ = <int****> self.__mem__.allocarray(n, sizeof(int ***))
            for x in range(n):
                self.__above_below_cache__[x] = <int***> self.__mem__.allocarray(self._n_extensions(x), sizeof(int **))
                for y in range(self._n_extensions(x)):
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
        if self.chi[c1][d1][a] == 1:
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
            if self.chi[a][b][i1] == self.chi[a][b][j1] or self.chi[i1][j1][a] == self.chi[i1][j1][b]:
                # The lines don't even intersect.
                continue

            # Orient i,j,a counter-clock-wise.
            if self.chi[i1][j1][a] == -1:
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
                if self.chi[c][d][j] == 1:
                    # j lies above c,d, (i is c or d) therefore
                    # i-j intersects a,b, above c,d.
                    above_below[index] = 1
                else:
                    above_below[index] = -1
            elif j in (c,d):
                if self.chi[c][d][i] == 1:
                    # As above.
                    above_below[index] = 1
                else:
                    above_below[index] = -1
            else:
                # In all remaining cases i,j must both lie above or both below c,d,
                # otherwise the location of the section i,j with respect to the intersection point (a,b,c,d)
                # would be up to a choice and would therefore be needed to be in the dictionary.
                assert self.chi[c][d][i] == self.chi[c][d][j]
                if self.chi[c][d][i] == 1:
                    above_below[index] = 1
                else:
                    above_below[index] = -1
        return tuple(above_below), index_cd

    cdef Bitset_wrapper partitions_one_bitset(self):
        r"""
        Return a bitset that has each color index set to 1, for which there is Tverberg partition of type 3, 3, 3, 1.

        The trailing bits are also set to 1.
        """
        cdef Bitset_wrapper foo
        if self._partitions_one_bitset is None:
            # Set ``foo`` to a bitset with all trailing bits set.
            cap = ((n_colors() - 1) // 256 + 1) * 256
            foo = Bitset_wrapper(cap, range(n_colors(), cap))

            for part in self._partitions_one():
                # Set the indices the 1 that correspond to color partitions,
                # for which ``part`` is rainbow.
                foo.union_assign(foo, colors_for_partition(part))

            self._partitions_one_bitset = foo

        return self._partitions_one_bitset

    def _partitions_one(self):
        r"""
        Obtain all Tverberg partitions of type (3, 3, 3, 1) determined by this chirotope.
        """
        for x in range(10):
            others = tuple(i for i in range(10) if i != x)
            a = others[0]
            for b, c in combinations(others[1:], 2):
                # We apply Lemma 2.7 to check whether x in conv(a, b, c).
                # TODO: Relabel according to \ref{Lem:Restriction}.
                if self.chi[a][b][x] == self.chi[b][c][x] == self.chi[c][a][x]:
                    others1 = tuple(i for i in range(10) if not i in (x,a,b,c))
                    d = others1[0]
                    for e, f in combinations(others1[1:], 2):
                        # We apply Lemma 2.7 to check whether x in conv(d, e, f).
                        # TODO: Relabel according to \ref{Lem:Restriction}.
                        if self.chi[d][e][x] == self.chi[e][f][x] == self.chi[f][d][x]:
                            g, h, i = tuple(i for i in range(10) if not i in (x,a,b,c,d,e,f))
                            # We apply Lemma 2.7 to check whether x in conv(g, h, i).
                            # TODO: Relabel according to \ref{Lem:Restriction}.
                            if self.chi[g][h][x] == self.chi[h][i][x] == self.chi[i][g][x]:
                                yield ((x,), (a,b,c), (d,e,f), (g,h,i))

    def partitions_two(self, dict dic, int a, int b, int c, int d):
        """
        Yield all Tverberg partitions given by a certain dictionary.
        """
        regions, regions_op, others = self.regions(a,b,c,d)
        intersection = ((a,b), (c,d))
        cdef int x,y,w,e,f,g

        def is_triangle(x,y,z):
            """
            Check if x,y,z is a triangle containing the intersection of ab and cd.
            """
            for i,j in ((0,3),(1,2)):
                if (any(w in regions[i] for w in (x,y,z)) and any(w in regions[j] for w in (x,y,z))):
                    k,l = (0,3) if i == 1 else (1,2)
                    # Now the opposite regions i,j contain one point of x,y,z.

                    if any(w in regions[k] or w in regions[l] for w in (x,y,z)):
                        # the vertices are in three regions, so it is all up to one line
                        v1,v2 = tuple(w for w in (x,y,z) if not (w in regions[k] or w in regions[l]))
                        v3 = tuple(w for w in (x,y,z) if not w in (v1,v2))[0]
                        if (v1,v2,intersection) in dic:
                            return dic[(v1,v2,intersection)] == self.chi[v1][v2][v3]
                        return dic[(v2,v1,intersection)] == self.chi[v2][v1][v3]
                    else:
                        # the vertices are only in two sections
                        one = tuple(w for w in (x,y,z) if w in regions[i])
                        two = tuple(w for w in (x,y,z) if w in regions[j])
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

    cdef Bitset_wrapper rainbow_partitions(self, int y, int z):
        r"""
        Return a bitset for the ``z``-th choice for the intersection point ``y`.

        The bitset indicates which color partition
        does not have a rainbow Tverberg partition
        of type (3, 3, 3, 1) or of type (3, 3, 2, 2) in ``y``.

        It is set to 0 if this color partition has a rainbow Tverberg
        partition of type (3, 3, 3, 1) or of type (3, 3, 2, 2) in ``y``.
        """
        cdef dict chi
        cdef int a, b, c, d
        chi, a, b, c, d = self.obtain_dic(y, z)
        cdef Bitset_wrapper foo = self.partitions_one_bitset().copy()

        for part in self.partitions_two(chi, a, b, c, d):
            foo.union_assign(foo, colors_for_partition(part))

        # Note that the trailing bits are set to 1 at this place.
        foo.flip_inplace()
        return foo

    def poss_graph(self):
        """
        Used to obtain the sample graphs
        in ``https://github.com/kliem/PyKPartiteKClique``.
        """
        from sage.graphs.graph import Graph
        G = Graph()
        ls = self.valid_intersection_points()
        n = len(ls)
        for i,j in combinations(range(n), 2):
            x = ls[i]
            y = ls[j]
            for k,a in enumerate(self.extensions(*x[0], *x[1])):
                for l,b in enumerate(self.extensions(*y[0], *y[1])):
                    if self.check_for_consistency(i,k,j,l):
                        G.add_edge((i,k), (j,l))

        # So far we have added all the vertices corresponding to valid intersection points
        # and orientations.
        # For each such vertex, we add an edge to those color partions, for which it does
        # not induce a rainbow partition.
        V = G.vertices()
        cdef int num_colors = n_colors()
        for v in V:
            rainbows = self.rainbow_partitions(v[0], v[1])
            l = rainbows.first(0)
            while l < num_colors:
                G.add_edge(v, (n, l))
                l = rainbows.first(l+1)
        return G


def poss_color_finder(Chirotope O):
    r"""
    Determine all colors for which this Chirotope has a solution:

    There exists some setup with the intersections so that the colors don't have a Tverberg Partition.

    This is done by finding a k-clique in a k-partite graph.

    The last part contains one node for each color.

    All the other parts represent an intersection point each.
    The oriented matroid does not decode where this intersection point lies with respect to other lines (between two points each).
    However, there is only a finite number of possibilities to place such an intersection point.

    An edge is placed wherever such a choice can be made consistenly.
    (Or rather no edge is place, if there is a contradiction.)
    Further each choice is connected to all colors, for which it does not yield a colorful Tverberg partition.

    A maximal clique in this graph would correspond to a consistent choice, for which there is no colorful Tverberg partition.
    If those choices could be realized, then this would give a counter example.

    However, at least for 10 points we can show that there never exists such a maximal clique, which means
    that already on the level of oriented matroids we can show that there exists always at least on colorful Tverberg partition.
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef Bitset_wrapper somea  # type awareness of cython
    cdef int i, j

    cdef int k = O.n_valid_intersection_points() + 1
    cdef int* first_per_part = <int*> mem.allocarray(k, sizeof(int))
    cdef int counter = 0
    first_per_part[0] = 0
    for i in range(k-1):
        counter += O._n_extensions(i)
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
    cdef int n_color_bits

    for i in range(k-1):
        offset_i = first_per_part[i]
        for ind_i in range(O._n_extensions(i)):
            somea = O.rainbow_partitions(i, ind_i)
            ind_color = somea.first(0)
            n_color_bits = somea.n_bits
            while ind_color < n_color_bits:
                incidences[offset_i + ind_i][offset_colors + ind_color] = True
                incidences[offset_colors + ind_color][offset_i + ind_i] = True
                ind_color = somea.first(ind_color + 1)

    cdef int l

    # All connections
    # For each constellation we store, which other constellations are still possible, if this is picked.
    cdef int ind_j, offset_j
    for i in range(k-1):
        offset_i = first_per_part[i]
        for ind_i in range(O._n_extensions(i)):
            for j in range(i):
                offset_j = first_per_part[j]
                for ind_j in range(O._n_extensions(j)):
                    # There is an arc between those vertices if and only if i,ind_i and j,ind_j are consistent.
                    if O.check_for_consistency(i, ind_i, j, ind_j):
                        incidences[offset_i + ind_i][offset_j+ ind_j] = True
                        incidences[offset_j + ind_j][offset_i+ ind_i] = True

    cdef KPartiteKClique* K = new KPartiteKClique(incidences, n, first_per_part, k)
    try:
        foo = K.next()
        if foo:
            # There is a counter example.
            raise ValueError
    finally:
        del K

    return []


def OrderTypesIterator(filename=None):
    r"""
    Iterate through all order types in general position for 10 points.

    See http://www.ist.tugraz.at/staff/aichholzer/research/rp/triangulations/ordertypes/
    """
    if filename is None:
        filename = "/srv/public/kliem/OrderTypes/otypes10.b16"
    with open(filename, "rb") as f:
        data = f.read(40)
        while len(data) == 40:
            data = f.read(40)
            yield tuple((int.from_bytes(data[4*i : 4*i + 2], "little"), int.from_bytes(data[4*i + 2:4*i + 4], "little")) for i in range(10))


def OrderType(n):
    r"""
    Return the n-th order type.

    See ``OrderTypesIterator``.
    """
    filename = "/srv/public/kliem/OrderTypes/otypes10.b16"
    with open(filename, "rb") as f:
        f.seek(40*n)
        data = f.read(40)
        if len(data) != 40:
            raise IndexError
        return tuple((int.from_bytes(data[4*i:4*i+2], "little"), int.from_bytes(data[4*i+2:4*i+4], "little")) for i in range(10))


def check_all_colors(start=0, end=2**20):
    r"""
    Check for each order set the possible counter example color indices, if any.
    """
    return _check_all_colors(start, end, OrderTypesIterator())


def check_all_colors_pseudo(start=0, end=2**20):
    r"""
    Check for each pseoduo order set the possible counter example color indices, if any.
    """
    from pseudo_order_types.pseudo_order_types import pseudo_order_type_iterator
    return _check_all_colors(start, end, pseudo_order_type_iterator(10, "/srv/public/kliem/OrderSets/10"))


def _check_all_colors(start, end, iterator):
    """
    See above.
    """
    for i, O1 in enumerate(iterator):
        if i < start:
            continue
        if i >= end:
            break
        O = Chirotope(O1)
        if i % 1000 == 0:
            print("currently doing: ", i)
        try:
            _ = poss_color_finder(O)
        except ValueError:
            # There appears to be a counter example.
            raise ValueError("found counter example at index {}".format(i))
