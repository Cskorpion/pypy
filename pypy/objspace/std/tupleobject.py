"""The builtin tuple implementation"""

import sys

from pypy.interpreter.baseobjspace import W_Root
from pypy.interpreter.error import OperationError, oefmt
from pypy.interpreter.gateway import (
    WrappedDefault, interp2app, interpindirect2app, unwrap_spec)
from pypy.interpreter.typedef import TypeDef
from pypy.objspace.std.objectobject import ObjectObjectDocstrings
from pypy.objspace.std.sliceobject import (W_SliceObject, unwrap_start_stop,
    normalize_simple_slice)
from pypy.objspace.std.util import negate, IDTAG_SPECIAL, IDTAG_SHIFT
from rpython.rlib import jit
from rpython.rlib.debug import make_sure_not_resized
from rpython.rlib.rarithmetic import intmask


UNROLL_CUTOFF = 10


def _unroll_condition_cmp(self, space, other):
    return self._unroll_condition() or other._unroll_condition()


def get_printable_location(tp):
    return "tuple.contains [%s]" % (tp.getname(tp.space), )

contains_driver = jit.JitDriver(greens = ['tp'], reds = 'auto',
                             name = 'tuple.contains',
                             get_printable_location=get_printable_location)

def get_printable_location(w_type):
    return "tuple.hash [%s]" % (w_type.getname(w_type.space), )

hash_driver = jit.JitDriver(
    name='tuple.hash',
    greens=['w_type'],
    reds='auto',
    get_printable_location=get_printable_location
    )

class W_AbstractTupleObject(W_Root):
    __slots__ = ()

    def is_w(self, space, w_other):
        if not isinstance(w_other, W_AbstractTupleObject):
            return False
        if self is w_other:
            return True
        if self.user_overridden_class or w_other.user_overridden_class:
            return False
        # empty tuples are unique-ified
        return 0 == w_other.length() == self.length()

    def immutable_unique_id(self, space):
        if self.user_overridden_class or self.length() > 0:
            return None
        # empty tuple: base value 258
        uid = (258 << IDTAG_SHIFT) | IDTAG_SPECIAL
        return space.newint(uid)

    def __repr__(self):
        """representation for debugging purposes"""
        reprlist = [repr(w_item) for w_item in self.tolist()]
        return "%s(%s)" % (self.__class__.__name__, ', '.join(reprlist))

    def unwrap(self, space):
        items = [space.unwrap(w_item) for w_item in self.tolist()]
        return tuple(items)

    def tolist(self):
        """Returns the items, as a fixed-size list."""
        raise NotImplementedError

    def getitems_copy(self):
        """Returns a copy of the items, as a resizable list."""
        raise NotImplementedError

    def length(self):
        raise NotImplementedError

    def getitem(self, space, item):
        raise NotImplementedError

    def descr_len(self, space):
        result = self.length()
        return space.newint(result)

    def descr_iter(self, space):
        from pypy.objspace.std import iterobject
        return iterobject.W_FastTupleIterObject(self, self.tolist())

    @staticmethod
    def descr_new(space, w_tupletype, w_sequence=None):
        if w_sequence is None:
            tuple_w = []
        elif (space.is_w(w_tupletype, space.w_tuple) and
              space.is_w(space.type(w_sequence), space.w_tuple)):
            return w_sequence
        else:
            tuple_w = space.fixedview(w_sequence)
        w_obj = space.allocate_instance(W_TupleObject, w_tupletype)
        W_TupleObject.__init__(w_obj, tuple_w)
        return w_obj

    def descr_repr(self, space):
        items = self.tolist()
        if len(items) == 1:
            return space.newtext("(" + space.text_w(space.repr(items[0])) + ",)")
        tmp = ", ".join([space.text_w(space.repr(item)) for item in items])
        return space.newtext("(" + tmp + ")")

    def descr_hash(self, space):
        raise NotImplementedError

    def descr_eq(self, space, w_other):
        raise NotImplementedError

    def descr_ne(self, space, w_other):
        raise NotImplementedError

    def _make_tuple_comparison(name):
        import operator
        op = getattr(operator, name)

        def compare_tuples(self, space, w_other):
            if not isinstance(w_other, W_AbstractTupleObject):
                return space.w_NotImplemented
            return _compare_tuples(self, space, w_other)

        @jit.look_inside_iff(_unroll_condition_cmp)
        def _compare_tuples(self, space, w_other):
            items1 = self.tolist()
            items2 = w_other.tolist()
            ncmp = min(len(items1), len(items2))
            # Search for the first index where items are different
            for p in range(ncmp):
                if not space.eq_w(items1[p], items2[p]):
                    return getattr(space, name)(items1[p], items2[p])
            # No more items to compare -- compare sizes
            return space.newbool(op(len(items1), len(items2)))

        compare_tuples.__name__ = 'descr_' + name
        return compare_tuples

    descr_lt = _make_tuple_comparison('lt')
    descr_le = _make_tuple_comparison('le')
    descr_gt = _make_tuple_comparison('gt')
    descr_ge = _make_tuple_comparison('ge')

    def descr_contains(self, space, w_obj):
        if self._unroll_condition():
            return self._descr_contains_unroll_safe(space, w_obj)
        else:
            return self._descr_contains_jmp(space, w_obj)

    @jit.unroll_safe
    def _descr_contains_unroll_safe(self, space, w_obj):
        for w_item in self.tolist():
            if space.eq_w(w_obj, w_item):
                return space.w_True
        return space.w_False

    def _descr_contains_jmp(self, space, w_obj):
        tp = space.type(w_obj)
        list_w = self.tolist()
        i = 0
        while i < len(list_w):
            contains_driver.jit_merge_point(tp=tp)
            w_item = list_w[i]
            if space.eq_w(w_obj, w_item):
                return space.w_True
            i += 1
        return space.w_False

    def descr_add(self, space, w_other):
        if not isinstance(w_other, W_AbstractTupleObject):
            return space.w_NotImplemented
        items1 = self.tolist()
        items2 = w_other.tolist()
        return space.newtuple(items1 + items2)

    def descr_mul(self, space, w_times):
        try:
            times = space.getindex_w(w_times, space.w_OverflowError)
        except OperationError as e:
            if e.match(space, space.w_TypeError):
                return space.w_NotImplemented
            raise
        if times == 1 and space.type(self) == space.w_tuple:
            return self
        items = self.tolist()
        return space.newtuple(items * times)

    def descr_getitem(self, space, w_index):
        if isinstance(w_index, W_SliceObject):
            return self._getslice(space, w_index)
        index = space.getindex_w(w_index, space.w_IndexError, "tuple index")
        return self.getitem(space, index)

    def _getslice(self, space, w_index):
        items = self.tolist()
        length = len(items)
        start, stop, step, slicelength = w_index.indices4(space, length)
        if slicelength == 0:
            subitems = []
        elif step == 1:
            assert 0 <= start <= stop
            subitems = items[start:stop]
        else:
            subitems = self._getslice_advanced(items, start, step, slicelength)
        return space.newtuple(subitems)

    @staticmethod
    def _getslice_advanced(items, start, step, slicelength):
        assert slicelength >= 0
        subitems = [None] * slicelength
        for i in range(slicelength):
            subitems[i] = items[start]
            start += step
        return subitems

    def descr_getslice(self, space, w_start, w_stop):
        length = self.length()
        start, stop = normalize_simple_slice(space, length, w_start, w_stop)
        return space.newtuple(self.tolist()[start:stop])

    def descr_getnewargs(self, space):
        return space.newtuple([space.newtuple(self.tolist())])

    @jit.look_inside_iff(lambda self, _1, _2: self._unroll_condition())
    def descr_count(self, space, w_obj):
        """count(obj) -> number of times obj appears in the tuple"""
        count = 0
        for w_item in self.tolist():
            if space.eq_w(w_item, w_obj):
                count += 1
        return space.newint(count)

    @unwrap_spec(w_start=WrappedDefault(0), w_stop=WrappedDefault(sys.maxint))
    @jit.look_inside_iff(lambda self, _1, _2, _3, _4: self._unroll_condition())
    def descr_index(self, space, w_obj, w_start, w_stop):
        """index(obj, [start, [stop]]) -> first index that obj appears in the
        tuple
        """
        length = self.length()
        start, stop = unwrap_start_stop(space, length, w_start, w_stop)
        for i in range(start, min(stop, length)):
            w_item = self.tolist()[i]
            if space.eq_w(w_item, w_obj):
                return space.newint(i)
        raise oefmt(space.w_ValueError, "tuple.index(x): x not in tuple")

    def _unroll_condition(self):
        raise NotImplementedError("abstract base class")

class TupleDocstrings:

    def __eq__():
        """x.__eq__(y) <==> x==y'"""

    def __ne__():
        """x.__ne__(y) <==> x!=y"""

    def __lt__():
        """x.__lt__(y) <==> x<y"""

    def __le__():
        """x.__le__(y) <==> x<=y"""

    def __gt__():
        """x.__gt__(y) <==> x>y"""

    def __ge__():
        """x.__ge__(y) <==> x>=y"""

    def __len__():
        """x.__len__() <==> len(x)"""

    def __iter__():
        """x.__iter__() <==> iter(x)"""

    def __contains__():
        """x.__contains__(y) <==> y in x"""

    def __add__():
        """x.__add__(y) <==> x+y"""

    def __mul__():
        """x.__mul__(n) <==> x*n"""

    def __rmul__():
        """x.__rmul__(n) <==> n*x"""

    def __getitem__():
        """x.__getitem__(y) <==> x[y]"""

    def __getslice__():
        """x.__getslice__(i, j) <==> x[i:j]

        Use of negative indices is not supported."""

    def __getnewargs__():
        """""" # empty

    def count():
        """T.count(value) -> integer -- return number of occurrences of value"""

    def index():
        """"T.index(value, [start, [stop]]) -> integer -- return first index of value.
        Raises ValueError if the value is not present."""

W_AbstractTupleObject.typedef = TypeDef(
    "tuple",
    __doc__ = """tuple() -> an empty tuple
tuple(sequence) -> tuple initialized from sequence's items

If the argument is a tuple, the return value is the same object.""",
    __new__ = interp2app(W_AbstractTupleObject.descr_new,
                         doc=ObjectObjectDocstrings.__new__.__doc__),
    __repr__ = interp2app(W_AbstractTupleObject.descr_repr,
                          doc=ObjectObjectDocstrings.__repr__.__doc__),
    __hash__ = interpindirect2app(W_AbstractTupleObject.descr_hash),

    __eq__ = interpindirect2app(W_AbstractTupleObject.descr_eq),
    __ne__ = interpindirect2app(W_AbstractTupleObject.descr_ne),
    __lt__ = interp2app(W_AbstractTupleObject.descr_lt,
                        doc=TupleDocstrings.__lt__.__doc__),
    __le__ = interp2app(W_AbstractTupleObject.descr_le,
                        doc=TupleDocstrings.__le__.__doc__),
    __gt__ = interp2app(W_AbstractTupleObject.descr_gt,
                        doc=TupleDocstrings.__gt__.__doc__),
    __ge__ = interp2app(W_AbstractTupleObject.descr_ge,
                        doc=TupleDocstrings.__ge__.__doc__),

    __len__ = interp2app(W_AbstractTupleObject.descr_len,
                         doc=TupleDocstrings.__len__.__doc__),
    __iter__ = interp2app(W_AbstractTupleObject.descr_iter,
                          doc=TupleDocstrings.__iter__.__doc__),
    __contains__ = interp2app(W_AbstractTupleObject.descr_contains,
                              doc=TupleDocstrings.__contains__.__doc__),

    __add__ = interp2app(W_AbstractTupleObject.descr_add,
                         doc=TupleDocstrings.__add__.__doc__),
    __mul__ = interp2app(W_AbstractTupleObject.descr_mul,
                         doc=TupleDocstrings.__mul__.__doc__),
    __rmul__ = interp2app(W_AbstractTupleObject.descr_mul,
                          doc=TupleDocstrings.__rmul__.__doc__),

    __getitem__ = interp2app(W_AbstractTupleObject.descr_getitem,
                             doc=TupleDocstrings.__getitem__.__doc__),
    __getslice__ = interp2app(W_AbstractTupleObject.descr_getslice,
                              doc=TupleDocstrings.__getslice__.__doc__),

    __getnewargs__ = interp2app(W_AbstractTupleObject.descr_getnewargs,
                                doc=TupleDocstrings.__getnewargs__.__doc__),
    count = interp2app(W_AbstractTupleObject.descr_count,
                       doc=TupleDocstrings.count.__doc__),
    index = interp2app(W_AbstractTupleObject.descr_index,
                       doc=TupleDocstrings.index.__doc__)
)
W_AbstractTupleObject.typedef.flag_sequence_bug_compat = True


class W_TupleObject(W_AbstractTupleObject):
    _immutable_fields_ = ['wrappeditems[*]']

    def __init__(self, wrappeditems):
        make_sure_not_resized(wrappeditems)
        self.wrappeditems = wrappeditems

    def tolist(self):
        return self.wrappeditems

    def getitems_copy(self):
        return self.wrappeditems[:]  # returns a resizable list

    def length(self):
        return len(self.wrappeditems)

    def descr_hash(self, space):
        if self._unroll_condition():
            res = self._descr_hash_unroll(space)
        else:
            res = self._descr_hash_jitdriver(space)
        return space.newint(res)

    @jit.unroll_safe
    def _descr_hash_unroll(self, space):
        mult = 1000003
        x = 0x345678
        z = len(self.wrappeditems)
        for w_item in self.wrappeditems:
            y = space.hash_w(w_item)
            x = (x ^ y) * mult
            z -= 1
            mult += 82520 + z + z
        x += 97531
        return intmask(x)

    def _descr_hash_jitdriver(self, space):
        mult = 1000003
        x = 0x345678
        z = len(self.wrappeditems)
        w_type = space.type(self.wrappeditems[0])
        wrappeditems = self.wrappeditems
        i = 0
        while i < len(wrappeditems):
            hash_driver.jit_merge_point(w_type=w_type)
            w_item = wrappeditems[i]
            y = space.hash_w(w_item)
            x = (x ^ y) * mult
            z -= 1
            mult += 82520 + z + z
            i += 1
        x += 97531
        return intmask(x)

    def descr_eq(self, space, w_other):
        if not isinstance(w_other, W_AbstractTupleObject):
            return space.w_NotImplemented
        return self._descr_eq(space, w_other)

    @jit.look_inside_iff(_unroll_condition_cmp)
    def _descr_eq(self, space, w_other):
        items1 = self.wrappeditems
        items2 = w_other.tolist()
        lgt1 = len(items1)
        lgt2 = len(items2)
        if lgt1 != lgt2:
            return space.w_False
        # XXX do we need a jit driver?
        for i in range(lgt1):
            item1 = items1[i]
            item2 = items2[i]
            if not space.eq_w(item1, item2):
                return space.w_False
        return space.w_True

    descr_ne = negate(descr_eq)

    def getitem(self, space, index):
        try:
            return self.wrappeditems[index]
        except IndexError:
            raise oefmt(space.w_IndexError, "tuple index out of range")

    def _unroll_condition(self):
        return jit.loop_unrolling_heuristic(
                self.wrappeditems, self.length(), UNROLL_CUTOFF)


def wraptuple(space, list_w):
    if space.config.objspace.std.withspecialisedtuple:
        from specialisedtupleobject import makespecialisedtuple, NotSpecialised
        try:
            return makespecialisedtuple(space, list_w)
        except NotSpecialised:
            pass
    return W_TupleObject(list_w)

def wraptuple2(space, w_a, w_b):
    if space.config.objspace.std.withspecialisedtuple:
        from specialisedtupleobject import makespecialisedtuple2, NotSpecialised
        try:
            return makespecialisedtuple2(space, w_a, w_b)
        except NotSpecialised:
            pass
    return W_TupleObject([w_a, w_b])

