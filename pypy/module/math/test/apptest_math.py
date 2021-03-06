import math

def test_product():
    assert math.prod([1, 2, 3]) == 6
    assert math.prod([1, 2, 3], start=0.5) == 3.0
    assert math.prod([]) == 1.0
    assert math.prod([], start=5) == 5

def test_julians_weird_test():
    class A:
        def __mul__(self, other):
                return 12
        def __imul__(self, other):
                return 13

    # check that the implementation doesn't use *=
    assert math.prod([1, 2], start=A())

def test_more_weird():
    start = [4]
    assert math.prod([2], start=start) == [4, 4]
    assert start == [4]
    start =  object()
    assert math.prod([], start=start) is start