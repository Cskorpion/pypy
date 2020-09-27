## DO NOT EDIT THIS FILE, IT IS AUTOGENERATED
from pypy.module._hpy_universal import llapi, handles
from pypy.module._hpy_universal.interp_slot import W_SlotWrapper
from pypy.module._hpy_universal.state import State

class W_SlotWrapper_unaryfunc(W_SlotWrapper):
    def call(self, space, __args__):
        self.check_args(space, __args__, 1)
        func = llapi.cts.cast("HPyFunc_unaryfunc", self.cfuncptr)
        ctx = space.fromcache(State).ctx
        w0 = __args__.arguments_w[0]
        with handles.using(space, w0) as c0:
            c_result = func(ctx, c0)
            return handles.consume(space, c_result)

class W_SlotWrapper_binaryfunc(W_SlotWrapper):
    def call(self, space, __args__):
        self.check_args(space, __args__, 2)
        func = llapi.cts.cast("HPyFunc_binaryfunc", self.cfuncptr)
        ctx = space.fromcache(State).ctx
        w0 = __args__.arguments_w[0]
        w1 = __args__.arguments_w[1]
        with handles.using(space, w0) as c0, handles.using(space, w1) as c1:
            c_result = func(ctx, c0, c1)
            return handles.consume(space, c_result)

class W_SlotWrapper_ssizeargfunc(W_SlotWrapper):
    def call(self, space, __args__):
        self.check_args(space, __args__, 2)
        func = llapi.cts.cast("HPyFunc_ssizeargfunc", self.cfuncptr)
        ctx = space.fromcache(State).ctx
        w0 = __args__.arguments_w[0]
        w1 = __args__.arguments_w[1]
        c1 = space.int_w(space.index(w1))
        with handles.using(space, w0) as c0:
            c_result = func(ctx, c0, c1)
            return handles.consume(space, c_result)

class W_SlotWrapper_reprfunc(W_SlotWrapper):
    def call(self, space, __args__):
        self.check_args(space, __args__, 1)
        func = llapi.cts.cast("HPyFunc_reprfunc", self.cfuncptr)
        ctx = space.fromcache(State).ctx
        w0 = __args__.arguments_w[0]
        with handles.using(space, w0) as c0:
            c_result = func(ctx, c0)
            return handles.consume(space, c_result)

