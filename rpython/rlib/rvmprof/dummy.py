from rpython.rlib.objectmodel import specialize
from rpython.rlib.rvmprof import cintf

class DummyVMProf(object):
    is_enabled = False
    cintf = cintf.CInterface({"vmprof_say_hi": lambda: None,
                              "vmprof_sample_stack_now_gc_triggered": lambda: None})

    def __init__(self):
        self._unique_id = 0

    def register_code_object_class(self, CodeClass, full_name_func):
        CodeClass._vmprof_unique_id = self._unique_id
        self._unique_id += 1

    @specialize.argtype(1)
    def register_code(self, code, full_name_func):
        pass

    def enable(self, fileno, interval, memory=0, native=0, real_time=0):
        pass

    def enable_allocation_triggered(self, fileno, sample_n_bytes=1024, interval=0.0, native=0):
        pass

    def sample_stack_now(self):
        pass

    def disable(self):
        pass

    def start_sampling(self):
        pass

    def stop_sampling(self):
        return -1
    
    def vmprof_report_minor_gc_objs(self, time, array, array_size):
        pass
