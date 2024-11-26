import pytest
import sys
from rpython.tool.udir import udir

class AppTestVMProf(object):
    spaceconfig = {'usemodules': ['_vmprof', 'struct']}

    def setup_class(cls):
        cls.w_tmpfilename = cls.space.wrap(str(udir.join('test__vmprof.1')))
        cls.w_tmpfilename2 = cls.space.wrap(str(udir.join('test__vmprof.2')))
        cls.w_plain = cls.space.wrap(not cls.runappdirect and
            '__pypy__' not in sys.builtin_module_names)

    def test_import_vmprof(self):
        tmpfile = open(self.tmpfilename, 'wb')
        tmpfileno = tmpfile.fileno()
        tmpfile2 = open(self.tmpfilename2, 'wb')
        tmpfileno2 = tmpfile2.fileno()

        import struct, sys, gc

        WORD = struct.calcsize('l')

        def count(s):
            i = 0
            count = 0
            i += 5 * WORD # header
            assert s[i    ] == '\x05'    # MARKER_HEADER
            assert s[i + 1] == '\x00'    # 0
            assert s[i + 2] == '\x07'    # VERSION_SAMPLE_TIMEOFFSET
            assert s[i + 3] == '\x08'    # PROFILE_RPYTHON
            assert s[i + 4] == chr(4)    # len('pypy')
            assert s[i + 5: i + 9] == 'pypy'
            i += 9
            while i < len(s):
                if s[i] == '\x03':
                    break
                elif s[i] == '\x01':
                    i += 1
                    _, size = struct.unpack("ll", s[i:i + 2 * WORD])
                    i += 2 * WORD + size * struct.calcsize("P")
                    i += WORD    # thread id
                elif s[i] == '\x02':
                    i += 1
                    _, size = struct.unpack("ll", s[i:i + 2 * WORD])
                    count += 1
                    i += 2 * WORD + size
                elif s[i] == '\x06':
                    print(s[i:i+24])
                    i += 1+8+8+8
                elif s[i] == '\x07':
                    i += 1
                    # skip string
                    size, = struct.unpack("l", s[i:i + WORD])
                    i += WORD+size
                    # skip string
                    size, = struct.unpack("l", s[i:i + WORD])
                    i += WORD+size
                else:
                    raise AssertionError(ord(s[i]))
            return count

        import _vmprof
        gc.collect()  # try to make the weakref list deterministic
        gc.collect()  # by freeing all dead code objects
        _vmprof.enable(tmpfileno, 0.01, 0, 0, 0, 0)
        _vmprof.disable()
        s = open(self.tmpfilename, 'rb').read()
        no_of_codes = count(s)
        assert no_of_codes > 10
        d = {}

        exec """def foo():
            pass
        """ in d

        gc.collect()
        gc.collect()
        _vmprof.enable(tmpfileno2, 0.01, 0, 0, 0, 0)

        exec """def foo2():
            pass
        """ in d

        _vmprof.disable()
        s = open(self.tmpfilename2, 'rb').read()
        no_of_codes2 = count(s)
        assert "py:foo:" in s
        assert "py:foo2:" in s
        assert no_of_codes2 >= no_of_codes + 2 # some extra codes from tests

    def test_enable_ovf(self):
        import _vmprof

        # Disabled, bacause allocation based sampling allows interval == 0
        #raises(_vmprof.VMProfError, _vmprof.enable, 2, 0, 0, 0, 0, 0) 
        
        raises(_vmprof.VMProfError, _vmprof.enable, 2, -2.5, 0, 0, 0, 0)
        raises(_vmprof.VMProfError, _vmprof.enable, 2, 1e300, 0, 0, 0, 0)
        raises(_vmprof.VMProfError, _vmprof.enable, 2, 1e300 * 1e300, 0, 0, 0, 0)
        NaN = (1e300*1e300) / (1e300*1e300)
        raises(_vmprof.VMProfError, _vmprof.enable, 2, NaN, 0, 0, 0, 0)

    def test_is_enabled(self):
        import _vmprof
        tmpfile = open(self.tmpfilename, 'wb')
        assert _vmprof.is_enabled() is False
        _vmprof.enable(tmpfile.fileno(), 0.01, 0, 0, 0, 0)
        assert _vmprof.is_enabled() is True
        _vmprof.disable()
        assert _vmprof.is_enabled() is False

    @pytest.mark.xfail(sys.platform.startswith('freebsd'), reason = "not implemented")
    def test_get_profile_path(self):
        import _vmprof
        with open(self.tmpfilename, "wb") as tmpfile:
            assert _vmprof.get_profile_path() is None
            _vmprof.enable(tmpfile.fileno(), 0.01, 0, 0, 0, 0)
            path = _vmprof.get_profile_path()
            _vmprof.disable()

        if path != tmpfile.name:
            with open(path, "rb") as fd1:
                with open(self.tmpfilename, "rb") as fd2:
                    assert fd1.read() == fd2.read()

        assert _vmprof.get_profile_path() is None

    def test_stop_sampling(self):
        if not self.plain:
            skip("unreliable test except on CPython without -A")
        import os
        import _vmprof
        tmpfile = open(self.tmpfilename, 'wb')
        native = 1
        def f():
            import sys
            import math
            j = sys.maxsize
            for i in range(500):
                j = math.sqrt(j)
        _vmprof.enable(tmpfile.fileno(), 0.01, 0, native, 0, 0)
        # get_vmprof_stack() always returns 0 here!
        # see vmprof_common.c and assume RPYTHON_LL2CTYPES is defined!
        f()
        fileno = _vmprof.stop_sampling()
        pos = os.lseek(fileno, 0, os.SEEK_CUR)
        f()
        pos2 = os.lseek(fileno, 0, os.SEEK_CUR)
        assert pos == pos2
        _vmprof.start_sampling()
        f()
        fileno = _vmprof.stop_sampling()
        pos3 = os.lseek(fileno, 0, os.SEEK_CUR)
        assert pos3 > pos
        _vmprof.disable()

    def test_enable_allocation_triggered(self):
        # Only works on Unix
        import _vmprof

        tmpfile = open(self.tmpfilename, 'wb')
        fd = tmpfile.fileno()
        MARKER_GC_STACKTRACE = b'\x09'

        #_vmprof.enable_allocation_triggered(fd, 0)# prepare everything but dont sample in gc
        _vmprof.enable(fd, 0, 0, 0, 0, 0) # workarround for constant folding error

        _vmprof.sample_stack_now()# manually trigger some samples
        _vmprof.sample_stack_now()
        _vmprof.sample_stack_now()
        _vmprof.sample_stack_now()

        _vmprof.disable()

        s = open(self.tmpfilename, "rb").read()

        assert s.count(MARKER_GC_STACKTRACE) == 4 or s.count(MARKER_GC_STACKTRACE) == 5 # sometimes there is a 0x09 in there that is not a marker

    def test_allocation_and_time_sampling(self):
        # Only works on Unix
        import _vmprof

        tmpfile = open(self.tmpfilename, 'wb')
        fd = tmpfile.fileno()
        MARKER_GC_STACKTRACE = b'\x09'
        MARKER_STACKTRACE = b'\x01'

        #_vmprof.enable_allocation_triggered(fd, 0, 0.0001, 0)# prepare everything but dont sample in gc
        _vmprof.enable(fd, 0.0001, 0, 0, 0, 0) # workarround for constant folding error
        _vmprof.start_sampling()

        counter = 0

        for i in range(100):
            # manually trigger some (gc) samples
            _vmprof.sample_stack_now()
            
            for _ in range(10000):
                counter += 1
            if counter == -1: 
                break

        _vmprof.disable()

        s = open(self.tmpfilename, "rb").read()

        assert s.count(MARKER_GC_STACKTRACE) > 0
        assert s.count(MARKER_STACKTRACE) > 0
