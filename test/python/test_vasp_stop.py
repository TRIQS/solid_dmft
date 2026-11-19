#!/usr/bin/env python
# vasp_manager.stop at the end of a CSC run, with forked processes in place of VASP:
# one that behaves like VASP in the CSC loop (waits for vasp.lock, reads LABORT from
# STOPCAR and exits), and one that hangs and has to be killed after the timeout.

import os
import tempfile
import time
import unittest

from solid_dmft.dft_managers import vasp_manager


def _fork(child):
    pid = os.fork()
    if pid == 0:
        try:
            child()
        finally:
            os._exit(0)
    return pid


def _fake_vasp():
    while True:
        if os.path.isfile('vasp.lock') and os.path.isfile('STOPCAR'):
            with open('STOPCAR') as stopcar:
                if 'LABORT = .TRUE.' in stopcar.read():
                    os.remove('vasp.lock')
                    return
        time.sleep(0.05)


def _hanging_process():
    time.sleep(60)


class TestVaspStop(unittest.TestCase):

    def setUp(self):
        self.cwd = os.getcwd()
        self.tmp = tempfile.TemporaryDirectory()
        os.chdir(self.tmp.name)

    def tearDown(self):
        os.chdir(self.cwd)
        self.tmp.cleanup()

    def test_vasp_exits_by_itself(self):
        pid = _fork(_fake_vasp)
        self.assertTrue(vasp_manager.stop(pid, timeout=20))
        self.assertFalse(os.path.isfile('vasp.lock'))

    def test_hanging_vasp_is_killed(self):
        pid = _fork(_hanging_process)
        start = time.time()
        self.assertFalse(vasp_manager.stop(pid, timeout=2))
        self.assertLess(time.time() - start, 30)
        with self.assertRaises(ChildProcessError):
            os.waitpid(pid, os.WNOHANG)


if __name__ == '__main__':
    unittest.main()
