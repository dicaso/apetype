import os
import sys
import multiprocessing as mp

if sys.platform.lower() == "win32":
    # Command to activate term colors in windows
    os.system('')

# ANSI codes for different styles
class termstyle():
    "Info: https://en.wikipedia.org/wiki/ANSI_escape_code"
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    UNDERLINE = '\033[4m'
    BOLD = '\033[1m'
    RESET = '\033[0m'

# Counting semaphore
class CountSemaphore(object):
    def __init__(self, value, manager):
        self._value = value
        self._semaphore = manager.BoundedSemaphore(value)
        self._condition = manager.Condition()
        #self._counter = mp.Value('i', 0)

    def acquire(self, procnr):
        #with self._counter.get_lock():
        #    self._setvalue = self._counter.value
        #    self._counter.value += 1
        # Above strategy failed due to racing
        if procnr >= self._value: return False
        else: return self._semaphore.acquire(False)

    def release(self):
        return self._semaphore.release()
