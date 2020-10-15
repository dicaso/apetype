import os
import sys

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
