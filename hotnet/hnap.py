import argparse
import shlex

class HotNetArgParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        if not arg_line.startswith('#'):
            for arg in shlex.split(arg_line):
                yield arg