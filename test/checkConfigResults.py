# -*- coding: iso-8859-1 -*-
import json
import shutil
import shlex
import sys
import tempfile
import os
path = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(path[:path.rindex('/')])
from hotnet2 import hnap
from bin import findComponents as fc

def parse_args(raw_args): 
    description = "Runs the runHotNet2 script with the given config files and compares the\
                   resulting components to those in the given results files."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-c', '--config_files', nargs='+',
                        help='Paths to config files to pass to runHotNet2. Note that\
                              output_directory parameter will be ignored. In addition, these\
                              config files must specify exactly one value for delta.')
    parser.add_argument('-r', '--results_files', nargs='+',
                        help='Paths to results.json files whose components the run output should\
                              be compared to. Note that the components section is the only one\
                              that is considered.')
    
    return parser.parse_args(raw_args)

def run(args):
    if len(args.config_files) != len(args.results_files):
        raise ValueError("An equal number of config files and results files are required.")
    
    for i, config in enumerate(args.config_files):
        tmp_dir = tempfile.mkdtemp()
        
        # add output_directory argument if not included
        with open(config) as f:
            run_args = f.read()
        if not "-o" in run_args or "--output_directory" in run_args:
            run_args = "--output_directory {}\n".format(tmp_dir) + run_args
        
        # parse args and ensure output directory is set to temporary directory
        run_args = fc.get_parser().parse_args(shlex.split(run_args))
        run_args.output_directory = tmp_dir
        
        # run HotNet2
        fc.run(run_args)
        
        # load results and compare to expected
        with open(args.results_files[i]) as f:
            expected = json.load(f)
        with open("{}/delta_{}/results.json".format(tmp_dir, run_args.deltas[0])) as f:
            actual = json.load(f)
            
        shutil.rmtree(tmp_dir)
        
        if set(expected) != set(actual):
            raise AssertionError("FAIL. Config {} did not produce the expected results.".format(config))
        
        print "Config {}: OK".format(config)
        
    print "PASS. All config files produced the expected results."
        

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
