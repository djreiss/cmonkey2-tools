## Code for running cmonkey2 from ipython, with setting arguments internally
## before running ipython: 
##  setenv PYTHONPATH `pwd`/cmonkey2
## Assumes cmonkey2/ is in current working directory that ipython was started in

import sys

import cmonkey.cmonkey_run as cmr
import cmonkey.config as conf
import cmonkey.meme as meme

## capture logging output in ipython --
##   https://stackoverflow.com/questions/18786912/get-output-from-the-logging-module-in-ipython-notebook
##   use logging.DEBUG to see all output including debugging
import logging
reload(logging) ## is this really necessary? after testing - YES
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO, datefmt='%I:%M:%S')

conf.USER_INI_PATH = 'cmonkey2/config/default.ini'
meme.USER_TEST_FASTA_PATH = 'cmonkey2/config/fasta_test.fa'
cmr.USER_KEGG_FILE_PATH = 'cmonkey2/config/KEGG_taxonomy'
cmr.USER_GO_FILE_PATH = 'cmonkey2/config/proteome2taxid'
cmr.PIPELINE_USER_PATHS = {
    'default': 'cmonkey2/config/default_pipeline.json',
    'rows': 'cmonkey2/config/rows_pipeline.json',
    'rowsandmotifs': 'cmonkey2/config/rows_and_motifs_pipeline.json',
    'rowsandnetworks': 'cmonkey2/config/rows_and_networks_pipeline.json'
}

def arg_ext(parser):
    parser.set_defaults( organism='eco', out='eco-test', ##config='eco-out-001/final.ini', 
                         ratios='ratios.tsv.gz' ) ##, resume=True ) ##, interactive=True )
    ##parser.add_argument('args', nargs=argparse.REMAINDER)
    return parser

sys.argv = [sys.argv[0]] ## HACK to remove extra ipython args, like --pylab --pdb

## Mostly copied from cmonkey2/cmonkey.py:

args, params, ratios = conf.setup(arg_ext)
cm = cmr.CMonkeyRun(ratios, params)

cm.config_params['multiprocessing'] = True
cm.config_params['num_cores'] = 2

cm.prepare_run()

cm.run_iteration(1) ## run the 1st iteration only

cm.run_iterations(2, 10) ## run 10 iterations starting at 2

cm.run()   ## run the rest
cmonkey_run.cleanup() ## close the db

#### While cm.run() is running, can also do this to query the run:

from cmonkeyobj import cMonkey2 as cm2
b = cm2(cm['out_database']) ##'eco-test/cmonkey_run.db')

b.plot_motif(11,1)
b.plot_stats()

#### Can also try to spawn the cmviewer!!!

## THis doesn't work just yet -- but try the below instead...
sys.argv=[sys.argv[0]]
import os
import cmviewer.main as cmv
import cherrypy

cmv.outdir = os.path.join(os.getcwd(), 'eco-test')
conf = {'/': {'request.dispatch': cmv.setup_routes()},
        '/static': {'tools.staticdir.on': True,
                    'tools.staticdir.dir': os.path.join(cmv.current_dir, 'static')}}
cherrypy.config.update(conf)
app = cherrypy.tree.mount(None, config=conf)
cherrypy.server.socket_host = '0.0.0.0'
cherrypy.server.socket_port = cmv.args.port
cherrypy.quickstart(app)


## in ipython
%run -i cmonkey2/cmviewer/main.py --out eco-test
