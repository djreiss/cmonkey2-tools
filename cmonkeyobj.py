import os
import cStringIO
import sqlite3 as sql3
import gzip,bz2
import cPickle as pickle
import ConfigParser
import pandas as pd
import numpy as np
from numpy import nan as NA

from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import weblogolib as wl

import utils as ut

"""
from cmonkeyobj import cMonkey2 as cm2
b = cm2('eco-out-001/cmonkey_run.db')
pd.Series([b.get_cluster_info(k)['residual'] for k in range(1,b.k_clust)]).plot(kind='hist',bins=20)
pd.DataFrame([b.get_cluster_info(k)['pclusts'] for k in range(1,b.k_clust)]).plot(kind='hist',bins=20,stacked=True)
"""

## TBD: plotting motif locations relative to gene start. See
##  http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc211
##  https://www.biostars.org/p/96470/#97713

class cMonkey2:
    dbfile = '' #None
    tables = {} #None
    iteration = 2001
    k_clust = 999 #None
    organism = '' #'eco'
    species = '' #'Escherichia_coli_K12' #None
    taxon_id = None
    ratios = pd.DataFrame()
    config = ConfigParser.ConfigParser() #None
    stats = None

    def __init__( self, dbfile ):
        self.dbfile = dbfile
        conn = sql3.connect( dbfile )
        tmp = pd.read_sql('select max(iteration) from iteration_stats', conn) ##last_iteration from run_infos', conn)
        conn.close()
        self.iteration = tmp.max()[0] ## get iteration
        print 'iteration =', self.iteration
        self.tables = self.__read_all_tables( dbfile, iteration=self.iteration )
        #self.iteration = max(self.tables['motif_infos'].iteration)
        self.k_clust = self.tables['run_infos'].num_clusters[0] ##max(self.tables['row_members'].cluster)
        self.organism = self.tables['run_infos'].organism[0]
        self.species = self.tables['run_infos'].species[0]
        self.config = self.load_config()

    def __read_all_tables( self, dbfile, iteration=2000 ): #limit=None ):
        """read out all tables in the sql3 db file into a dict of pandas dataframes"""
        conn = sql3.connect( dbfile )
        tnames = pd.read_sql("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name", conn)
        tables = {}
        for tname in tnames.name.values:
            #print tname
            tmp = pd.read_sql( 'select * from %s limit 3' % tname, conn )
            if tname != 'motif_infos' and 'iteration' in tmp.columns.values.tolist():
                query = 'select * from %s where iteration=' + str(iteration)
            else:
                query = 'select * from %s'
            table = pd.read_sql(query % tname, conn)
            if tname == 'motif_infos':
                table = table[ table.iteration == iteration ]
            tables[ tname ] = table

        conn.close()
        table = tables[ 'meme_motif_sites' ]
        table = table.ix[ np.in1d( table.motif_info_id, tables[ 'motif_infos' ].index.values ) ]
        tables[ 'meme_motif_sites' ] = table
        return tables

    def reload( self ):
        conn = sql3.connect( self.dbfile )
        tmp = pd.read_sql('select max(iteration) from iteration_stats',conn)
        conn.close()
        self.iteration = tmp.max()[0] ## get iteration
        print 'iteration =', self.iteration
        self.tables = self.__read_all_tables( self.dbfile, iteration=self.iteration )        
        self.stats = None

    def get_feature_names( self ):
        feature_names_file = './cache/' + self.species + '_feature_names'
        feature_names = pd.read_table( feature_names_file, sep='\t', header=None, skiprows=4 )
        feature_names.columns = ['id','names','type']
        #feature_names = feature_names.set_index( 'names' )
        return feature_names

    def get_features( self ):
        features_file = './cache/' + self.species + '_features'
        features = pd.read_table( features_file, sep='\t', header=0, skiprows=16 )
        cols = features.columns.values; cols[0] = 'id'; features.columns = cols
        #features = features.set_index( 'od' )
        return features

    def get_genome_seqs( self ):
        features = self.get_features()
        contigs = np.unique( features.contig )
        seqs = {}
        for contig in contigs:
            genome_file = './cache/' + self.species + '_' + contig
            seq = ut.readLines( genome_file )[0].strip().upper()
            seqs[contig] = seq
        return seqs

    def get_networks( self, include_operons=True, include_string=True ):
        networks = {}
        taxid = self.load_taxon_id()
        if include_operons and os.path.exists('./cache/gnc' + str(taxid) + '.named'):
            op_file = './cache/gnc' + str(taxid) + '.named'
            operons = pd.read_table( op_file )
            networks['operons'] = operons
        if include_string and os.path.exists('./cache/' + str(taxid) + '.gz'):
            string_file = './cache/' + str(taxid) + '.gz'
            string = pd.read_table( gzip.GzipFile( string_file ), header=None )
            networks['string'] = string
        return networks

    ## see http://www.kegg.jp/kegg/rest/keggapi.html
    ## and http://biopython.org/DIST/docs/api/Bio.KEGG.REST-module.html
    ## another option: http://www.genome.jp/kegg-bin/show_organism?org=eco
    def load_taxon_id( self, in_code=None ):
        ''' lets try getting it directly from KEGG based on inputted organism 3-letter code
            a bit hairy but it works!  TODO: cache the org_table and gen_table in cache/'''
        if self.taxon_id is not None:
            return self.taxon_id
        import Bio.KEGG.REST as kegg ## requires BioPython 1.65 or later!
        if in_code is None:
            in_code = self.tables['run_infos'].organism[0]

        org_table = kegg.kegg_list('organism').readlines()
        org_table = ''.join( org_table )
        buf = cStringIO.StringIO( org_table )
        org_table = pd.read_table( buf, sep='\t', header=None )
        #full_org_name = org_table.ix[org_table[1]==in_code][2].values[0]
        buf.close()
        kegg_code = org_table.ix[org_table[1]==in_code][0].values[0]

        gen_table = kegg.kegg_list('genome').readlines()
        gen_table = ''.join( gen_table )
        buf = cStringIO.StringIO( gen_table )
        gen_table = pd.read_table( buf, sep='\t', header=None )
        buf.close()
        taxon_id = int(gen_table.ix[ gen_table[0] == 'genome:'+kegg_code ][1].values[0].split(', ')[2].split('; ')[0])
        self.taxon_id = taxon_id
        return taxon_id

    def load_ratios( self, ratios_file=None ):
        if ratios_file is None:
            ratios_file = os.path.dirname(self.dbfile) + '/ratios.tsv.gz'
        if self.ratios is None:
            self.ratios = pd.read_table( gzip.GzipFile( ratios_file ), sep='\t' )
        return self.ratios

    def load_config( self, config_file=None ):
        """then can do e.g., b.config.getfloat('Rows', 'scaling_constant')
           or simply, dict(b.config.items('Rows'))"""
        if config_file is None:
            config_file = os.path.dirname(self.dbfile) + '/final.ini'
        config_parser = ConfigParser.ConfigParser()
        config_parser.read( config_file )
        self.config = config_parser
        return self.config

    def pickle_all( self, outfile=None, include_genome=False, include_networks=False ):
        '''Try to pickle up ALL relevant info from the cmonkey run
           can load it via b = pickle.load(gzip.GzipFile(outfile)) '''
        ## another thing to try is to load the 
        feature_names = self.get_feature_names()
        features = self.get_features()
        genome = None
        if include_genome:
            genome = self.get_genome_seqs()
        networks = None
        if include_networks:
            networks = self.get_networks()
        self.load_ratios()
        self.load_config()
        self.get_stats()
        ## do pickling here
        if outfile is None:
            outfile = gzip.GzipFile( os.path.dirname(self.dbfile) + '/dump.pkl.gz', 'wb' )
        obj = { 'b': self, 
                'feature_names': feature_names,
                'features': features,
                'genome': genome,
                'networks': networks }
        print outfile
        pickle.dump( obj, outfile )
        outfile.close()

    def get_rows( self, k ):
        t1 = self.tables['row_members']
        t1 = t1[ t1.iteration == self.iteration ]
        t1 = t1[ t1.cluster == k ]
        t2 = self.tables['row_names']
        t2 = pd.merge( t1, t2, on='order_num' )
        return t2.name.values

    def get_cols( self, k ):
        t1 = self.tables['column_members']
        t1 = t1[ t1.iteration == self.iteration ]
        t1 = t1[ t1.cluster == k ]
        t2 = self.tables['column_names']
        t2 = pd.merge( t1, t2, on='order_num' )
        return t2.name.values

    def get_ratios( self, k=None, rows=None, cols=None, included=True ):
        """Extract submatrix of ratios for cluster or rows/cols. 
        If ~included, extract submatrix of ratios for conditions NOT in cluster."""
        if self.ratios is None:
            ratios = self.load_ratios()
        if k is not None:
            if rows is None:
                rows = self.get_rows( k )
            if cols is None:
                cols = self.get_cols( k )
        if not included:
            cols = ratios.columns.values[ np.in1d( ratios.columns.values, cols, invert=True ) ]
        rats = self.ratios.ix[ rows, cols ]
        return rats

    def plot_ratios( self, k=None, rows=None, cols=None, included=True, kind='line' ):
        ## see http://pandas.pydata.org/pandas-docs/version/0.15.0/visualization.html -- cool!
        ## can use kind = 'box' too!
        rats = self.get_ratios( k, rows, cols, included )
        rats = rats.transpose()

        if kind == 'box': ## sort by mean of columns
            means = rats.mean(1)
            tmp = pd.concat( [rats, means], 1 )
            cols = tmp.columns.values; cols[-1] = 'MEANS'; tmp.columns = cols
            tmp = tmp.sort( ['MEANS'] )
            tmp = tmp.drop( 'MEANS', 1 )
            rats = tmp.transpose()
            rats.plot(kind=kind, use_index=False, title='Cluster %d'%(k), legend=False, sym='.')
        else:
            rats.plot(kind=kind, use_index=False, title='Cluster %d'%(k), legend=False)
        ## use plt.close() to close the window

    def get_cluster_info( self, k ):
        t1 = self.tables['cluster_stats']
        t1 = t1[ t1.cluster == k ]
        #t1 = t1.drop( ['iteration', 'cluster'], 1 )

        t2 = self.tables['motif_infos']
        t2 = t2[ t2.cluster == k ]
        #t2 = t2.drop( ['iteration', 'cluster'], 1 )

        ## Extract it.
        out = {'residual':t1.residual.values[0],
               'nrows':t1.num_rows.values[0],
               'ncols':t1.num_cols.values[0],
               'e_values':t2.evalue.values}

        ## Also get p-clust
        pclusts = np.array([self.get_motif_pclust(k,i) for i in range(1,t2.shape[0]+1)])
        out['pclusts'] = pclusts
        
        return out

    def get_cluster_networks( self, k ):
        networks = self.get_networks()
        genes = self.get_rows( k )
        out_nets = {}
        if 'string' in networks.keys():
            string = networks['string']
            string = string.ix[ np.in1d(string[0], genes) ] ## slow!
            string = string.ix[ np.in1d(string[1], genes) ]
            out_nets['string'] = string
        if 'operons' in networks.keys():
            ops = networks['operons']
            ops = ops.ix[ np.in1d(ops.SysName1, genes) | np.in1d(ops.SysName2, genes) ]
            ops = ops.ix[ ops.bOp == True ]
            out_nets['operons'] = ops
        return out_nets

    ## see https://www.udacity.com/wiki/creating-network-graphs-with-python
    def plot_cluster_networks( self, k ):
        import networkx as nx
        out_nets = self.get_cluster_networks( k )
        gr = nx.Graph()
        if 'string' in out_nets.keys():
            strng = out_nets[ 'string' ]
            buf = cStringIO.StringIO()  ## round-about way to do it but wtf?
            strng.to_csv( buf, sep='\t', header=False, index=False )
            buf.flush(); buf.seek(0)
            gr = nx.read_weighted_edgelist( buf )
            buf.close()
        if 'operons' in out_nets.keys():
            ops = out_nets[ 'operons' ]
            ops = ops.ix[ ops.bOp == True ]
            ops = ops[ ['SysName1','SysName2','pOp'] ]
            ops.pOp = ops.pOp * 1000.
            buf = cStringIO.StringIO()  ## round-about way to do it but wtf?
            ops.to_csv( buf, sep='\t', header=False, index=False )
            buf.flush(); buf.seek(0)
            gr2 = nx.read_weighted_edgelist( buf )
            buf.close()
            #gr2 = nx.Graph( [ tuple(x) for x in ops[['SysName1','SysName2']].to_records(index=False) ], 
            #                weight=ops.pOp.values*1000, typ='operons' )
            ## from https://stackoverflow.com/questions/11758774/merging-two-network-maps-in-networkx-by-unique-labels :
            gr.add_nodes_from(gr2.nodes(data=True))
            gr.add_edges_from(gr2.edges(data=True)) #, weight=gr2.graph['weight'], type=gr2.graph['type'])

        pos = nx.spring_layout(gr, k=0.9, iterations=2000)
        ## requires installation of graphviz-dev and pygraphviz:
        ##from networkx import graphviz_layout
        ##pos = nx.graphviz_layout( gr, prog='neato'
        pos2 = { i:k for i,k in pos.items() if i in gr2.nodes() }
        nx.draw_networkx_edges(gr2, pos2, edge_color='r', width=4, alpha=0.5)
        nx.draw_networkx(gr, pos, node_size=50, node_color='b', edge_color='b', font_size=7, width=2, alpha=0.3)

    def clusters_w_genes( self, genes ):
        t1 = self.tables['row_members']
        t1 = t1[ (t1.iteration == self.iteration) ]
        t2 = self.tables['row_names']
        t2 = t2[ np.in1d(t2.name, genes) ]
        t2 = pd.merge( t1, t2, on='order_num' )
        t2 = t2.drop( ['iteration', 'order_num'], 1 )
        return t2

    def clusters_w_conds( self, conds ):
        t1 = self.tables['column_members']
        t1 = t1[ (t1.iteration == self.iteration) ]
        t2 = self.tables['column_names']
        t2 = t2[ np.in1d(t2.name, conds) ]
        t2 = pd.merge( t1, t2, on='order_num' )
        t2 = t2.drop( ['iteration', 'order_num'], 1 )
        return t2

    def cluster_summary( self ):
        tab = self.tables['cluster_stats']
        infos = { k: self.get_cluster_info(k+1) for k in range(self.k_clust) }
        tab[ 'e_value1' ] = pd.Series( [ infos[k]['e_values'][0] if 
                               len(infos[k]['e_values']) > 0 else NA for k in range(self.k_clust) ] )
        tab[ 'e_value2' ] = pd.Series( [ infos[k]['e_values'][1] if 
                               len(infos[k]['e_values']) > 1 else NA for k in range(self.k_clust) ] )
        tab[ 'p_clust1' ] = pd.Series( [ infos[k]['pclusts'][0] if 
                               len(infos[k]['pclusts']) > 0 else NA for k in range(self.k_clust) ] )
        tab[ 'p_clust2' ] = pd.Series( [ infos[k]['pclusts'][1] if 
                               len(infos[k]['pclusts']) > 1 else NA for k in range(self.k_clust) ] )
        tab = tab.set_index( tab.cluster )
        tab = tab.drop( ['iteration', 'cluster'], axis=1 )
        return tab

    def get_stats( self ):
        if self.stats is not None:
            return self.stats
        conn = sql3.connect( self.dbfile )
        table = pd.read_sql('select * from iteration_stats', conn)
        conn.close()
        tmp = self.tables['statstypes'].copy()
        tmp.index = tmp.index + 1
        table = pd.merge(table,tmp,left_on='statstype',right_index=True)
        tmp = table.groupby( 'name' )
        tmp = { name:df for name,df in tmp }
        for name in tmp.keys():
            tmp2 = tmp[name]
            tmp2.index = tmp2.iteration
            tmp2 = tmp2.drop( ['statstype', 'category', 'name', 'iteration'], axis=1 )
            tmp2.columns=[name]
            tmp[name] = tmp2
        #if 'SetEnrichment' in tmp.keys():
        #    pvs = pd.read_csv( os.path.dirname(self.dbfile) + '/setEnrichment_pvalue.csv', index_col=0 )
        #    pvs = pvs.fillna( 1.0 )
        #    tmp['SetEnrichment'] = np.log10(pvs+1e-30).median(1) ##.plot()
        table = pd.concat( tmp, axis=1 )
        table.columns = [i[0] for i in table.columns.values]
        self.stats = table
        return table

    def plot_stats( self ):
        table = self.get_stats()
        ut.setup_text_plots( usetex=False )
        if 'SetEnrichment' in table.columns.values:
            table.SetEnrichment.replace( 0, NA, inplace=True )
        table.plot( subplots=True, layout=[3,-1], sharex=True, legend=True, fontsize=8 )
        #fig, axes = plt.subplots(nrows=3, ncols=3, sharex=True)
        #for i, c in enumerate(table.columns):
        #    table[c].plot( ax=axes[i/3][i%3], title=c )

    def __get_motif_id(self, cluster_num, motif_num):
        motif_infos = self.tables['motif_infos']
        rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
                            (motif_infos.cluster==cluster_num) & 
                            (motif_infos.motif_num==motif_num)].index.values[0]+1
        return rowid
        #motif_id = self.tables['meme_motif_sites'].ix[rowid].motif_info_id
        #return motif_id

    def get_motif_pssm(self, cluster_num, motif_num):
        """export the specified motif to a pandas dataframe
        Parameters:
        - cluster_num: bicluster number
        - motif_num: motif number
        """
        #conn = sql3.connect(self.dbfile)
        #cursor = conn.cursor()
        #cursor.execute('select max(iteration) from motif_infos')
        #iteration = cursor.fetchone()[0]

        #query = 'select rowid from motif_infos where iteration=? and cluster=? and motif_num=?'
        #params = [self.iteration, cluster_num, motif_num]
        #cursor.execute(query, params)
        #rowid = cursor.fetchone()[0]
        #motif_infos = self.tables['motif_infos']
        #rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
        #                    (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        rowid = self.__get_motif_id(cluster_num, motif_num)

        #query = 'select a,c,g,t from motif_pssm_rows where iteration=? and motif_info_id=?'
        #params = [self.iteration, rowid]
        #pssm = pd.read_sql( query, conn, params=params )
        motif_pssm_rows = self.tables['motif_pssm_rows']
        pssm = motif_pssm_rows[(motif_pssm_rows.iteration==self.iteration) & (motif_pssm_rows.motif_info_id==rowid)]
        pssm.drop( ['motif_info_id', 'iteration', 'row'], 1, inplace=True )
        return pssm

    def get_motif_sites(self, cluster_num, motif_num=None):
        #motif_infos = self.tables['motif_infos']
        #rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
        #                    (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        rowid = self.__get_motif_id(cluster_num, motif_num)
        print rowid

        sites = self.tables['meme_motif_sites']
        sites = sites[ sites.motif_info_id == rowid ]
        sites = sites.drop( ['motif_info_id'], 1 )

        feature_names = self.get_feature_names()
        tmp = pd.merge( sites, feature_names, left_on='seq_name', right_on='id' )
        tmp = tmp[ np.in1d( tmp.names.values, self.tables['row_names'].name.values ) ]
        tmp = tmp.drop( ['seq_name', 'type'], 1 )
        tmp = tmp.drop_duplicates()

        return tmp ## need to update genes based on synonyms

    def plot_motif_sites(self, cluster_num, motif_num):
        """THIS NEEDS MORE WORK but has the beginnings of something...
        TODO: multiple motifs on same tracks, include ALL genes (i.e. in operons that were not included),
              do reverse-complement positioning correctly (based on gene strand), 
              use MAST scan output (from b.tables['motif_annotations'])
        """
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Graphics import GenomeDiagram
        from reportlab.lib.units import cm
        from reportlab.lib import colors

        """To get this to work: download http://www.reportlab.com/ftp/fonts/pfbfer.zip
           and unzip it into /usr/lib/python2.7/dist-packages/reportlab/fonts/
        """

        motif_sites = self.get_motif_sites(cluster_num, motif_num)
        pv_range = np.max(-np.log10(motif_sites.pvalue.values)) - 4 ## divide -log10(pval) by this to get alpha to use
        len_range = np.max(motif_sites.start.values) + 10

        gdd = GenomeDiagram.Diagram('Motif sites: %d, %d' % (cluster_num, motif_num))

        for i in range(motif_sites.shape[0]):
            gdt_features = gdd.new_track(1, start=0, end=len_range, greytrack=True, greytrack_labels=1,
                                         name=motif_sites.names.values[i], scale=True, greytrack_fontsize=4)
            gds_features = gdt_features.new_set()
            col = colors.red.clone()
            col.alpha = ( -np.log10(motif_sites.pvalue.values[i]) - 4 ) / pv_range
            m_start = motif_sites.start.values[i]
            m_len = len(motif_sites.seq.values[i])
            m_strand = motif_sites.reverse.values[i]
            if m_strand == 0:
                m_strand = -1
            feature = SeqFeature(FeatureLocation(m_start, m_start+m_len-1), strand=m_strand)
            gds_features.add_feature(feature, name=str(i+1), label=False, color=col)

        gdd.draw(format='linear', pagesize=(15*cm,motif_sites.shape[0]*cm/2), fragments=1, start=0, end=len_range+10)
        ##gdd.write("GD_labels_default.pdf", "pdf") ## looks like only output is to file, so do this:
        #output = cStringIO.StringIO()
        #gdd.write(output, 'png', dpi=300)
        #output.seek(0)
        output = gdd.write_to_string(output='png', dpi=300)
        output = cStringIO.StringIO(output)
        img = mpimg.imread(output)
        plt.axis('off')
        imgplot = plt.imshow( img, interpolation='bicubic' )
        output.close()
        return gdd

    def get_motif_pclust(self, cluster_num, motif_num):
        rowid = self.__get_motif_id(cluster_num, motif_num)
        sites = self.tables['meme_motif_sites']
        sites = sites[ sites.motif_info_id == rowid ]
        #sites = sites.drop( ['motif_info_id'], 1 )
        return np.mean( np.log10(sites.pvalue.values) )

    def get_biop_motif(self, cluster_num, motif_num, option='sites'):
        ##import egrin2.export_motifs as em
        """export the specified motif to a biopython motif object
        Parameters:
        - cluster_num: bicluster number
        - motif_num: motif number
        - option of how to translate - sites: jaspar 'sites' file; pfm: jaspar 'pfm' file
        """
        #conn = sql3.connect(self.dbfile)
        #cursor = conn.cursor()
        #cursor.execute('select max(iteration) from motif_infos')
        #iteration = cursor.fetchone()[0]

        #query = 'select rowid from motif_infos where iteration=? and cluster=? and motif_num=?'
        #params = [self.iteration, cluster_num, motif_num]
        #cursor.execute(query, params)
        #rowid = cursor.fetchone()[0]
        #motif_infos = self.tables['motif_infos']
        #rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
        #                    (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        rowid = self.__get_motif_id(cluster_num, motif_num)
        #mot_info = pd.read_sql('select * from motif_infos where rowid=?', conn, params=[rowid])
        #mot_sites = pd.read_sql('select * from meme_motif_sites where motif_info_id=?', conn, params=[rowid])
        mot_sites = self.tables['meme_motif_sites'][self.tables['meme_motif_sites'].motif_info_id == rowid]
            
        output = cStringIO.StringIO()
        ## ONE WAY TO TRY -- but Bio.motifs cant parse the incomplete MEME file
        ##output.write(em.MEME_FILE_HEADER % (0.25, 0.25, 0.25, 0.25))
        ##em.write_pssm(output, cursor, os.path.dirname(self.dbfile), cluster_num, rowid,
        ##              motif_num, mot_info['evalue'][0], 10)
        ##output.seek(0)
        ##mot = motifs.read( output, 'meme' )
            
        ## Second way - create a jaspar 'pfm' file from the pssm
        if option == 'pfm':
            #query = 'select a,c,g,t from motif_pssm_rows where iteration=? and motif_info_id=?'
            #params = [self.iteration, rowid]
            #pssm = pd.read_sql( query, conn, params=params )

            motif_pssm_rows = self.tables['motif_pssm_rows']
            pssm = motif_pssm_rows[(motif_pssm_rows.iteration==self.iteration) & (motif_pssm_rows.motif_info_id==rowid)]
            pssm = pssm.drop( ['motif_info_id', 'iteration', 'row'], 1 )

            counts = np.round( pssm * mot_sites.shape[0] ).transpose()
            counts.to_string(output, header=False, index=False )
            output.seek(0)
            mot = motifs.read( output, 'pfm' )

            ## Third way - create a jaspar 'sites' file
        elif option == 'sites':
            seqs = {}
            for i in mot_sites.index.values:
                name = mot_sites.ix[i].seq_name
                flank_left = mot_sites.ix[i].flank_left
                flank_left = Seq(flank_left if flank_left is not None else "", IUPAC.IUPACAmbiguousDNA()).lower()
                seq = Seq(mot_sites.ix[i].seq, IUPAC.IUPACAmbiguousDNA())
                flank_right = mot_sites.ix[i].flank_right
                flank_right = Seq(flank_right if flank_right is not None else "", IUPAC.IUPACAmbiguousDNA()).lower()
                full_seq = flank_left + seq + flank_right
                bs = SeqRecord( full_seq, id=name )
                seqs[i] = bs
                    
            SeqIO.write(seqs.values(), output, 'fasta')
            output.seek(0)
            mot = motifs.read( output, 'sites' )
            
        output.close()
        ## Note Bio.motifs.weblogo() uses the weblogo server (slow? requires connection.)
        #kwargs = dict(color_scheme='classic')
        #mot.weblogo('file.png', color_scheme='color_classic') ## note, can use format='PDF'
        #img = mpimg.imread('file.png')
        #imgplot = plt.imshow( img )
        #plt.show()
        return mot

    ## This uses weblogolib package to create files directly (installed as weblogo via pip)
    ## https://code.google.com/p/weblogo/
    def plot_motif( self, cluster_num, motif_num, img_format='png' ):
        #conn = sql3.connect(self.dbfile)
        #cursor = conn.cursor()
        #cursor.execute('select max(iteration) from motif_infos')
        #iteration = cursor.fetchone()[0]

        #query = 'select rowid from motif_infos where iteration=? and cluster=? and motif_num=?'
        #params = [self.iteration, cluster_num, motif_num]
        #cursor.execute(query, params)
        #rowid = cursor.fetchone()[0]
        #mot_info = pd.read_sql('select * from motif_infos where rowid=?', conn, params=[rowid])
        #mot_sites = pd.read_sql('select * from meme_motif_sites where motif_info_id=?', conn, params=[rowid])

        #motif_infos = self.tables['motif_infos']
        #rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
        #                    (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        rowid = self.__get_motif_id(cluster_num, motif_num)
        mot_sites = self.tables['meme_motif_sites'][self.tables['meme_motif_sites'].motif_info_id == rowid]

        ldata = wl.LogoData.from_seqs(wl.SeqList(mot_sites.seq.values.tolist(), wl.unambiguous_dna_alphabet))
        options = wl.LogoOptions()
        options.fineprint = os.path.dirname(self.dbfile) + ' %03d %03d' % ( cluster_num, motif_num )
        format = wl.LogoFormat(ldata, options) 
        format.color_scheme = wl.classic
        format.resolution = 150
        if img_format == 'png':
            tmp = wl.png_formatter( ldata, format )
            output = cStringIO.StringIO(tmp)
            img = mpimg.imread(output)
            plt.axis('off')
            imgplot = plt.imshow( img )
            #plt.show()
            return plt
        elif img_format == 'svg':
            tmp = wl.svg_formatter( ldata, format )
            return tmp
            ## note then can do e.g. ut.writeLines(svg.split('\n'),'test.svg')
            

